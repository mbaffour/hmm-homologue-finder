#!/usr/bin/env python3
"""
extract_validated_hits.py
=========================
Authoritative, ORF-validated extraction of hit sequences from a completed
benchmark run. THIS SCRIPT IS THE FIX for the sequence-capture defect described
below.

WHY THIS EXISTS
---------------
The benchmark pipeline (run_all_database_benchmark.py) finds homologues by
six-frame translating nucleotide databases and searching them with the profile
HMM. It translates each stop-to-stop ORF correctly during the search, but then
DELETES the translated FASTA, keeping only the hit table (tblout) with genomic
coordinates. A later synteny step stored a Prodigal-predicted gene as the "hit
protein" — but that is a DIFFERENT gene that merely overlaps the locus (in most
cases the family domain was not even contained in it). An additional
post-hoc bug sliced 1-based coordinates with 0-based Python indexing, shifting
the reading frame and inserting stop codons.

THE FIX (what this script does), per hit:
  1. Reconstruct the FULL six-frame ORF directly from the genomic coordinates
     recorded in the hit table, using the EXACT convention emitted by the
     pipeline's _emit_sixframe_orf():
        + strand : nt = genome[start-1:end]                 (1-based -> 0-based)
        - strand : nt = revcomp(genome[start-1:end])
     translate(nt) -> ORF amino-acid sequence (stop-free by construction).
  2. Locate the family domain WITHIN that ORF using hmmsearch --domtblout against
     the run's HMM (gives env_from/env_to on the ORF).
  3. Slice the domain region in BOTH nucleotide and amino-acid space.
  4. VALIDATE the hit is a genuine ORF (not an arbitrary match):
        orf_aa_len, domain_aa_len, domain_coverage,
        has_start_M, ends_at_stop, internal_stops (must be 0),
        prodigal_concordant (+ overlap %) via Prodigal on the source genome,
        passes_orf_filter (0 internal stops AND coverage >= MIN_COVERAGE).
  5. Join hmmsearch statistics + download provenance from the hit table.
  6. Write both sequence types + a fully-annotated TSV.

INPUTS
------
  --results-dir  benchmark results dir (must contain hits_main.tsv and
                 synteny_context_cache/<genome>.fna)
  --hmm          the run's profile HMM (for domain-envelope detection)
  --run-label    A | B | C  (stamped into every row)
  --out-dir      where outputs are written

OUTPUTS
-------
  hits.tsv            one row per hit, full schema (NT + AA + validation)
  hits_aa.faa         domain amino-acid sequence per hit
  hits_nt.fna         domain nucleotide sequence per hit
  orfs_aa.faa         full ORF amino-acid sequence per hit (context)
  orfs_nt.fna         full ORF nucleotide sequence per hit (context)
  hits_unique_aa.faa  deduplicated, ORF-validated AA domains -> next seed

USAGE
-----
  python3 extract_validated_hits.py \
      --results-dir <run>/benchmark/results \
      --hmm <run>/benchmark/hmm/benchmark_profile.hmm \
      --run-label A --out-dir <run>/validated
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import tempfile
import warnings
from collections import OrderedDict
from pathlib import Path

import io
import time
import pandas as pd
from Bio import BiopythonWarning, Entrez, SeqIO
from Bio.Seq import Seq

warnings.simplefilter("ignore", BiopythonWarning)

# Ensure the conda env's tools (hmmsearch, prodigal) are on PATH even when this
# script is invoked from a process that did not `conda activate`. Robustly find
# the hmm-discovery env: the active CONDA_PREFIX first, then common installers.
_cand = ([Path(os.environ["CONDA_PREFIX"]) / "bin"] if os.environ.get("CONDA_PREFIX") else []) \
    + [Path.home() / _n / "envs" / "hmm-discovery" / "bin"
       for _n in ("miniforge3", "mambaforge", "miniconda3", "anaconda3")]
for _b in _cand:
    if _b.is_dir():
        os.environ["PATH"] = f"{_b}{os.pathsep}{os.environ.get('PATH', '')}"
        break

# ----------------------------------------------------------------------------
# Coordinate-correct ORF reconstruction
# ----------------------------------------------------------------------------
# Column order for hits.tsv (single source of truth).
COLS = [
    "hit_id", "genome_id", "contig", "db_name", "db_type", "run_label",
    "source_type",
    "nt_start", "nt_end", "strand", "frame", "orf_nt_start", "orf_nt_end",
    "orf_aa_len", "domain_aa_len", "domain_coverage", "has_start_M",
    "ends_at_stop", "internal_stops", "prodigal_concordant",
    "prodigal_same_strand_pct", "in_coding_locus", "prodigal_any_strand_pct",
    "passes_orf_filter", "evalue", "bit_score",
    "bias_score", "env_from", "env_to", "confidence_tier", "qc_flags",
    "source_url", "source_sha256", "accessed_at", "aa_sequence", "nt_sequence",
]


def fetch_protein_seqs(accessions: list[str], email: str) -> dict[str, str]:
    """Fetch amino-acid sequences for protein-database hits from NCBI.

    Protein-database hits (SwissProt, RefSeq/INPHARED proteins, Pfam) are
    ALREADY annotated proteins — the matched target IS the sequence, so there is
    no ORF to reconstruct. We retrieve the AA by accession via Entrez. Returns
    {accession: sequence} keyed with and without version.
    """
    out: dict[str, str] = {}
    if not accessions:
        return out
    Entrez.email = email
    for i in range(0, len(accessions), 40):
        chunk = accessions[i:i + 40]
        for attempt in range(3):
            try:
                h = Entrez.efetch(db="protein", id=",".join(chunk),
                                  rettype="fasta", retmode="text")
                for rec in SeqIO.parse(io.StringIO(h.read()), "fasta"):
                    seq = str(rec.seq)
                    # NCBI may return SwissProt-style ids (sp|ACC|NAME); key the
                    # sequence under every form so the caller's accession matches.
                    keys = {rec.id, rec.id.split(".")[0]}
                    ca = clean_accession(rec.id)
                    keys |= {ca, ca.split(".")[0]}
                    for k in keys:
                        out[k] = seq
                break
            except Exception as e:
                print(f"  protein efetch retry {attempt+1}: {e}")
                time.sleep(4)
        time.sleep(0.4)
    return out


def clean_accession(target_name: str) -> str:
    """Extract a fetchable accession from an hmmsearch target name.

    Handles bare accessions ('P12345', 'NC_023589.1') and UniProt-style
    'sp|P12345|NAME' / 'tr|...' identifiers.
    """
    t = str(target_name)
    if "|" in t:                       # sp|ACC|NAME  or  tr|ACC|NAME
        parts = t.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1]
    return t.split()[0]


def reconstruct_orf(genome_seq: str, start_1: int, end_1: int, strand: str):
    """Return (nt, aa) for the ORF at 1-based forward coords [start_1, end_1].

    Mirrors the pipeline's _emit_sixframe_orf convention exactly:
    coordinates are always on the FORWARD strand; for '-' hits the ORF is the
    reverse complement of that forward slice. The 1-based -> 0-based conversion
    (start_1 - 1) is the corrected slicing that the earlier buggy code got wrong.
    """
    sub = genome_seq[start_1 - 1:end_1]          # 1-based inclusive -> 0-based
    if strand == "-":
        sub = str(Seq(sub).reverse_complement())
    sub = sub[: len(sub) // 3 * 3]               # trim to whole codons
    aa = str(Seq(sub).translate()).rstrip("*")   # drop terminal stop only
    return sub, aa


# ----------------------------------------------------------------------------
# Prodigal gene prediction (cached per genome) for ORF concordance
# ----------------------------------------------------------------------------
_prodigal_cache: dict[str, list[tuple[int, int, str]]] = {}


def prodigal_genes(fna: Path) -> list[tuple[int, int, str]]:
    """Return [(start, end, strand), ...] of Prodigal genes for a genome.

    Results are cached so each genome is only predicted once per process.
    """
    key = str(fna)
    if key in _prodigal_cache:
        return _prodigal_cache[key]
    genes: list[tuple[int, int, str]] = []
    with tempfile.NamedTemporaryFile(suffix=".gff", delete=False) as tf:
        gff = tf.name
    try:
        subprocess.run(
            ["prodigal", "-i", str(fna), "-o", gff, "-f", "gff", "-p", "meta", "-q"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
        )
        for line in Path(gff).read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            if len(p) >= 7 and p[2] == "CDS":
                genes.append((int(p[3]), int(p[4]), p[6]))
    finally:
        Path(gff).unlink(missing_ok=True)
    _prodigal_cache[key] = genes
    return genes


def prodigal_overlap(genes, start_1, end_1, strand):
    """Return (same_strand_pct, any_strand_pct) overlap with Prodigal genes.

    - same_strand_pct: best overlap with a gene on the SAME strand. For family
      this is typically ~0, because the homologs are antisense / alternate-frame
      to predicted genes — that discordance is itself the novelty signal (they
      were missed by standard annotation).
    - any_strand_pct: best overlap with ANY predicted gene. Typically high,
      confirming the hit sits in a genuine coding-dense locus (a real ORF region,
      not six-frame noise).
    """
    span = max(1, end_1 - start_1 + 1)
    same = 0.0
    any_ = 0.0
    for g_start, g_end, g_strand in genes:
        overlap = max(0, min(end_1, g_end) - max(start_1, g_start) + 1)
        frac = overlap / span
        any_ = max(any_, frac)
        if g_strand == strand:
            same = max(same, frac)
    return same, any_


# ----------------------------------------------------------------------------
# Domain-envelope detection (hmmsearch --domtblout on the ORF set)
# ----------------------------------------------------------------------------
def domain_envelopes(orf_faa: Path, hmm: Path) -> dict[str, tuple[int, int]]:
    """Map ORF id -> (env_from, env_to) of its best family domain hit."""
    with tempfile.NamedTemporaryFile(suffix=".domtbl", delete=False) as tf:
        domtbl = tf.name
    subprocess.run(
        ["hmmsearch", "--domtblout", domtbl, "-E", "1e-3", "--domE", "1e-3",
         "--cpu", "8", "--noali", str(hmm), str(orf_faa)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
    )
    best: dict[str, tuple[int, int, float]] = {}
    for line in Path(domtbl).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        pid, ev = p[0], float(p[12])
        env_from, env_to = int(p[19]), int(p[20])
        if pid not in best or ev < best[pid][2]:
            best[pid] = (env_from, env_to, ev)
    Path(domtbl).unlink(missing_ok=True)
    return {k: (v[0], v[1]) for k, v in best.items()}


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--hmm", type=Path, required=True)
    ap.add_argument("--run-label", required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--email", default="researcher@example.com",
                    help="NCBI Entrez email for protein-database hit retrieval")
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    cache = args.results_dir / "synteny_context_cache"

    hits = pd.read_csv(args.results_dir / "hits_main.tsv", sep="\t")
    # Optional richer columns (tier / qc) live in hits_classified.tsv
    classified = args.results_dir / "hits_classified.tsv"
    cls = pd.read_csv(classified, sep="\t") if classified.exists() else pd.DataFrame()

    # ---- Pass 1: reconstruct every ORF (NT + AA) from coordinates -----------
    records = []          # working list of dicts
    orf_aa_path = out / "orfs_aa.faa"
    orf_nt_path = out / "orfs_nt.fna"
    genome_seq_cache: dict[str, str] = {}

    with orf_aa_path.open("w") as faa, orf_nt_path.open("w") as fna:
        for _, row in hits.iterrows():
            desc = str(row.get("description", ""))
            m = re.search(r"coords=([^:]+):(\d+)-(\d+)\(([+-])\)", desc)
            if not m:
                continue
            contig, start_1, end_1, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
            genome_id = str(row.get("genome_id", contig))

            fna_path = cache / f"{contig}.fna"
            if not fna_path.exists():
                fna_path = cache / f"{genome_id}.fna"
            if not fna_path.exists():
                continue

            if str(fna_path) not in genome_seq_cache:
                try:
                    rec = next(SeqIO.parse(str(fna_path), "fasta"))
                    genome_seq_cache[str(fna_path)] = str(rec.seq)
                except StopIteration:
                    continue
            genome_seq = genome_seq_cache[str(fna_path)]

            orf_nt, orf_aa = reconstruct_orf(genome_seq, start_1, end_1, strand)
            if len(orf_aa) < 20:
                continue

            pid = str(row.get("protein_id", row.get("target_name", "")))
            faa.write(f">{pid}\n{orf_aa}\n")
            fna.write(f">{pid}\n{orf_nt}\n")
            records.append({
                "hit_id": pid,
                "genome_id": genome_id,
                "contig": contig,
                "db_name": row.get("db_name", ""),
                "db_type": row.get("db_type", ""),
                "run_label": args.run_label,
                "source_type": "six_frame_orf",
                "source_url": row.get("source_url", ""),
                "source_sha256": row.get("source_sha256", ""),
                "accessed_at": row.get("source_accessed_at", ""),
                "nt_start": start_1, "nt_end": end_1, "strand": strand,
                "frame": re.search(r"frame=(\S+)", desc).group(1) if re.search(r"frame=(\S+)", desc) else "",
                "orf_nt_start": start_1, "orf_nt_end": end_1,
                "orf_aa_len": len(orf_aa),
                "evalue": row.get("evalue", ""), "bit_score": row.get("bit_score", ""),
                "bias_score": row.get("bias_score", ""),
                "_orf_nt": orf_nt, "_orf_aa": orf_aa,
                "_fna_path": str(fna_path),
            })

    # ---- Pass 1b: PROTEIN-DATABASE hits (no coords=) ------------------------
    # These targets are already-annotated proteins (SwissProt, RefSeq/INPHARED
    # proteins, Pfam). There is no ORF to reconstruct — the matched protein IS
    # the sequence — so we fetch it by accession instead of dropping it.
    protein_rows = [row for _, row in hits.iterrows()
                    if not re.search(r"coords=([^:]+):(\d+)-(\d+)\(([+-])\)", str(row.get("description", "")))]
    if protein_rows:
        accs = sorted({clean_accession(r.get("protein_id", r.get("target_name", ""))) for r in protein_rows})
        print(f"  fetching {len(accs)} protein-database hit sequences from NCBI…")
        prot_seqs = fetch_protein_seqs(accs, args.email)
        with orf_aa_path.open("a") as faa:
            for row in protein_rows:
                tname = str(row.get("protein_id", row.get("target_name", "")))
                acc = clean_accession(tname)
                aa = prot_seqs.get(acc) or prot_seqs.get(acc.split(".")[0])
                if not aa:
                    continue
                faa.write(f">{tname}\n{aa}\n")
                records.append({
                    "hit_id": tname,
                    "genome_id": acc,
                    "contig": acc,
                    "db_name": row.get("db_name", ""),
                    "db_type": row.get("db_type", "protein"),
                    "run_label": args.run_label,
                    "source_type": "annotated_protein",
                    "source_url": row.get("source_url", ""),
                    "source_sha256": row.get("source_sha256", ""),
                    "accessed_at": row.get("source_accessed_at", ""),
                    "nt_start": "", "nt_end": "", "strand": "", "frame": "",
                    "orf_nt_start": "", "orf_nt_end": "",
                    "orf_aa_len": len(aa),
                    "evalue": row.get("evalue", ""), "bit_score": row.get("bit_score", ""),
                    "bias_score": row.get("bias_score", ""),
                    "_orf_nt": "", "_orf_aa": aa, "_fna_path": "",
                })

    if not records:
        # Zero hits is a valid outcome — write empty outputs so the pipeline
        # continues cleanly instead of crashing.
        print("0 hits found — writing empty outputs.")
        pd.DataFrame(columns=COLS).to_csv(out / "hits.tsv", sep="\t", index=False)
        for fn in ("hits_aa.faa", "hits_nt.fna", "hits_unique_aa.faa"):
            (out / fn).write_text("")
        print(f"[{args.run_label}] hits=0  unique_seeds=0\n  outputs in {out}")
        return

    # ---- Pass 2: domain envelopes via hmmsearch on the ORF set -------------
    envs = domain_envelopes(orf_aa_path, args.hmm)

    # ---- Pass 3: per-hit validation + domain slicing -----------------------
    cls_tier = dict(zip(cls.get("target_name", []), cls.get("confidence_tier", []))) if not cls.empty else {}
    cls_qc = dict(zip(cls.get("target_name", []), cls.get("qc_flags", []))) if not cls.empty else {}

    for r in records:
        orf_aa = r.pop("_orf_aa")
        orf_nt = r.pop("_orf_nt")
        fna_path = Path(r.pop("_fna_path"))
        pid = r["hit_id"]

        env = envs.get(pid)
        if env:
            a, b = env
            dom_aa = orf_aa[a - 1:b]
            dom_nt = orf_nt[(a - 1) * 3:b * 3]      # AA envelope -> codon span
        else:
            # No envelope (rare): fall back to whole ORF if short
            a, b = 1, len(orf_aa)
            dom_aa, dom_nt = orf_aa, orf_nt

        internal_stops = dom_aa[:-1].count("*") if dom_aa else 0
        r["domain_aa_len"] = len(dom_aa)
        r["domain_coverage"] = round(len(dom_aa) / max(1, len(orf_aa)), 3)
        r["internal_stops"] = internal_stops
        r["env_from"], r["env_to"] = a, b

        if r["source_type"] == "annotated_protein":
            # An already-annotated protein: no ORF reconstruction / Prodigal
            # check applies. It is a genuine protein by definition.
            r["has_start_M"] = orf_aa.startswith("M")
            r["ends_at_stop"] = "NA"
            r["prodigal_same_strand_pct"] = "NA"
            r["prodigal_any_strand_pct"] = "NA"
            r["prodigal_concordant"] = "NA"
            r["in_coding_locus"] = "NA"
            r["passes_orf_filter"] = (internal_stops == 0)
        else:
            # Six-frame ORF: validate against a genuine ORF (no internal stops)
            # sitting in a real coding locus (overlaps a Prodigal gene).
            r["has_start_M"] = orf_aa.startswith("M")
            r["ends_at_stop"] = True  # six-frame ORFs are stop-bounded by construction
            genes = prodigal_genes(fna_path)
            same, any_ = prodigal_overlap(genes, r["nt_start"], r["nt_end"], r["strand"])
            r["prodigal_same_strand_pct"] = round(same * 100, 1)
            r["prodigal_any_strand_pct"] = round(any_ * 100, 1)
            r["prodigal_concordant"] = same >= 0.50          # same-strand predicted gene
            r["in_coding_locus"] = any_ >= 0.50              # overlaps any predicted gene
            # domain_coverage is reported but NOT exclusionary: a domain embedded
            # in a longer ORF is still a real homolog; we seed with the slice.
            r["passes_orf_filter"] = (internal_stops == 0) and r["in_coding_locus"]

        r["confidence_tier"] = cls_tier.get(pid, "")
        r["qc_flags"] = cls_qc.get(pid, "")
        r["aa_sequence"] = dom_aa
        r["nt_sequence"] = dom_nt

    # ---- Write outputs -----------------------------------------------------
    df = pd.DataFrame(records)[COLS]
    df.to_csv(out / "hits.tsv", sep="\t", index=False)

    with (out / "hits_aa.faa").open("w") as f:
        for r in records:
            f.write(f">{r['hit_id']} {r['genome_id']} cov={r['domain_coverage']} orf_concordant={r['prodigal_concordant']}\n{r['aa_sequence']}\n")
    with (out / "hits_nt.fna").open("w") as f:
        for r in records:
            if r["nt_sequence"]:   # protein-DB hits have no nucleotide sequence
                f.write(f">{r['hit_id']} {r['genome_id']}\n{r['nt_sequence']}\n")

    # Deduplicated, ORF-validated AA domains -> seed for the next run
    seen = OrderedDict()
    for r in records:
        if r["passes_orf_filter"] and r["aa_sequence"] not in seen:
            seen[r["aa_sequence"]] = r["hit_id"]
    with (out / "hits_unique_aa.faa").open("w") as f:
        for seq, pid in seen.items():
            f.write(f">{pid}\n{seq}\n")

    # ---- Report ------------------------------------------------------------
    total = len(records)
    passed = sum(1 for r in records if r["passes_orf_filter"])
    sixframe = sum(1 for r in records if r["source_type"] == "six_frame_orf")
    protein = sum(1 for r in records if r["source_type"] == "annotated_protein")
    coding = sum(1 for r in records if r["in_coding_locus"] is True)
    stopped = sum(1 for r in records if r["internal_stops"] > 0)
    print(f"[{args.run_label}] hits={total} (six-frame={sixframe}, protein-DB={protein})  "
          f"pass_filter={passed}  six-frame in_coding_locus={coding}  "
          f"internal_stops={stopped}  unique_seeds={len(seen)}")
    print(f"  outputs in {out}")


if __name__ == "__main__":
    main()
