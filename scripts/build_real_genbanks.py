#!/usr/bin/env python3
"""
build_real_genbanks.py
=====================
Rebuild family gene-neighbourhood GenBank files containing REAL nucleotide
sequence (not the N-placeholders used for clinker), so they open directly in
Artemis / Geneious / UGENE and can be viewed together with the GFF3.

For each hit genome:
  1. Obtain the genome sequence:
       - NCBI-accession genomes  -> Entrez efetch (nuccore)
       - metagenomic genomes     -> stream the source catalogue (GVD-AVrC / GPD)
  2. Gene-call with Prodigal for flanking-gene context.
  3. Cut a neighbourhood window (hit ORF + 5 flanking genes each side).
  4. Write a GenBank record with the REAL window sequence and CDS features:
       central CDS = the genuine family ORF (validated translation);
       flanking CDS = Prodigal genes.
  A genome containing >1 hit gets all its family ORFs annotated in one record.

INPUT
-----
  --hits-tsv : a run's hits.tsv (has genome_id, coords, strand, aa_sequence)
  --out-dir  : where the *.gbk files are written
  --email    : NCBI Entrez email (any valid address)

USAGE
-----
  python3 build_real_genbanks.py \
      --hits-tsv .../runA/hits.tsv --out-dir .../genbank_files_with_sequence
"""
from __future__ import annotations

import argparse
import gzip
import io
import os
import re
import subprocess
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

_cand = ([Path(os.environ["CONDA_PREFIX"]) / "bin"] if os.environ.get("CONDA_PREFIX") else []) \
    + [Path.home() / _n / "envs" / "hmm-discovery" / "bin"
       for _n in ("miniforge3", "mambaforge", "miniconda3", "anaconda3")]
for _b in _cand:
    if _b.is_dir():
        os.environ["PATH"] = f"{_b}{os.pathsep}{os.environ.get('PATH', '')}"
        break

FLANKS = 5
# Source catalogues for metagenomic genomes (id prefix -> URL)
CATALOGUES = {
    "GutCatV1_": "https://zenodo.org/records/11426065/files/AVrC_allrepresentatives.fasta.gz",  # GVD-AVrC
    "uvig_": "https://zenodo.org/records/6503062/files/GPD_sequences.fa.gz",                     # GPD
}


def is_ncbi(genome_id: str) -> bool:
    return bool(re.match(r"^[A-Z]{1,2}_?\d{5,8}", genome_id))


# ---------------------------------------------------------------------------
def _clean_name(title: str) -> str:
    """Turn an NCBI title into a phage name, e.g.
    'Escherichia phage vB_EcoP_G7C, complete genome' -> 'Escherichia phage vB_EcoP_G7C'."""
    t = re.split(r",\s*(complete|partial|whole|genome assembly|DNA)\b", title, maxsplit=1)[0]
    return t.strip().rstrip(",").strip()


def fetch_ncbi(ids: list[str], email: str) -> tuple[dict[str, str], dict[str, str]]:
    """Batch-fetch genome sequences AND phage names from NCBI nuccore.

    Returns (sequences, names): both keyed by accession (with and without
    version). `names` holds the organism/phage name parsed from the record
    title (empty for records without a usable title).
    """
    Entrez.email = email
    seqs: dict[str, str] = {}
    names: dict[str, str] = {}
    batch = 40
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        # --- phage names (esummary: fast, no sequence) ---
        for attempt in range(3):
            try:
                h = Entrez.esummary(db="nuccore", id=",".join(chunk))
                for rec in Entrez.read(h):
                    acc = rec.get("AccessionVersion", "")
                    nm = _clean_name(rec.get("Title", ""))
                    if acc and nm:
                        names[acc] = nm
                        names[acc.split(".")[0]] = nm
                break
            except Exception as e:
                print(f"  esummary retry {attempt+1}: {e}")
                time.sleep(5)
        # --- sequences (efetch fasta) ---
        for attempt in range(3):
            try:
                h = Entrez.efetch(db="nuccore", id=",".join(chunk),
                                  rettype="fasta", retmode="text")
                for rec in SeqIO.parse(io.StringIO(h.read()), "fasta"):
                    seqs[rec.id.split(".")[0]] = str(rec.seq)
                    seqs[rec.id] = str(rec.seq)
                break
            except Exception as e:
                print(f"  efetch retry {attempt+1} ({chunk[0]}…): {e}")
                time.sleep(5)
        time.sleep(0.4)  # be polite to NCBI
        print(f"  NCBI fetched {min(i+batch, len(ids))}/{len(ids)}")
    return seqs, names


def fetch_catalogue(url: str, wanted: set[str]) -> dict[str, str]:
    """Stream a gzipped catalogue and pull just the wanted contigs."""
    found: dict[str, str] = {}
    print(f"  streaming {url.split('/')[-1]} for {len(wanted)} contigs…")
    proc = subprocess.Popen(["bash", "-lc", f"curl -sS -L --retry 10 {url!r} | gunzip -c"],
                            stdout=subprocess.PIPE)
    cur_id, cur_seq, keep = None, [], False
    for raw in io.TextIOWrapper(proc.stdout, encoding="utf-8", errors="replace"):
        if raw.startswith(">"):
            if cur_id and keep:
                found[cur_id] = "".join(cur_seq)
                if len(found) >= len(wanted):
                    break
            cur_id = raw[1:].split()[0]
            keep = cur_id in wanted
            cur_seq = []
        elif keep:
            cur_seq.append(raw.strip())
    if cur_id and keep and cur_id not in found:
        found[cur_id] = "".join(cur_seq)
    proc.terminate()
    print(f"    recovered {len(found)}/{len(wanted)}")
    return found


# ---------------------------------------------------------------------------
_pg: dict[str, list] = {}


def prodigal_genes(genome_id: str, seq: str):
    if genome_id in _pg:
        return _pg[genome_id]
    with tempfile.NamedTemporaryFile("w", suffix=".fna", delete=False) as tf:
        tf.write(f">{genome_id}\n{seq}\n")
        fna = tf.name
    faa, gff = fna + ".faa", fna + ".gff"
    subprocess.run(["prodigal", "-i", fna, "-a", faa, "-o", gff, "-f", "gff", "-p", "meta", "-q"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    prots = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(faa, "fasta")} if Path(faa).exists() else {}
    genes = []
    if Path(gff).exists():
        idx = 0
        for line in Path(gff).read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            if len(p) >= 7 and p[2] == "CDS":
                idx += 1
                genes.append((int(p[3]), int(p[4]), 1 if p[6] == "+" else -1,
                              prots.get(f"{p[0]}_{idx}", "")))
    for f in (fna, faa, gff):
        Path(f).unlink(missing_ok=True)
    genes.sort()
    _pg[genome_id] = genes
    return genes


def build(genome_id: str, seq: str, hits: pd.DataFrame, out_dir: Path,
          phage_name: str = "") -> Path | None:
    genes = prodigal_genes(genome_id, seq)
    if not genes:
        return None
    # window = span of (all hit ORFs in this genome) + nearest flanking genes
    h_lo = int(hits.nt_start.min())
    h_hi = int(hits.nt_end.max())
    centre = (h_lo + h_hi) // 2
    nearby = sorted(sorted(genes, key=lambda g: abs((g[0]+g[1])//2 - centre))[:FLANKS*2+2],
                    key=lambda g: g[0])
    lo = max(1, min([h_lo] + [g[0] for g in nearby]) - 100)
    hi = min(len(seq), max([h_hi] + [g[1] for g in nearby]) + 100)
    window = seq[lo-1:hi]

    # Prefer the phage name in the human-readable fields; keep the accession in
    # the identifier. Metagenomic genomes have no name -> fall back to the id.
    label = phage_name or genome_id
    rec = SeqRecord(Seq(window), id=genome_id[:16], name=genome_id[:16],
                    description=f"{label} ({genome_id}) family neighbourhood (real sequence)",
                    annotations={"molecule_type": "DNA", "topology": "linear",
                                 "organism": label, "source": label})
    feats = []
    for (gs, ge, gst, aa) in nearby:
        feats.append(SeqFeature(FeatureLocation(max(0, gs-lo), max(1, ge-lo), strand=gst),
                                type="CDS", qualifiers={"product": ["flanking CDS"],
                                                        "translation": [aa or "X"]}))
    for _, hrow in hits.iterrows():
        hs, he = int(hrow.nt_start), int(hrow.nt_end)
        strand = 1 if hrow.strand == "+" else -1
        feats.append(SeqFeature(FeatureLocation(max(0, hs-lo), max(1, he-lo), strand=strand),
                                type="CDS", qualifiers={"gene": ["family"],
                                "product": ["homologue (HMM hit)"],
                                "translation": [str(hrow.aa_sequence).rstrip("*") or "X"]}))
    rec.features = sorted(feats, key=lambda f: int(f.location.start))
    # Filename: "<PhageName>_<accession>.gbk" when a name is known, else the id.
    stem = f"{label}_{genome_id}" if phage_name else genome_id
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", stem)[:70].strip("_")
    out = out_dir / f"{safe}.gbk"
    SeqIO.write(rec, str(out), "genbank")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits-tsv", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--email", default="researcher@example.com")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.hits_tsv, sep="\t")
    by_genome = {gid: sub for gid, sub in df.groupby("genome_id")}
    genome_ids = list(by_genome)
    print(f"{len(genome_ids)} genomes ({len(df)} hits)")

    # 1. NCBI genomes (sequences + phage names)
    ncbi_ids = [g for g in genome_ids if is_ncbi(g)]
    print(f"Fetching {len(ncbi_ids)} NCBI genomes…")
    seqs, names = fetch_ncbi(ncbi_ids, args.email)

    # 2. metagenomic genomes, grouped by catalogue prefix (uncultured -> no name)
    meta_ids = [g for g in genome_ids if not is_ncbi(g)]
    for prefix, url in CATALOGUES.items():
        want = {g for g in meta_ids if g.startswith(prefix)}
        if want:
            seqs.update(fetch_catalogue(url, want))

    # 3. build GenBanks (named by phage where known)
    built, named, missing = 0, 0, []
    for gid in genome_ids:
        seq = seqs.get(gid) or seqs.get(gid.split(".")[0])
        if not seq:
            missing.append(gid)
            continue
        nm = names.get(gid) or names.get(gid.split(".")[0], "")
        if build(gid, seq, by_genome[gid], args.out_dir, phage_name=nm):
            built += 1
            named += bool(nm)
    print(f"\nBuilt {built} real-sequence GenBank files in {args.out_dir} "
          f"({named} with phage names; {built-named} uncultured/metagenomic)")
    if missing:
        print(f"Could not retrieve {len(missing)} genomes: {missing}")


if __name__ == "__main__":
    main()
