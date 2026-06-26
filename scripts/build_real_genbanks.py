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
import socket
import subprocess
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO

socket.setdefaulttimeout(60)  # bound NCBI calls so unattended runs can't hang
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from env_paths import ensure_env_on_path  # noqa: E402  (sibling helper in scripts/)
ensure_env_on_path()

FLANKS = 7  # ORFs each side of the family gene (publication default for phage neighbourhoods)
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


def fetch_ncbi(ids: list[str], email: str) -> tuple[dict[str, str], dict[str, str], dict[str, list]]:
    """Batch-fetch genome sequences, phage names, AND the genome's REAL gene
    annotations from NCBI nuccore — by retrieving the full GenBank record
    (``rettype='gbwithparts'``) rather than bare FASTA. This is what lets the
    neighbourhood maps label flanking genes with their genuine ``/product`` names
    (terminase, major capsid, tail fibre, …) instead of re-predicting anonymous
    Prodigal ORFs that all render as 'hypothetical / unknown'.

    Returns (sequences, names, feats), each keyed by accession with AND without
    version: ``seqs[acc]`` = nucleotide sequence; ``names[acc]`` = organism/phage
    name; ``feats[acc]`` = ``[(start, end, strand, product, gene, translation)]``
    with 1-based inclusive coordinates (empty list if the record had no CDS).
    """
    Entrez.email = email
    seqs: dict[str, str] = {}
    names: dict[str, str] = {}
    feats: dict[str, list] = {}
    batch = 20  # full GenBank records are larger than FASTA → smaller batches
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        for attempt in range(3):
            try:
                h = Entrez.efetch(db="nuccore", id=",".join(chunk),
                                  rettype="gbwithparts", retmode="text")
                for rec in SeqIO.parse(io.StringIO(h.read()), "genbank"):
                    acc = rec.id
                    base = acc.split(".")[0]
                    s = str(rec.seq)
                    seqs[acc] = s
                    seqs[base] = s
                    nm = (rec.annotations.get("organism") or _clean_name(rec.description) or "").strip()
                    if nm and nm.lower() not in ("", "."):
                        names[acc] = nm
                        names[base] = nm
                    cds = []
                    for f in rec.features:
                        if f.type != "CDS" or f.location is None:
                            continue
                        try:
                            st = int(f.location.start) + 1   # 0-based half-open -> 1-based inclusive
                            en = int(f.location.end)
                        except Exception:
                            continue
                        strand = 1 if (f.location.strand or 1) >= 0 else -1
                        prod = (f.qualifiers.get("product") or [""])[0]
                        gene = (f.qualifiers.get("gene") or [""])[0]
                        transl = (f.qualifiers.get("translation") or [""])[0]
                        cds.append((st, en, strand, prod, gene, transl))
                    cds.sort()
                    feats[acc] = cds
                    feats[base] = cds
                break
            except Exception as e:
                print(f"  efetch gb retry {attempt+1} ({chunk[0]}…): {e}")
                time.sleep(5)
        time.sleep(0.4)  # be polite to NCBI
        print(f"  NCBI fetched {min(i+batch, len(ids))}/{len(ids)}")
    return seqs, names, feats


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
          phage_name: str = "", real_feats: list | None = None) -> Path | None:
    # Flanking-gene annotation: prefer the genome's REAL CDS features (with their
    # /product names) when available (NCBI GenBank records); fall back to a Prodigal
    # gene-call only for metagenomic catalogue genomes, which ship no annotation.
    if real_feats:
        genes = list(real_feats)                                       # (s,e,strand,product,gene,aa)
    else:
        genes = [(s, e, st, "", "", aa) for (s, e, st, aa) in prodigal_genes(genome_id, seq)]
    if not genes:
        return None
    # window = span of (all hit ORFs in this genome) + nearest flanking genes
    h_lo = int(hits.nt_start.min())
    h_hi = int(hits.nt_end.max())
    centre = (h_lo + h_hi) // 2

    def _is_family_call(g):
        # a real gene that IS the family ORF (reciprocal overlap) — drawn separately
        # as the gene of interest, so don't also list it as a flank
        ov = max(0, min(g[1], h_hi) - max(g[0], h_lo))
        return ov > 0.6 * max(1, g[1] - g[0]) and ov > 0.6 * max(1, h_hi - h_lo)
    flank_pool = [g for g in genes if not _is_family_call(g)]
    nearby = sorted(sorted(flank_pool, key=lambda g: abs((g[0]+g[1])//2 - centre))[:FLANKS*2+2],
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
    for (gs, ge, gst, prod, gene, aa) in nearby:
        q = {"product": [prod or "hypothetical protein"]}
        if gene:
            q["gene"] = [gene]
        if aa:
            q["translation"] = [aa]
        feats.append(SeqFeature(FeatureLocation(max(0, gs-lo), max(1, ge-lo), strand=gst),
                                type="CDS", qualifiers=q))
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
    # genome-map figure marking the gene of interest (the HMM hit) among its neighbours,
    # now carrying each flank's REAL product/gene so the map colours + labels them by
    # function instead of collapsing every neighbour to 'hypothetical / unknown'.
    try:
        import genome_map as GM
        a_st = 1 if str(hits.iloc[0].strand) == "+" else -1
        fk = {(g[0], g[1]) for g in nearby}
        GM.draw(GM.build_genes((h_lo, h_hi, a_st),
                               [(g[0], g[1], g[2], {"product": g[3], "gene": g[4]}) for g in nearby],
                               flank_keys=fk),
                (h_lo, h_hi), out_dir / f"{safe}_genome_map",
                "gene of interest (HMM hit) + neighbours",
                log=print, track_name=f"{label}\n{genome_id}",
                tool=os.environ.get("GENOME_MAP_TOOL", "dfv"),
                palette=os.environ.get("GENOME_MAP_PALETTE", "default"),
                functional_labels=os.environ.get("GENOME_MAP_FUNCTIONAL", "") == "1",
                module_brackets=os.environ.get("GENOME_MAP_BRACKETS", "") == "1",
                genbank=out)
    except Exception as e:
        print(f"  (genome map skipped for {genome_id}: {e})")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits-tsv", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--email", default=None,
                    help="NCBI Entrez email for genome retrieval. Never assumed: if omitted "
                         "(and $NCBI_EMAIL unset) NCBI genomes are skipped — no address is "
                         "ever sent to NCBI. Metagenomic-catalogue genomes are unaffected.")
    args = ap.parse_args()
    email = args.email or (os.environ.get("NCBI_EMAIL") or "").strip() or None

    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.hits_tsv, sep="\t")
    # Genome neighbourhoods apply only to six-frame ORF hits (which carry genomic
    # coordinates). Protein-DB hits (e.g. RefSeq YP_ accessions) have no genome
    # here — skip them rather than failing to fetch them as nuccore. (Their parent
    # genome could be resolved via the protein record's coded_by; future work.)
    if "source_type" in df.columns:
        before = len(df)
        df = df[df["source_type"] == "six_frame_orf"].copy()
        skipped = before - len(df)
        if skipped:
            print(f"Skipping {skipped} protein-DB hit(s) without genomic coordinates "
                  "(no neighbourhood to draw)")
    if df.empty:
        print("No six-frame hits with coordinates; no GenBank neighbourhoods to build.")
        return
    by_genome = {gid: sub for gid, sub in df.groupby("genome_id")}
    genome_ids = list(by_genome)
    print(f"{len(genome_ids)} genomes ({len(df)} hits)")

    # 1. NCBI genomes (sequences + phage names). Skipped offline (no real email):
    #    never send a placeholder address to NCBI; metagenomic catalogues still run.
    ncbi_ids = [g for g in genome_ids if is_ncbi(g)]
    if ncbi_ids and email:
        print(f"Fetching {len(ncbi_ids)} NCBI genomes (with annotations)…")
        seqs, names, feats = fetch_ncbi(ncbi_ids, email)
    else:
        if ncbi_ids:
            print(f"(offline: skipping {len(ncbi_ids)} NCBI genome(s) — no --email/$NCBI_EMAIL; "
                  "their neighbourhoods are omitted)")
        seqs, names, feats = {}, {}, {}

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
        rf = feats.get(gid) or feats.get(gid.split(".")[0])
        if build(gid, seq, by_genome[gid], args.out_dir, phage_name=nm, real_feats=rf):
            built += 1
            named += bool(nm)
    print(f"\nBuilt {built} real-sequence GenBank files in {args.out_dir} "
          f"({named} with phage names; {built-named} uncultured/metagenomic)")
    if missing:
        print(f"Could not retrieve {len(missing)} genomes: {missing}")


if __name__ == "__main__":
    main()
