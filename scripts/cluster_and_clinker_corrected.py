#!/usr/bin/env python3
"""
cluster_and_clinker_corrected.py
================================
Cluster the CORRECT family domains and build clinker gene-neighbourhood figures
grouped by those clusters, with the central "hit gene" being the genuine family
six-frame ORF (not the overlapping Prodigal gene that the earlier, defective
output used).

WHY: the original clustering/clinker were computed on the wrong protein file
(Prodigal CDS overlapping each locus). This script rebuilds both on the
ORF-validated family domains produced by extract_validated_hits.py.

STEPS
-----
  1. CD-HIT cluster <validated>/hits_unique_aa.faa (40% id, 80% cov).
  2. Map every hit (from hits.tsv) to a cluster by sequence identity to a
     cluster representative.
  3. For each hit, build a GenBank neighbourhood from its genome:
        - Prodigal genes provide the flanking CDS context (5 up + 5 down)
        - the CENTRAL CDS is the real family ORF (coords + validated translation)
  4. Group GenBanks by cluster; run clinker per cluster (>=2 loci); write an
     index.html linking every cluster figure.

INPUT
-----
  --validated-dir : a run's validated/ dir (hits.tsv, hits_unique_aa.faa)
  --cache-dir     : synteny_context_cache/ with <genome>.fna files
  --out-dir       : output dir for clusters + clinker figures

USAGE
-----
  python3 cluster_and_clinker_corrected.py \
      --validated-dir runA/benchmark/validated \
      --cache-dir runA/benchmark/results/synteny_context_cache \
      --out-dir runA/downstream/clinker
"""
from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Ensure conda env tools (prodigal, cd-hit, clinker) are on PATH regardless of
# how this script was invoked.
_cand = ([Path(os.environ["CONDA_PREFIX"]) / "bin"] if os.environ.get("CONDA_PREFIX") else []) \
    + [Path.home() / _n / "envs" / "hmm-discovery" / "bin"
       for _n in ("miniforge3", "mambaforge", "miniconda3", "anaconda3")]
for _b in _cand:
    if _b.is_dir():
        os.environ["PATH"] = f"{_b}{os.pathsep}{os.environ.get('PATH', '')}"
        break

FLANKS = 5


# ----------------------------------------------------------------------------
# CD-HIT clustering of the correct family domains
# ----------------------------------------------------------------------------
def cdhit(faa: Path, out_prefix: Path) -> dict[int, list[tuple[str, bool]]]:
    """Run CD-HIT; return {cluster_id: [(member_id, is_representative), ...]}."""
    subprocess.run(
        ["cd-hit", "-i", str(faa), "-o", str(out_prefix),
         "-c", "0.4", "-n", "2", "-M", "0", "-T", "8", "-aL", "0.8", "-d", "0"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
    )
    clusters: dict[int, list[tuple[str, bool]]] = {}
    cur = None
    for line in Path(str(out_prefix) + ".clstr").read_text().splitlines():
        if line.startswith(">Cluster"):
            cur = int(line.split()[-1])
            clusters[cur] = []
        elif cur is not None and ">" in line:
            name = line.split(">")[1].split("...")[0]
            clusters[cur].append((name, "*" in line))
    return clusters


# ----------------------------------------------------------------------------
# Prodigal genes per genome (cached)
# ----------------------------------------------------------------------------
_pg: dict[str, list[tuple[int, int, int, str]]] = {}


def prodigal_genes(fna: Path):
    """Return [(start,end,strand_int,aa), ...] for a genome (cached)."""
    key = str(fna)
    if key in _pg:
        return _pg[key]
    genes = []
    with tempfile.NamedTemporaryFile(suffix=".faa", delete=False) as tg:
        gff = tg.name + ".gff"
        faa = tg.name
    subprocess.run(["prodigal", "-i", str(fna), "-a", faa, "-o", gff,
                    "-f", "gff", "-p", "meta", "-q"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    prots = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(faa, "fasta")} \
        if Path(faa).exists() else {}
    if Path(gff).exists():
        idx = 0
        for line in Path(gff).read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            if len(p) >= 7 and p[2] == "CDS":
                idx += 1
                # prodigal protein ids are like <contig>_<idx>
                pid = f"{p[0]}_{idx}"
                aa = prots.get(pid, "")
                genes.append((int(p[3]), int(p[4]), 1 if p[6] == "+" else -1, aa))
    Path(faa).unlink(missing_ok=True)
    Path(gff).unlink(missing_ok=True)
    genes.sort()
    _pg[key] = genes
    return genes


# ----------------------------------------------------------------------------
# Build a GenBank neighbourhood with the real family ORF as the central gene
# ----------------------------------------------------------------------------
def build_genbank(row, cache: Path, out_dir: Path) -> Path | None:
    contig = str(row["contig"])
    fna = cache / f"{contig}.fna"
    if not fna.exists():
        fna = cache / f"{row['genome_id']}.fna"
    if not fna.exists():
        return None
    genes = prodigal_genes(fna)
    if not genes:
        return None

    h_start, h_end = int(row["nt_start"]), int(row["nt_end"])
    h_strand = 1 if row["strand"] == "+" else -1
    centre = (h_start + h_end) // 2

    # nearest flanking Prodigal genes by midpoint distance, excluding ones that
    # essentially ARE the antisense overlap of the hit (keep them — they are
    # real neighbours/context)
    ordered = sorted(genes, key=lambda g: abs((g[0] + g[1]) // 2 - centre))
    nearby = sorted(ordered[: FLANKS * 2 + 2], key=lambda g: g[0])

    window_lo = min([h_start] + [g[0] for g in nearby])
    window_hi = max([h_end] + [g[1] for g in nearby])
    rec = SeqRecord(Seq("N" * (window_hi - window_lo + 1)),
                    id=str(row["genome_id"])[:16],
                    name=str(row["genome_id"])[:16],
                    description=f"{row['genome_id']} family neighbourhood",
                    annotations={"molecule_type": "DNA", "topology": "linear"})

    feats = []
    # flanking genes
    for (gs, ge, gst, aa) in nearby:
        loc = FeatureLocation(gs - window_lo, ge - window_lo, strand=gst)
        feats.append(SeqFeature(loc, type="CDS", qualifiers={
            "product": ["flanking CDS"],
            "translation": [aa or "X"],
        }))
    # central family gene (the real ORF)
    loc = FeatureLocation(h_start - window_lo, h_end - window_lo, strand=h_strand)
    feats.append(SeqFeature(loc, type="CDS", qualifiers={
        "product": ["homologue (HMM hit)"],
        "gene": ["family"],
        "translation": [str(row["aa_sequence"]).rstrip("*") or "X"],
    }))
    rec.features = sorted(feats, key=lambda f: int(f.location.start))

    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in str(row["genome_id"]))[:40]
    out = out_dir / f"{safe}.gbk"
    SeqIO.write(rec, str(out), "genbank")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--validated-dir", type=Path, required=True)
    ap.add_argument("--cache-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    gbk_dir = out / "genbank_files"
    gbk_dir.mkdir(exist_ok=True)

    hits = pd.read_csv(args.validated_dir / "hits.tsv", sep="\t")
    unique_faa = args.validated_dir / "hits_unique_aa.faa"

    # 1. cluster the correct domains
    clusters = cdhit(unique_faa, out / "clusters")
    print(f"CD-HIT: {len(clusters)} clusters from {sum(len(v) for v in clusters.values())} unique domains")

    # representative sequence -> cluster id
    uniq = {r.id: str(r.seq) for r in SeqIO.parse(str(unique_faa), "fasta")}
    repid_to_cluster = {}
    for cid, members in clusters.items():
        for mid, is_rep in members:
            repid_to_cluster[mid] = cid

    # 2. map every hit to a cluster (by exact sequence match to a unique rep,
    #    else by identical sequence)
    seq_to_cluster = {uniq[mid]: cid for cid, members in clusters.items()
                      for mid, _ in members if mid in uniq}
    # write cluster membership table
    mem_rows = []

    # 3. build GenBanks, grouped by cluster
    by_cluster: dict[int, list[Path]] = defaultdict(list)
    for _, row in hits.iterrows():
        seq = str(row["aa_sequence"])
        cid = seq_to_cluster.get(seq)
        if cid is None:
            # assign by membership of its own id if present
            cid = repid_to_cluster.get(str(row["hit_id"]))
        if cid is None:
            continue
        gbk = build_genbank(row, args.cache_dir, gbk_dir)
        if gbk:
            by_cluster[cid].append(gbk)
            mem_rows.append({"hit_id": row["hit_id"], "genome_id": row["genome_id"],
                             "cluster_id": cid, "db_name": row["db_name"]})

    pd.DataFrame(mem_rows).to_csv(out / "cluster_membership.tsv", sep="\t", index=False)

    # 4. clinker per cluster
    figdir = out / "clinker_figures"
    figdir.mkdir(exist_ok=True)
    produced = {}
    for cid in sorted(by_cluster):
        gbks = sorted(set(by_cluster[cid]))
        if len(gbks) < 2:
            continue
        html = figdir / f"cluster_{cid}.html"
        try:
            subprocess.run(["clinker", *[str(g) for g in gbks[:30]],
                            "-p", str(html), "-i", "0.3", "-j", "4"],
                           capture_output=True, text=True, timeout=600)
            if html.exists() and html.stat().st_size > 100:
                produced[cid] = html
                print(f"  cluster_{cid}: {len(gbks)} loci -> {html.stat().st_size//1024} KB")
        except Exception as e:
            print(f"  cluster_{cid}: clinker failed: {e}")

    # index
    idx = ["<!DOCTYPE html><html><head><meta charset='utf-8'>",
           "<title>family clinker (corrected)</title>",
           "<style>body{font-family:sans-serif;margin:2em}td,th{border:1px solid #ddd;padding:6px 12px}</style>",
           "</head><body><h1>family synteny by corrected cluster</h1><table>",
           "<tr><th>Cluster</th><th>Loci</th><th>Figure</th></tr>"]
    for cid in sorted(by_cluster):
        n = len(set(by_cluster[cid]))
        link = f'<a href="clinker_figures/cluster_{cid}.html">open</a>' if cid in produced else "&mdash;"
        idx.append(f"<tr><td>{cid}</td><td>{n}</td><td>{link}</td></tr>")
    idx.append("</table></body></html>")
    (out / "index.html").write_text("\n".join(idx))
    print(f"Done. {len(produced)} clinker figures. Index: {out/'index.html'}")


if __name__ == "__main__":
    main()
