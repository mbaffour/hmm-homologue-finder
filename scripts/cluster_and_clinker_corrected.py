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
import re
import socket
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

socket.setdefaulttimeout(60)  # bound NCBI coded_by/genome fetches; never hang a run
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Ensure conda env tools (prodigal, cd-hit, clinker) are on PATH regardless of
# how this script was invoked.
from env_paths import ensure_env_on_path  # noqa: E402  (sibling helper in scripts/)
ensure_env_on_path()

_BROWSER_WARNED = [False]


def _html_to_png(html: Path, png: Path) -> bool:
    """Render a clinker HTML plot to a static PNG via headless chromium (playwright),
    when available. clinker's plot is JS-rendered (clustermap.js), so a static image
    needs a browser. No-op + a one-time note when none is installed — the SVG/PNG/PDF
    synteny panels from synteny_figure.py remain the primary static figures.
    Enable with: pip install playwright && playwright install chromium."""
    try:
        from playwright.sync_api import sync_playwright
    except Exception:
        if not _BROWSER_WARNED[0]:
            print("  (static clinker PNG skipped: no headless browser — enable with "
                  "`pip install playwright && playwright install chromium`; the static "
                  "synteny panels in downstream/synteny/ are produced regardless)")
            _BROWSER_WARNED[0] = True
        return False
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_page(viewport={"width": 1600, "height": 1000})
            page.goto(html.resolve().as_uri())
            page.wait_for_timeout(2500)        # let clustermap.js draw the SVG
            page.screenshot(path=str(png), full_page=True)
            browser.close()
        return png.exists() and png.stat().st_size > 1000
    except Exception as e:
        if not _BROWSER_WARNED[0]:
            print(f"  (static clinker PNG skipped: {e})")
            _BROWSER_WARNED[0] = True
        return False

FLANKS = 7  # ORFs each side of the family gene (publication default for phage neighbourhoods)


# ----------------------------------------------------------------------------
# CD-HIT clustering of the correct family domains
# ----------------------------------------------------------------------------
# Flags for THIS step (neighbourhood grouping), pinned as constants so the figure
# grouping is reproducible and so a second caller cannot silently drift from it:
#   -c 0.4  remote-homology band — these are divergent phage proteins, not orthologs
#   -n 2    cd-hit's own required word size for the 0.4-0.5 identity band
#   -aL 0.8 the alignment must cover 80 % of the LONGER sequence, so a short domain
#           cannot be absorbed into an unrelated long protein at only 40 % identity
CDHIT_IDENT = 0.4
CDHIT_WORD = 2
CDHIT_ALIGN_LONGER = 0.8
CDHIT_THREADS = 8


def parse_clstr(text: str) -> list[list[tuple[str, bool]]]:
    """Parse cd-hit `.clstr` text -> one ``[(member_name, is_representative), ...]``
    list per cluster, in file order.

    Split out from `cdhit` because this is where the bugs live — cd-hit truncates
    every name with a literal ``...`` and marks the cluster representative with a
    trailing ``*`` — and parsing needs no cd-hit binary, so it can be unit-tested.

    cd-hit numbers clusters sequentially from 0 in file order, so a cluster's INDEX
    in the returned list is its cd-hit cluster id; callers can `enumerate` it.
    """
    clusters: list[list[tuple[str, bool]]] = []
    cur = None
    for line in str(text or "").splitlines():
        if line.startswith(">Cluster"):
            cur = []
            clusters.append(cur)
        elif cur is not None and ">" in line:
            name = line.split(">")[1].split("...")[0]
            cur.append((name, "*" in line))
    return clusters


def cdhit(faa: Path, out_prefix: Path,
          ident: float = CDHIT_IDENT, word: int = CDHIT_WORD,
          aL: float = CDHIT_ALIGN_LONGER, cpu: int = CDHIT_THREADS,
          ) -> dict[int, list[tuple[str, bool]]]:
    """Run CD-HIT; return {cluster_id: [(member_id, is_representative), ...]}.

    The defaults reproduce the flags this step has always used, so the clinker
    grouping is byte-identical to before. The parameters exist so a second caller
    (`family_census.py`, which clusters the seeds-plus-hits union at 100/95/90 %)
    goes through ONE cd-hit wrapper instead of growing a divergent copy.
    """
    subprocess.run(
        ["cd-hit", "-i", str(faa), "-o", str(out_prefix),
         "-c", f"{ident:g}", "-n", str(word), "-M", "0", "-T", str(cpu),
         "-aL", f"{aL:g}", "-d", "0"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
    )
    return dict(enumerate(parse_clstr(Path(str(out_prefix) + ".clstr").read_text())))


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
# Protein-DB hits: resolve parent genome + CDS coords via the GenPept coded_by,
# fetch the genome into the synteny cache so they get a real neighbourhood too.
# ----------------------------------------------------------------------------
def _parse_coded_by(genpept_text: str):
    """'…coded_by="complement(NC_019520.1:37102..37269)"…' -> (acc, start, end, strand)."""
    m = re.search(r'coded_by="?(complement\()?([A-Za-z0-9_.]+):[<>]?(\d+)\.\.[<>]?(\d+)',
                  genpept_text)
    if not m:
        return None
    return m.group(2), int(m.group(3)), int(m.group(4)), ("-" if m.group(1) else "+")


def resolve_protein_neighborhoods(hits, cache_dir: Path, email: str) -> int:
    """For protein-DB hits, look up the parent genome + CDS coords (coded_by) and
    fetch the genome into cache_dir, filling in coords so build_genbank can draw a
    neighbourhood. Mutates `hits` in place. Degrades gracefully (skips on any
    NCBI/parse failure); returns the count resolved."""
    if "source_type" not in hits.columns:
        return 0
    Entrez.email = email
    resolved = 0
    for i in hits.index[hits["source_type"] == "annotated_protein"]:
        acc = str(hits.at[i, "genome_id"]).strip()
        if not acc:
            continue
        try:
            gb = Entrez.efetch(db="protein", id=acc, rettype="gb", retmode="text").read()
        except Exception:
            continue
        cb = _parse_coded_by(gb)
        if not cb:
            continue
        gacc, gs, ge, strand = cb
        fna = cache_dir / f"{gacc}.fna"
        if not (fna.exists() and fna.stat().st_size > 0):
            try:
                fna.write_text(Entrez.efetch(db="nuccore", id=gacc,
                                             rettype="fasta", retmode="text").read())
            except Exception:
                continue
        if not (fna.exists() and fna.stat().st_size > 0):
            continue
        hits.at[i, "contig"] = gacc
        hits.at[i, "nt_start"] = gs
        hits.at[i, "nt_end"] = ge
        hits.at[i, "strand"] = strand
        resolved += 1
    if resolved:
        print(f"Resolved {resolved} protein-DB hit(s) to genome neighbourhoods via coded_by")
    return resolved


def dedup_synteny_loci(hits) -> set:
    """Indices to KEEP for synteny: collapse hits mapping to the SAME genomic
    locus (same parent accession, overlapping coordinates) — e.g. a six-frame ORF
    and the RefSeq protein of the same gene (resolved via coded_by). Prefer the
    six-frame ORF, then higher bit score. Genuine paralogs (non-overlapping loci
    on one genome) are kept separate. Hits without coordinates are left untouched
    (they build nothing anyway)."""
    rows = []
    for i in hits.index:
        contig = str(hits.at[i, "contig"]).split(".")[0].strip()
        try:
            s, e = int(hits.at[i, "nt_start"]), int(hits.at[i, "nt_end"])
        except (ValueError, TypeError):
            continue
        if not contig:
            continue
        try:
            bs = float(hits.at[i, "bit_score"])
        except (ValueError, TypeError):
            bs = 0.0
        rows.append((i, contig, min(s, e), max(s, e),
                     str(hits.at[i, "source_type"]), bs))
    keep = set()
    by_contig = defaultdict(list)
    for r in rows:
        by_contig[r[1]].append(r)
    for _contig, rs in by_contig.items():
        rs.sort(key=lambda r: r[2])
        clusters = []  # [lo, hi, [members]]
        for r in rs:
            lo, hi = r[2], r[3]
            for c in clusters:
                if lo <= c[1] and hi >= c[0]:        # overlapping -> same locus
                    c[0], c[1] = min(c[0], lo), max(c[1], hi)
                    c[2].append(r)
                    break
            else:
                clusters.append([lo, hi, [r]])
        for c in clusters:
            best = max(c[2], key=lambda r: (r[4] == "six_frame_orf", r[5]))
            keep.add(best[0])
    return keep


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
    genome_id = str(row["genome_id"])
    organism = str(row.get("organism", "") or "").strip()
    label = f"{organism} ({genome_id})" if organism else genome_id
    locus = "".join(c if c.isalnum() or c in "._-" else "_" for c in genome_id)[:16]  # GenBank LOCUS limit
    rec = SeqRecord(Seq("N" * (window_hi - window_lo + 1)),
                    id=locus,
                    name=locus,
                    description=label,
                    annotations={"molecule_type": "DNA", "topology": "linear",
                                 "organism": organism or genome_id})

    feats = []
    # flanking genes
    for (gs, ge, gst, aa) in nearby:
        loc = FeatureLocation(gs - window_lo, ge - window_lo, strand=gst)
        feats.append(SeqFeature(loc, type="CDS", qualifiers={
            "product": ["flanking CDS"],
            "translation": [aa or "X"],
        }))
    # central family gene (the real ORF) — consistent name so clinker colours and
    # links it identically across every locus, making the family gene easy to spot.
    loc = FeatureLocation(h_start - window_lo, h_end - window_lo, strand=h_strand)
    feats.append(SeqFeature(loc, type="CDS", qualifiers={
        "product": ["FAMILY HOMOLOGUE (HMM hit)"],
        "gene": ["family_homologue"],
        "translation": [str(row["aa_sequence"]).rstrip("*") or "X"],
    }))
    rec.features = sorted(feats, key=lambda f: int(f.location.start))

    # File name drives the clinker track label — use the organism so tracks read
    # as phage names, not bare accessions.
    safe = "".join(c if c.isalnum() or c in "._- " else "_" for c in label).strip()
    safe = "_".join(safe.split())[:60] or genome_id
    out = out_dir / f"{safe}.gbk"
    SeqIO.write(rec, str(out), "genbank")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--validated-dir", type=Path, required=True)
    ap.add_argument("--cache-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--email", default="",
                    help="NCBI email; if set, protein-DB hits are resolved to genome "
                         "neighbourhoods via coded_by")
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    gbk_dir = out / "genbank_files"
    gbk_dir.mkdir(exist_ok=True)

    hits = pd.read_csv(args.validated_dir / "hits.tsv", sep="\t")
    unique_faa = args.validated_dir / "hits_unique_aa.faa"

    # Give protein-DB hits a neighbourhood too: resolve their parent genome +
    # CDS coords via coded_by and fetch the genome into the synteny cache.
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    if args.email and "annotated_protein" in set(hits.get("source_type", [])):
        try:
            resolve_protein_neighborhoods(hits, args.cache_dir, args.email)
        except Exception as e:
            print(f"(protein coded_by resolution skipped: {e})")

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

    # 3. build GenBanks, grouped by cluster. First collapse hits that map to the
    #    same genomic locus (e.g. a six-frame ORF and the RefSeq protein of the
    #    same gene) so a locus is drawn once, not double-counted.
    keep_idx = dedup_synteny_loci(hits)
    n_coord = sum(1 for i in hits.index
                  if str(hits.at[i, "nt_start"]).strip() not in ("", "nan"))
    if keep_idx and len(keep_idx) < n_coord:
        print(f"Synteny dedup: {n_coord} coordinate hits -> {len(keep_idx)} distinct loci")
    by_cluster: dict[int, list[Path]] = defaultdict(list)
    for idx, row in hits.iterrows():
        if keep_idx and idx not in keep_idx:
            continue
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

    # 4. clinker per cluster. Cap the number of tracks per figure so it stays
    #    legible: a figure with 30-80 genome tracks is unreadable. When a cluster
    #    is larger, evenly sample MAX_LOCI loci across it (the representative is
    #    first, so it is always kept) for a readable, diverse figure.
    MAX_LOCI = 16
    figdir = out / "clinker_figures"
    figdir.mkdir(exist_ok=True)
    produced = {}
    static_png = {}
    shown_counts = {}
    for cid in sorted(by_cluster):
        gbks_all = sorted(set(by_cluster[cid]))
        if len(gbks_all) < 2:
            continue
        if len(gbks_all) > MAX_LOCI:
            step = len(gbks_all) / MAX_LOCI
            gbks = [gbks_all[int(i * step)] for i in range(MAX_LOCI)]
        else:
            gbks = gbks_all
        html = figdir / f"cluster_{cid}.html"
        try:
            subprocess.run(["clinker", *[str(g) for g in gbks],
                            "-p", str(html), "-i", "0.3", "-j", "4"],
                           capture_output=True, text=True, timeout=600)
            if html.exists() and html.stat().st_size > 100:
                produced[cid] = html
                shown_counts[cid] = (len(gbks), len(gbks_all))
                note = f"showing {len(gbks)} of {len(gbks_all)}" if len(gbks) < len(gbks_all) else f"{len(gbks)} loci"
                print(f"  cluster_{cid}: {note} -> {html.stat().st_size // 1024} KB")
                png = figdir / f"cluster_{cid}.png"
                if _html_to_png(html, png):
                    static_png[cid] = png
        except Exception as e:
            print(f"  cluster_{cid}: clinker failed: {e}")

    # index
    idx = ["<!DOCTYPE html><html><head><meta charset='utf-8'>",
           "<title>family clinker (corrected)</title>",
           "<style>body{font-family:sans-serif;margin:2em}td,th{border:1px solid #ddd;padding:6px 12px}</style>",
           "</head><body><h1>family synteny by corrected cluster</h1>",
           "<p>The central <b>FAMILY HOMOLOGUE</b> gene is named identically in every "
           "track so clinker colours/links it consistently. Tracks are labelled by "
           "organism. Large clusters show an evenly-sampled, readable subset.</p><table>",
           "<tr><th>Cluster</th><th>Loci shown / total</th><th>Interactive</th><th>Static PNG</th></tr>"]
    for cid in sorted(by_cluster):
        n = len(set(by_cluster[cid]))
        shown = shown_counts.get(cid, (0, n))[0]
        cell = f"{shown} / {n}" if cid in produced else f"&mdash; / {n}"
        link = f'<a href="clinker_figures/cluster_{cid}.html">open</a>' if cid in produced else "&mdash;"
        spng = f'<a href="clinker_figures/cluster_{cid}.png">PNG</a>' if cid in static_png else "&mdash;"
        idx.append(f"<tr><td>{cid}</td><td>{cell}</td><td>{link}</td><td>{spng}</td></tr>")
    idx.append("</table></body></html>")
    (out / "index.html").write_text("\n".join(idx))
    print(f"Done. {len(produced)} clinker figures ({len(static_png)} static PNG). "
          f"Index: {out/'index.html'}")


if __name__ == "__main__":
    main()
