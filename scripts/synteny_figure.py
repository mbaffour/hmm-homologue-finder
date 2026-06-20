#!/usr/bin/env python3
"""
synteny_figure.py — publication-quality static synteny panels (one per cluster).

Built the way synteny is shown in papers:
  * every locus is ANCHORED on the family gene (aligned column) and strand-
    normalised so the family gene always points right;
  * neighbourhood genes are FUNCTIONALLY ANNOTATED (VOGDB VFAM, hmmscan) and
    coloured by function, so each colour is a real gene function named in the
    legend (not an arbitrary number);
  * homologous genes are joined by shaded links between adjacent loci;
  * the family gene is highlighted and labelled inline; tracks are labelled by
    organism; a scale bar and a legend (outside the plot — no overlapping text)
    are included.
Outputs SVG + PNG + PDF per cluster, a per-gene annotation CSV, and index.html.
Needs matplotlib, Biopython, CD-HIT (always), and VOGDB for function names
(downloaded once into the shared db-cache; falls back to "hypothetical protein"
labels if unavailable).
"""
from __future__ import annotations

import argparse
import csv
import difflib
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
# Keep text EDITABLE in vector exports: SVG keeps real <text> elements (editable
# in Inkscape), PDF/PS embed TrueType (selectable/editable in Illustrator) rather
# than converting glyphs to outlined paths.
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
from Bio import SeqIO

from env_paths import ensure_env_on_path  # noqa: E402
import annotate_genes  # noqa: E402
ensure_env_on_path()

FAMILY_COLOR = "#ffd400"          # bold gold for the gene of interest (distinct from all categories)
HYPO_COLOR = "#dde2e8"
MAX_LOCI = 12

# Broad functional categories (consistent colours across every figure). Genes are
# binned by keyword from their VOGDB function description; first match wins.
CATEGORY_RULES = [
    ("lysis", ["holin", "endolysin", "spanin", "lysin", "lysis", "amidase",
               "peptidoglycan", "murein", "lysozyme"]),
    ("packaging", ["terminase", "packaging", "headful"]),
    # checked before replication so "(DNA-directed) RNA polymerase" -> transcription,
    # while "DNA polymerase" still falls through to replication
    ("transcription / regulation",
     ["rna polymerase", "rnap", "transcription", "regulator", "repressor",
      "sigma factor", "anti-sigma", "transcriptional"]),
    # specific structural terms only — NOT bare "virion" (too greedy: it caught
    # "virion DNA-directed RNA polymerase")
    ("structural", ["capsid", "tail", "baseplate", "neck", "collar", "fiber",
                    "fibre", "sheath", "tape measure", "decoration", "portal",
                    "prohead", "scaffold", "spike", "head-tail", "head completion",
                    "major head", "structural protein", "virion structural"]),
    ("replication / nucleotide metabolism",
     ["polymerase", "primase", "helicase", "ligase", "nuclease", "exonuclease",
      "endonuclease", "recombinase", "ribonucleotide", "thymidylate",
      "topoisomerase", "single strand", "single-strand", "ssb", "replication",
      "dna binding", "resolvase", "reductase"]),
    ("integration / mobile",
     ["integrase", "transposase", "excisionase", "mobile", "intron", "homing"]),
    ("host / metabolism / defense",
     ["methyltransferase", "transferase", "metabolism", "moron", "toxin",
      "restriction", "defense", "thoeris", "hydrolase", "oxidoreductase", "deaminase"]),
    ("RNA gene", ["trna", "tmrna", "ncrna", "ribosomal rna"]),
]
CATEGORY_COLORS = {
    "structural": "#3b7dd8",
    "packaging": "#8a5fd0",
    "replication / nucleotide metabolism": "#1d9e75",
    "transcription / regulation": "#e0962a",
    "lysis": "#d1495b",
    "integration / mobile": "#9b6a3a",
    "host / metabolism / defense": "#5dada0",
    "RNA gene": "#c46aa8",
    "hypothetical / unknown": HYPO_COLOR,
}
HYPO_CAT = "hypothetical / unknown"


_GENERIC_ORG = ("uncultured", "unclassified", "unknown", "environmental",
                "metagenom", "virus sp", "sp.")


def set_row_labels(loci: list[dict]) -> None:
    """Set loc['label']: organism, but append the accession when the organism is
    missing, generic (e.g. 'uncultured virus'), or duplicated — so rows are
    distinguishable in the figure."""
    from collections import Counter
    orgs = [(l.get("organism") or "").strip() for l in loci]
    dup = {o for o, c in Counter(orgs).items() if o and c > 1}
    for l, o in zip(loci, orgs):
        gid = l.get("genome_id", "")
        generic = (not o) or o in dup or any(g in o.lower() for g in _GENERIC_ORG)
        l["label"] = (f"{o} ({gid})" if (o and generic) else (o or gid))[:60]


# VOGDB FunctionalCategory single-letter codes -> our broad categories (used only
# as a fallback when the description keywords don't resolve). Ambiguous combos
# (e.g. "XhXpXrXs") and Xp/Xu (poorly characterised / unknown) stay hypothetical.
_VOG_CAT = {"Xs": "structural", "Xr": "replication / nucleotide metabolism",
            "Xh": "host / metabolism / defense"}


def categorize(func: str, vog_cat: str = "") -> str:
    f = (func or "").lower()
    if f and not any(w in f for w in ("hypothetical", "uncharacterized", "unknown", "duf")):
        for cat, kws in CATEGORY_RULES:
            if any(k in f for k in kws):
                return cat
    # keyword mapping didn't resolve -> fall back to VOGDB's own category if unambiguous
    return _VOG_CAT.get((vog_cat or "").strip(), HYPO_CAT)


def parse_locus(gbk: Path) -> dict | None:
    try:
        rec = next(SeqIO.parse(str(gbk), "genbank"))
    except Exception:
        return None
    genes = []
    for f in rec.features:
        if f.type != "CDS":
            continue
        gene_q = "".join(f.qualifiers.get("gene", [""]))
        prod = "".join(f.qualifiers.get("product", [""]))
        genes.append({
            "s": int(f.location.start), "e": int(f.location.end),
            "st": 1 if f.location.strand in (1, None) else -1,
            "aa": "".join(f.qualifiers.get("translation", [""])),
            "fam": (gene_q == "family_homologue") or ("FAMILY HOMOLOGUE" in prod),
            "og": None, "func": "hypothetical protein", "vfam": "",
            "vog_cat": "", "category": HYPO_CAT,
        })
    if not any(g["fam"] for g in genes):
        return None
    genes.sort(key=lambda g: g["s"])
    return {"genome_id": rec.id,
            "organism": rec.annotations.get("organism", rec.id),
            "genes": genes}


def assign_orthogroups(loci: list[dict]) -> None:
    """CD-HIT all neighbourhood proteins; tag each gene with g['og']."""
    with tempfile.TemporaryDirectory() as td:
        fasta = Path(td) / "genes.faa"
        idmap = {}
        with fasta.open("w") as fh:
            for li, loc in enumerate(loci):
                for gi, g in enumerate(loc["genes"]):
                    if g["fam"] or not g["aa"]:
                        continue
                    uid = f"{li}_{gi}"
                    idmap[uid] = (li, gi)
                    fh.write(f">{uid}\n{g['aa']}\n")
        if not idmap:
            return
        out = Path(td) / "og"
        try:
            subprocess.run(["cd-hit", "-i", str(fasta), "-o", str(out), "-c", "0.4",
                            "-n", "2", "-aL", "0.6", "-d", "0", "-M", "0", "-T", "4"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            cur = -1
            for line in (Path(str(out) + ".clstr")).read_text().splitlines():
                if line.startswith(">Cluster"):
                    cur = int(line.split()[-1])
                elif ">" in line:
                    uid = line.split(">")[1].split("...")[0]
                    if uid in idmap:
                        li, gi = idmap[uid]
                        loci[li]["genes"][gi]["og"] = cur
        except Exception:
            pass


def consensus_functions(loci: list[dict]) -> None:
    """Set each gene's func to its orthogroup's consensus function."""
    byog = defaultdict(list)
    for loc in loci:
        for g in loc["genes"]:
            if g["og"] is not None:
                byog[g["og"]].append(g)
    for og, genes in byog.items():
        named = [g["func"] for g in genes
                 if g["func"] and g["func"] != "hypothetical protein"]
        if named:
            cons = Counter(named).most_common(1)[0][0]
            for g in genes:
                g["func"] = cons


def anchor(loc: dict) -> dict | None:
    genes = loc["genes"]
    span = max(g["e"] for g in genes)
    fam = next((g for g in genes if g["fam"]), None)
    if fam is None:
        return None
    if fam["st"] == -1:
        for g in genes:
            g["s"], g["e"], g["st"] = span - g["e"], span - g["s"], -g["st"]
        genes.sort(key=lambda g: g["s"])
        fam = next(g for g in genes if g["fam"])
    off = -fam["s"]
    for g in genes:
        g["s"] += off
        g["e"] += off
    return loc


_SIM_CACHE: dict = {}


def _similarity(a: str, b: str) -> float:
    """Fast 0..1 sequence-similarity proxy (cached) for grading homology links."""
    if not a or not b:
        return 0.0
    key = (a, b) if a <= b else (b, a)
    if key not in _SIM_CACHE:
        _SIM_CACHE[key] = difflib.SequenceMatcher(None, a, b).ratio()
    return _SIM_CACHE[key]


def _conservation(loci: list[dict]) -> dict:
    """orthogroup -> fraction of loci that contain it (0..1)."""
    n = len(loci) or 1
    seen = defaultdict(set)
    for i, loc in enumerate(loci):
        for g in loc["genes"]:
            if g["og"] is not None:
                seen[g["og"]].add(i)
    return {og: len(idx) / n for og, idx in seen.items()}


def _arrow(ax, g, y, color, h=0.34, lw=0.4, ec="#33373d"):
    x0, x1, st = g["s"], g["e"], g["st"]
    head = min((x1 - x0) * 0.45, 220)
    if st >= 0:
        b = x1 - head
        pts = [(x0, y - h / 2), (b, y - h / 2), (x1, y), (b, y + h / 2), (x0, y + h / 2)]
    else:
        b = x0 + head
        pts = [(x1, y - h / 2), (b, y - h / 2), (x0, y), (b, y + h / 2), (x1, y + h / 2)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor=ec,
                         linewidth=lw, zorder=4 if lw > 1 else 3))


def draw_cluster(cid, loci, out_dir, color_by="function", suffix=""):
    loci = [a for a in (anchor(l) for l in loci) if a]
    if len(loci) < 2:
        return None
    set_row_labels(loci)
    n = len(loci)
    cons = _conservation(loci) if color_by == "conservation" else {}

    # colour by broad functional CATEGORY (default) or by CONSERVATION gradient
    def color_of(g):
        if g["fam"]:
            return FAMILY_COLOR
        if color_by == "conservation":
            return plt.cm.Blues(0.25 + 0.7 * cons.get(g["og"], 1.0 / n))
        return CATEGORY_COLORS.get(g["category"], HYPO_COLOR)

    cats_present = [c for c in CATEGORY_COLORS
                   if any((not g["fam"]) and g["category"] == c for l in loci for g in l["genes"])]

    big = n > 20
    row_h = 0.32 if big else 0.6        # shorter rows when there are many loci
    arrow_h = 0.24 if big else 0.34
    lab_fs = 6 if big else 8
    xmin = min(g["s"] for l in loci for g in l["genes"])
    xmax = max(g["e"] for l in loci for g in l["genes"])
    fig, ax = plt.subplots(figsize=(12, max(2.4, row_h * n + 1.8)))

    for i, l in enumerate(loci):
        y = n - i
        ax.hlines(y, xmin, xmax, color="#c9ced6", lw=0.6, zorder=0)
        for g in l["genes"]:
            fam = g["fam"]
            _arrow(ax, g, y, color_of(g), h=arrow_h * (1.25 if fam else 1.0),
                   lw=1.8 if fam else 0.4, ec="#1a1a1a" if fam else "#33373d")
        if not big:  # inline "family" tag clutters very tall figures; legend covers it
            fam = next(g for g in l["genes"] if g["fam"])
            ax.text((fam["s"] + fam["e"]) / 2, y + 0.34, "family", ha="center",
                    va="bottom", fontsize=6.5, color=FAMILY_COLOR)

    # homology links (same orthogroup) between adjacent loci
    for i in range(n - 1):
        top, bot = loci[i], loci[i + 1]
        yt, yb = n - i, n - (i + 1)
        bot_by_og = {}
        for g in bot["genes"]:
            if g["og"] is not None:
                bot_by_og.setdefault(g["og"], g)
        for g in top["genes"]:
            o = g["og"]
            if o is None or o not in bot_by_og:
                continue
            gb = bot_by_og[o]
            col = color_of(g)
            sim = _similarity(g.get("aa", ""), gb.get("aa", ""))  # darker link = more similar
            pts = [(g["s"], yt - 0.17), (g["e"], yt - 0.17),
                   (gb["e"], yb + 0.17), (gb["s"], yb + 0.17)]
            ax.add_patch(Polygon(pts, closed=True, facecolor=col,
                                 edgecolor="none", alpha=0.10 + 0.45 * sim, zorder=1))

    pad = (xmax - xmin) * 0.02
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(0.2, n + 1.0)
    ax.set_yticks([n - i for i in range(n)])
    ax.set_yticklabels([l.get("label", l["organism"])[:50] for l in loci], fontsize=lab_fs)
    ax.set_xlabel("position relative to family gene (bp)", fontsize=9)
    ax.set_title(f"Cluster {cid} — gene-neighbourhood synteny "
                 f"(anchored on family gene, n={n})", fontsize=11)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(left=False)

    bar = 1000
    ax.plot([xmin, xmin + bar], [0.45, 0.45], color="#33373d", lw=1.6)
    ax.text(xmin + bar / 2, 0.5, "1 kb", ha="center", va="bottom", fontsize=8)

    handles = [Patch(facecolor=FAMILY_COLOR, edgecolor="#1a1a1a", linewidth=1.4,
                     label="gene of interest")]
    if color_by == "conservation":
        handles += [Patch(facecolor=plt.cm.Blues(0.25 + 0.7 * f), edgecolor="#33373d", label=lab)
                    for f, lab in ((1.0, "conserved (all loci)"), (0.5, "intermediate"),
                                   (1.0 / n, "unique (1 locus)"))]
    else:
        handles += [Patch(facecolor=CATEGORY_COLORS[c], edgecolor="#33373d", label=c)
                    for c in cats_present if c != HYPO_CAT]
        handles += [Patch(facecolor=HYPO_COLOR, edgecolor="#33373d", label=HYPO_CAT)]
    ncol = 3 if len(handles) > 6 else 2
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13),
              ncol=ncol, fontsize=7.5, frameon=False, handlelength=1.4, columnspacing=1.2)

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / f"cluster_{cid}_synteny{suffix}"
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    return Path(f"{stem}.png")


def build_heatmap(by_cluster: dict, out_dir: Path):
    """Neighbourhood functional-composition heatmap: genomes x category counts.
    Shows whether the gene of interest consistently sits among the same kinds of
    genes (conserved module) across genomes. Also writes a conservation CSV."""
    import numpy as np
    rows = []
    for cid in sorted(by_cluster):
        set_row_labels(by_cluster[cid])
        for loc in by_cluster[cid]:
            cnt = Counter(g["category"] for g in loc["genes"] if not g["fam"])
            rows.append((cid, loc.get("label", loc["organism"]), loc["genome_id"], cnt))
    if len(rows) < 2:
        return None
    present = [c for c in CATEGORY_COLORS if any(r[3].get(c, 0) for r in rows)]
    if not present:
        return None
    M = np.array([[r[3].get(c, 0) for c in present] for r in rows], dtype=float)
    nrow = len(rows)
    fig, ax = plt.subplots(figsize=(max(6, 0.65 * len(present) + 3),
                                    max(3, 0.2 * nrow + 1.6)))
    im = ax.imshow(M, aspect="auto", cmap="Blues", vmin=0)
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(present, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(nrow))
    ax.set_yticklabels([r[1][:42] for r in rows], fontsize=5 if nrow > 40 else 7)
    ax.set_title("Neighbourhood functional composition\n(genes per category around the gene of interest)",
                 fontsize=10)
    cb = fig.colorbar(im, ax=ax, shrink=0.4)
    cb.set_label("genes in window", fontsize=8)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out_dir / f"neighbourhood_heatmap.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    with (out_dir / "neighbourhood_conservation.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cluster", "organism", "genome_id"] + present)
        for cid, org, gid, cnt in rows:
            w.writerow([cid, org, gid] + [cnt.get(c, 0) for c in present])
    return out_dir / "neighbourhood_heatmap.png"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clinker-dir", type=Path, required=True,
                    help="downstream/clinker dir (has genbank_files/ + cluster_membership.tsv)")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--annotation-cache", type=Path,
                    default=Path.home() / ".cache" / "hmm-homologue-finder",
                    help="shared cache holding the VOGDB VFAM annotation database")
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--max-loci", type=int, default=MAX_LOCI,
                    help="max loci drawn per cluster figure; 0 = all loci on one figure")
    ap.add_argument("--color-by", choices=("function", "conservation", "both"), default="function",
                    help="gene colouring: 'function' (category), 'conservation' "
                         "(core->unique gradient), or 'both'. Links shaded by similarity in all.")
    args = ap.parse_args()

    gbk_dir = args.clinker_dir / "genbank_files"
    membership = args.clinker_dir / "cluster_membership.tsv"
    if not gbk_dir.is_dir() or not membership.exists():
        print(f"  (synteny figures skipped: missing {gbk_dir} or {membership})")
        return

    gid_to_cluster = {}
    with membership.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid_to_cluster.setdefault(str(r.get("genome_id", "")), r.get("cluster_id", ""))

    by_cluster: dict[int, list[dict]] = {}
    for gbk in sorted(gbk_dir.glob("*.gbk")):
        loc = parse_locus(gbk)
        if not loc:
            continue
        cid = gid_to_cluster.get(loc["genome_id"])
        if cid in (None, ""):
            continue
        by_cluster.setdefault(int(cid), []).append(loc)

    # subsample each cluster to a readable, representative set (0 = keep all)
    cap = args.max_loci
    for cid in list(by_cluster):
        loci = by_cluster[cid]
        if cap and len(loci) > cap:
            step = len(loci) / cap
            by_cluster[cid] = [loci[int(i * step)] for i in range(cap)]

    # functional annotation (VOGDB VFAM) of every neighbourhood protein, once
    ann_ready = annotate_genes.is_ready(args.annotation_cache)
    if ann_ready:
        proteins = {}
        for cid, loci in by_cluster.items():
            for li, loc in enumerate(loci):
                for gi, g in enumerate(loc["genes"]):
                    if g["aa"] and not g["fam"]:
                        proteins[f"{cid}|{li}|{gi}"] = g["aa"]
        print(f"  annotating {len(proteins)} neighbourhood proteins with VOGDB VFAM…")
        hits = annotate_genes.annotate(proteins, args.annotation_cache, cpu=args.cpu)
        for cid, loci in by_cluster.items():
            for li, loc in enumerate(loci):
                for gi, g in enumerate(loc["genes"]):
                    h = hits.get(f"{cid}|{li}|{gi}")
                    if h:
                        g["func"], g["vfam"] = h["function"], h["vfam"]
                        g["vog_cat"] = h.get("category", "")
    else:
        print("  (VOGDB not present — genes shown as 'hypothetical'; run with the "
              "annotation DB to get function names)")

    produced, ann_rows = [], []
    for cid in sorted(by_cluster):
        loci = by_cluster[cid]
        if len(loci) < 2:
            continue
        assign_orthogroups(loci)
        consensus_functions(loci)
        for loc in loci:
            for g in loc["genes"]:
                g["category"] = ("gene of interest" if g["fam"]
                                 else categorize(g["func"], g.get("vog_cat", "")))
        modes = ["function", "conservation"] if args.color_by == "both" else [args.color_by]
        drew = False
        for mode in modes:
            suf = "" if mode == "function" else f"_{mode}"
            if draw_cluster(cid, loci, args.out_dir, color_by=mode, suffix=suf):
                drew = True
        if drew:
            produced.append(cid)
            print(f"  cluster_{cid}: synteny figure ({len(loci)} loci; {'/'.join(modes)})")
        for loc in loci:
            for g in loc["genes"]:
                ann_rows.append({
                    "cluster": cid, "organism": loc["organism"], "genome_id": loc["genome_id"],
                    "start": g["s"], "end": g["e"], "strand": "+" if g["st"] >= 0 else "-",
                    "is_family": g["fam"], "orthogroup": g["og"], "category": g["category"],
                    "vfam": g["vfam"], "function": "family homologue" if g["fam"] else g["func"],
                })

    heatmap = build_heatmap(by_cluster, args.out_dir)

    if ann_rows:
        keys = ["cluster", "organism", "genome_id", "start", "end", "strand",
                "is_family", "orthogroup", "category", "vfam", "function"]
        with (args.out_dir / "neighbour_gene_annotations.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(ann_rows)

    idx = ["<!DOCTYPE html><html><head><meta charset='utf-8'>",
           "<title>publication synteny figures</title>",
           "<style>body{font-family:sans-serif;margin:2em}img{max-width:100%;border:1px solid #ddd;"
           "border-radius:6px;margin:8px 0}h2{margin-top:1.5em}</style></head><body>",
           "<h1>Publication synteny figures</h1>",
           "<p>Each locus is anchored on the family gene (red) and strand-normalised. "
           "Genes are coloured by function (VOGDB VFAM); same colour = same function, named in "
           "the legend. Shaded links join homologous genes between adjacent rows. "
           "Per-gene functions are in <code>neighbour_gene_annotations.csv</code>; "
           "SVG/PDF versions sit alongside each PNG.</p>"]
    if heatmap:
        idx.append("<h2>Neighbourhood conservation</h2>"
                   "<p>Functional composition of each genome's neighbourhood "
                   "(<code>neighbourhood_conservation.csv</code>).</p>"
                   "<img src='neighbourhood_heatmap.png' alt='neighbourhood heatmap'>")
    for cid in produced:
        idx.append(f"<h2>Cluster {cid}</h2>")
        for suf, lab in (("", "coloured by function"), ("_conservation", "coloured by conservation")):
            if (args.out_dir / f"cluster_{cid}_synteny{suf}.png").exists():
                idx.append(f"<h3>{lab}</h3>"
                           f"<img src='cluster_{cid}_synteny{suf}.png' alt='cluster {cid} {lab}'>")
    idx.append("</body></html>")
    (args.out_dir / "index.html").write_text("\n".join(idx))
    print(f"Done. {len(produced)} synteny figures. Index: {args.out_dir / 'index.html'}")


if __name__ == "__main__":
    main()
