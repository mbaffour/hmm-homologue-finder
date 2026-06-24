#!/usr/bin/env python3
"""genome_map.py — a linear genome map that marks a gene of interest among its neighbours.

Shared by the single-genome scan (scan_genome.py) and the database discovery run
(build_real_genbanks.py). Renders with **pyGenomeViz** (already a dependency) — genes drawn
as strand arrows labelled with the genome's own annotation / gp numbers — and falls back to
a basic matplotlib renderer if pyGenomeViz is unavailable. Genes are coloured by broad
FUNCTIONAL CATEGORY using the same scheme as the synteny figures (structural, packaging,
replication, transcription/regulation, lysis, …); the gene of interest is bold gold and
labelled by the phage/organism name from the record.
"""
from __future__ import annotations

from pathlib import Path

# Layout-role -> lane only (colour comes from functional category, below).
_LANE = {"overlap": 0.55}


def _scheme():
    """(CATEGORY_COLORS, FAMILY_COLOR, HYPO_COLOR, categorize) — reuse the synteny scheme
    so genome maps and synteny figures colour genes identically. Safe fallback."""
    try:
        from synteny_figure import CATEGORY_COLORS, FAMILY_COLOR, HYPO_COLOR, categorize
        return CATEGORY_COLORS, FAMILY_COLOR, HYPO_COLOR, categorize
    except Exception:
        return {}, "#ffd400", "#dde2e8", (lambda f, v="": "hypothetical")


def build_genes(anchor, called, flank_keys=None, label_of=None) -> list:
    """Build the gene list for draw() in ABSOLUTE genome coordinates. anchor = (a_s, a_e,
    a_strand) of the gene of interest; called = [(s, e, strand, meta)]. flank_keys = set of
    (s,e) to mark 'flank' (labelled). A called gene that IS the gene of interest (reciprocal
    overlap) is dropped, but a nested overprint partner (a long gene that merely spans it) is
    kept. Each gene gets a functional `category` (via the synteny categorizer on its product).
    Returns [{start, end, strand, role, label, category}]."""
    a_s, a_e, a_st = anchor
    flank_keys = flank_keys or set()
    label_of = label_of or (lambda m: m.get("gene") or m.get("product") or m.get("locus_tag") or "")
    categorize = _scheme()[3]
    genes = [{"start": a_s, "end": a_e, "strand": a_st, "role": "anchor",
              "label": "gene of interest", "category": "gene of interest"}]
    for (s, e, st, m) in called:
        ov = max(0, min(e, a_e) - max(s, a_s))
        if ov > 0.6 * (e - s) and ov > 0.6 * (a_e - a_s):
            continue                                   # the gene of interest's own call
        role = "overlap" if ov > 0 else ("flank" if (s, e) in flank_keys else "other")
        cat = m.get("category") or categorize(m.get("product") or m.get("gene") or "", m.get("vfam") or "")
        genes.append({"start": s, "end": e, "strand": st, "role": role,
                      "label": (label_of(m) if role in ("overlap", "flank") else ""),
                      "category": cat})
    return genes


MAP_TOOLS = ("pub", "pygenomeviz", "matplotlib", "easyfig", "auto")


def draw(genes: list, anchor, out_base, title: str, log=print, track_name="genome",
         tool="pub", genbank=None, labels=True):
    """Draw a linear genome map (PNG + SVG + PDF) coloured by functional category, your gene
    gold, the track labelled `track_name` (phage name; may be two lines name\\naccession).
    Genes are strand arrows (direction = strand); overlapping genes are packed onto separate
    lanes so nothing is hidden. `labels` toggles the gene-name labels. `tool`: 'pub' (default
    — the publication genome diagram), 'pygenomeviz', or 'easyfig' (needs a `genbank` + an
    installed Easyfig; falls back). Returns the base path or None."""
    genes = [g for g in genes if g.get("start") is not None and g.get("end") is not None]
    if not genes:
        return None
    tool = (tool or "pub").lower()
    if tool == "easyfig":
        try:
            return _draw_easyfig(genbank, out_base, title, log)
        except Exception as e:
            log(f"  (Easyfig unavailable: {e}; using the publication renderer)")
            tool = "pub"
    if tool == "pygenomeviz":
        try:
            return _draw_pgv(genes, out_base, title, track_name, log)
        except Exception as e:
            log(f"  (pyGenomeViz unavailable: {e}; using the publication renderer)")
            tool = "pub"
    try:
        return _draw_pub(genes, anchor, out_base, title, track_name, log, labels=labels)
    except Exception as e:
        log(f"  (genome map skipped: {e})")
        return None


def _nice_bar(span):
    import math
    raw = max(span / 5.0, 1)
    mag = 10 ** int(math.log10(raw))
    for m in (1, 2, 5, 10):
        if mag * m >= raw:
            return int(mag * m)
    return int(mag * 10)


def _gene_arrow(ax, s, e, st, ymid, h, head, color, ec, lw):
    from matplotlib.patches import Polygon
    L = e - s
    hd = min(head, L) if L > 0 else head
    if st >= 0:
        xh = e - hd
        pts = [(s, ymid - h), (xh, ymid - h), (xh, ymid - h * 1.7),
               (e, ymid), (xh, ymid + h * 1.7), (xh, ymid + h), (s, ymid + h)]
    else:
        xh = s + hd
        pts = [(e, ymid - h), (xh, ymid - h), (xh, ymid - h * 1.7),
               (s, ymid), (xh, ymid + h * 1.7), (xh, ymid + h), (e, ymid + h)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color, edgecolor=ec, lw=lw, zorder=3))


def _pack(genes):
    """Assign each gene to a lane so OVERLAPPING genes land on different lanes (greedy
    interval packing). Returns ({id(gene): lane_index}, n_lanes). Most genomes need 1-2
    lanes; a gene nested in/over another (e.g. an overprint partner) gets its own lane."""
    lane_end = []                              # lane_end[i] = end coord of last gene in lane i
    assign = {}
    for g in sorted(genes, key=lambda x: (x["start"], x["end"])):
        for i, end in enumerate(lane_end):
            if g["start"] >= end:              # fits after the last gene in lane i
                lane_end[i] = g["end"]
                assign[id(g)] = i
                break
        else:
            lane_end.append(g["end"])
            assign[id(g)] = len(lane_end) - 1
    return assign, max(len(lane_end), 1)


def _draw_pub(genes, anchor, out_base, title, track_name, log, labels=True):
    """Publication genome diagram (works for any genome): genes as strand arrows (direction
    = strand) coloured by functional category; the gene of interest is bold gold. Overlapping
    genes are PACKED onto separate lanes so nothing is hidden. Gene-name labels (toggle with
    `labels`) are stacked into collision-free rows with leader lines, using real label widths.
    Figure size + margins scale to the genome and the content. PNG (300 dpi) + SVG + PDF."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    CC, FAM, HYPO, _ = _scheme()
    lo = min(g["start"] for g in genes)
    hi = max(g["end"] for g in genes)
    span = max(hi - lo, 1)
    H, head = 0.36, max(span * 0.006, 1)
    lane_h, rowstep = 0.95, 0.6
    width = min(30.0, max(9.0, span / 2200.0))
    xlo, xhi = lo - span * 0.10, hi + span * 0.03
    xspan, axes_w_in = xhi - xlo, width * 0.94
    assign, nlanes = _pack(genes)
    top_y = 0.0                                # lane 0 at top; lanes stack downward
    bot_y = top_y - (nlanes - 1) * lane_h

    def _hw(text, fs):
        return 0.5 * len(text) * (fs * 0.62 / 72.0) * (xspan / axes_w_in) + xspan * 0.004
    # label rows (above the top lane), collision-free, only if requested
    place, label_rows = {}, 0
    if labels:
        occ = []
        for g in sorted([g for g in genes if g.get("label")
                         and g["role"] in ("anchor", "overlap", "flank")],
                        key=lambda x: (x["start"] + x["end"]) / 2):
            cx = (g["start"] + g["end"]) / 2
            hw = _hw(str(g["label"])[:30], 8.5 if g["role"] == "anchor" else 6.6)
            r = 0
            while r < len(occ) and (cx - hw) < occ[r]:
                r += 1
            if r == len(occ):
                occ.append(-1e18)
            occ[r] = cx + hw
            place[id(g)] = r
        label_rows = len(occ)
    height = max(4.2, 2.8 + 0.5 * nlanes + 0.42 * label_rows)
    fig, ax = plt.subplots(figsize=(width, height))
    # arrows on their packed lanes
    for g in sorted(genes, key=lambda x: x["start"]):
        st = 1 if int(g.get("strand", 1)) >= 0 else -1
        is_a = g["role"] == "anchor"
        _gene_arrow(ax, g["start"], g["end"], st, top_y - assign[id(g)] * lane_h,
                    H * (1.3 if is_a else 1.0), head,
                    FAM if is_a else CC.get(g.get("category", ""), HYPO),
                    "#1a1a1a" if is_a else "#333333", 1.8 if is_a else 0.5)
    # labels above, with leader lines down to each gene's lane
    lbase = top_y + H + 0.25
    if labels:
        for g, r in ((g, place[id(g)]) for g in genes if id(g) in place):
            cx = (g["start"] + g["end"]) / 2
            is_a = g["role"] == "anchor"
            gy = top_y - assign[id(g)] * lane_h + H + 0.05
            yt = lbase + rowstep * (r + 1)
            ax.plot([cx, cx], [gy, yt], color="#bbb", lw=0.4, zorder=1)
            ax.text(cx, yt, str(g["label"])[:30], ha="center", va="bottom",
                    fontsize=(8.5 if is_a else 6.6), fontweight=("bold" if is_a else "normal"),
                    color=("#806000" if is_a else "#1a2230"), zorder=4)
    # scale bar below the lowest lane
    bar = _nice_bar(span)
    ybar = bot_y - (H + 0.7)
    ax.plot([lo, lo + bar], [ybar, ybar], color="#333", lw=2.0)
    ax.text(lo + bar / 2.0, ybar - 0.12, (f"{bar // 1000} kb" if bar >= 1000 else f"{bar} bp"),
            ha="center", va="top", fontsize=7.5)
    ytop = (lbase + rowstep * (label_rows + 1)) if labels else (top_y + H + 0.4)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ybar - 0.5, ytop + 0.3)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    # title/caption/legend positioned in ABSOLUTE inches from the edges, so spacing is
    # consistent whatever the figure height (which varies with lanes/label rows).
    nm = str(track_name).split("\n")
    fig.subplots_adjust(top=1 - 1.25 / height, bottom=1.05 / height, left=0.04, right=0.99)
    fig.suptitle(nm[0], fontsize=12, fontweight="bold", y=1 - 0.38 / height)
    if len(nm) > 1:
        fig.text(0.5, 1 - 0.80 / height, nm[1], ha="center", fontsize=9, color="#555")
    fig.text(0.5, 0.72 / height, title, ha="center", fontsize=8.5, color="#444")
    fig.legend(handles=_legend_handles(genes, CC, FAM, HYPO), loc="lower center",
               bbox_to_anchor=(0.5, 0.06 / height), ncol=5, fontsize=7.5, frameon=False)
    out_base = Path(out_base)
    try:
        fig.savefig(f"{out_base}.png", dpi=300)
        fig.savefig(f"{out_base}.svg")
        fig.savefig(f"{out_base}.pdf")
    finally:
        plt.close(fig)
    log(f"  genome map -> {out_base.name}.png / .svg / .pdf")
    return out_base


def write_locus_genbank(genes, contig_seq, organism, accession, out_path):
    """Write a GenBank of the locus (the window spanning `genes`, with the contig sequence),
    CDS features carrying the annotation/gp-number labels and the gene of interest marked
    /gene=gene_of_interest. This is the tool-agnostic deliverable — open it in **Easyfig**,
    Artemis, clinker, pyGenomeViz, etc. Returns the path (or None if the sequence is absent)."""
    if not contig_seq:
        return None
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    wlo = min(g["start"] for g in genes)
    whi = max(g["end"] for g in genes)
    sub = contig_seq[max(0, wlo - 1):whi]
    acc = (accession or "locus")
    rec = SeqRecord(Seq(sub), id=acc[:16], name=acc[:16],
                    description=f"{organism} ({acc}) — gene-of-interest neighbourhood",
                    annotations={"molecule_type": "DNA", "topology": "linear",
                                 "organism": organism or acc})
    feats = []
    for g in genes:
        if g["role"] == "anchor":
            q = {"gene": ["gene_of_interest"], "product": ["gene of interest (HMM hit)"]}
        else:
            q = {"product": [g.get("label") or "CDS"]}
        feats.append(SeqFeature(FeatureLocation(max(0, g["start"] - wlo), max(1, g["end"] - wlo),
                     strand=int(g.get("strand", 1))), type="CDS", qualifiers=q))
    rec.features = sorted(feats, key=lambda f: int(f.location.start))
    SeqIO.write(rec, str(out_path), "genbank")
    return out_path


def _draw_easyfig(genbank, out_base, title, log):
    """Render with Easyfig (github.com/mjsull/Easyfig) — a Python-2 standalone script, NOT a
    pip/conda package. Enable by installing it and setting $EASYFIG_PY to Easyfig.py (and
    $EASYFIG_PYTHON to a python2 if needed). Raises (-> caller falls back) when unavailable."""
    import os
    import shutil
    import subprocess
    if not genbank or not Path(genbank).exists():
        raise RuntimeError("Easyfig needs a locus GenBank (none provided)")
    ef = os.environ.get("EASYFIG_PY") or shutil.which("Easyfig.py") or shutil.which("Easyfig")
    if not ef:
        raise RuntimeError("Easyfig not found — install github.com/mjsull/Easyfig and set "
                           "$EASYFIG_PY to Easyfig.py (the locus GenBank is written for it)")
    py = (os.environ.get("EASYFIG_PYTHON") or shutil.which("python2")
          or shutil.which("python2.7") or shutil.which("python"))
    png = f"{out_base}.png"
    subprocess.run([py, ef, "-f", "-i", str(genbank), "-o", png, "-legend", "single",
                    "-leg_name", "product"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not Path(png).exists() or Path(png).stat().st_size < 1000:
        raise RuntimeError("Easyfig produced no usable output")
    log(f"  genome map (Easyfig) -> {Path(out_base).name}.png")
    return Path(out_base)


def _legend_handles(genes, CC, FAM, HYPO):
    from matplotlib.patches import Patch
    present = []
    for g in genes:
        c = g.get("category")
        if g["role"] != "anchor" and c in CC and c not in present:
            present.append(c)
    handles = [Patch(facecolor=FAM, edgecolor="#1a1a1a", linewidth=1.2, label="gene of interest")]
    handles += [Patch(facecolor=CC[c], edgecolor="#33373d", label=c) for c in present]
    handles += [Patch(facecolor=HYPO, edgecolor="#33373d", label="hypothetical / other")]
    return handles


def _color_of(g, CC, FAM, HYPO):
    return FAM if g["role"] == "anchor" else CC.get(g.get("category", ""), HYPO)


def _draw_pgv(genes, out_base, title, track_name, log):
    from pygenomeviz import GenomeViz
    import matplotlib.pyplot as plt
    CC, FAM, HYPO, _ = _scheme()
    wlo = min(g["start"] for g in genes)
    whi = max(g["end"] for g in genes)
    L = max(whi - wlo, 1)
    gv = GenomeViz(fig_width=min(22.0, max(9.0, L / 3500.0)), fig_track_height=0.6, show_axis=True)
    track = gv.add_feature_track((str(track_name) or "genome")[:45], L)
    seg = track.get_segment()
    for g in genes:
        seg.add_feature(g["start"] - wlo, g["end"] - wlo, int(g.get("strand", 1)),
                        plotstyle="arrow", label=g.get("label", "") or "",
                        fc=_color_of(g, CC, FAM, HYPO),
                        ec=("#1a1a1a" if g["role"] == "anchor" else "#33373d"),
                        lw=(1.4 if g["role"] == "anchor" else 0.4),
                        text_kws={"rotation": 30, "size": 7, "vpos": "top", "hpos": "left"})
    fig = gv.plotfig()
    handles = _legend_handles(genes, CC, FAM, HYPO)
    fig.legend(handles=handles, loc="lower center", ncol=min(4, len(handles)), fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, -0.5))
    fig.text(0.5, -0.62, title, ha="center", va="top", fontsize=10)
    out_base = Path(out_base)
    try:
        fig.savefig(f"{out_base}.png", dpi=200, bbox_inches="tight")
        fig.savefig(f"{out_base}.svg", bbox_inches="tight")
    finally:
        plt.close(fig)
    log(f"  genome map -> {out_base.name}.png / .svg")
    return out_base


def _draw_mpl(genes, anchor, out_base, title, track_name, log):
    """Fallback renderer (no pyGenomeViz): arrows along the position axis, relative to the
    gene of interest, coloured by functional category."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrow
    CC, FAM, HYPO, _ = _scheme()
    a_s = anchor[0]
    rel = [{**g, "s": g["start"] - a_s, "e": g["end"] - a_s} for g in genes]
    lo, hi = min(g["s"] for g in rel), max(g["e"] for g in rel)
    span = max(hi - lo, 1)
    fig, ax = plt.subplots(figsize=(min(20, max(9, span / 3500)), 3.2))
    for g in rel:
        y = _LANE.get(g["role"], 0.0)
        right = int(g.get("strand", 1)) >= 0
        x0, dx = (g["s"], g["e"] - g["s"]) if right else (g["e"], g["s"] - g["e"])
        ax.add_patch(FancyArrow(x0, y, dx or 1, 0, width=0.16, length_includes_head=True,
                                head_width=0.3, head_length=max(span * 0.005, 1),
                                color=_color_of(g, CC, FAM, HYPO),
                                ec=("#1a1a1a" if g["role"] == "anchor" else "#33373d"),
                                lw=(1.4 if g["role"] == "anchor" else 0.4), zorder=3))
        if g.get("label") and g["role"] in ("anchor", "overlap", "flank"):
            ax.text((g["s"] + g["e"]) / 2, y + 0.27, str(g["label"])[:24], ha="center",
                    va="bottom", fontsize=6.5, rotation=28, color="#1a2230")
    ax.axhline(0, color="#c9ced6", lw=0.6, zorder=0)
    ax.set_xlim(lo - span * 0.03, hi + span * 0.03)
    ax.set_ylim(-0.7, 1.6)
    ax.set_yticks([])
    ax.set_xlabel(f"position relative to your gene (bp) — {track_name}", fontsize=9)
    ax.set_title(title, fontsize=9)
    ax.legend(handles=_legend_handles(genes, CC, FAM, HYPO), loc="upper center",
              bbox_to_anchor=(0.5, -0.18), ncol=4, fontsize=7, frameon=False)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    out_base = Path(out_base)
    try:
        fig.savefig(f"{out_base}.png", dpi=200, bbox_inches="tight")
        fig.savefig(f"{out_base}.svg", bbox_inches="tight")
    finally:
        plt.close(fig)
    log(f"  genome map -> {out_base.name}.png / .svg")
    return out_base
