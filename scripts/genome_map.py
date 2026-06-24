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


def draw(genes: list, anchor, out_base, title: str, log=print, track_name="genome"):
    """Draw a linear genome map (PNG + SVG) coloured by functional category, your gene gold,
    the track labelled `track_name` (the phage/organism name). Tries pyGenomeViz, falls back
    to matplotlib. Returns the base path or None."""
    genes = [g for g in genes if g.get("start") is not None and g.get("end") is not None]
    if not genes:
        return None
    try:
        return _draw_pgv(genes, out_base, title, track_name, log)
    except Exception as e:
        log(f"  (pyGenomeViz unavailable: {e}; basic renderer)")
        try:
            return _draw_mpl(genes, anchor, out_base, title, track_name, log)
        except Exception as e2:
            log(f"  (genome map skipped: {e2})")
            return None


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
