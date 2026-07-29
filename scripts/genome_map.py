#!/usr/bin/env python3
"""genome_map.py — a linear genome map that marks a gene of interest among its neighbours.

Shared by the single-genome scan (scan_genome.py) and the database discovery run
(build_real_genbanks.py). Renders by default with **DNA Features Viewer** (the Edinburgh
Genome Foundry library, tool key 'dfv') — genes drawn as clean strand arrows labelled with
the genome's own annotation / gp numbers, overlapping genes auto-stacked onto their own
level, and a real genome-coordinate axis spanning the whole locus — and falls back to the
built-in matplotlib 'pub' renderer (arrow-direction strand, packed lanes, a full-length
coordinate ruler) when DNA Features Viewer is unavailable. 'pygenomeviz' and 'easyfig' are
also selectable. Genes are coloured by broad FUNCTIONAL CATEGORY using the same scheme as the
synteny figures (structural, packaging, replication, transcription/regulation, lysis, …); the
gene of interest is bold gold and labelled by the phage/organism name from the record.
"""
from __future__ import annotations

from pathlib import Path

# Layout-role -> lane only (colour comes from functional category, below).
_LANE = {"overlap": 0.55}

# ONE definition of "overprint" for this whole file. Phage genes routinely share a base or
# two (the ATGA start/stop overlap is the commonest arrangement in a phage operon), so a
# bare "overlaps by >=1 bp" is an operon neighbour, not an overprint. Both the row-packing
# (_dfv_tolerant_levels) and the legend (_legend_handles) key off this same threshold, so a
# gene that is stacked onto its own row because it genuinely overprints the gene of interest
# is exactly the gene the legend names — the figure cannot say two different things.
OVERPRINT_MIN_OVERLAP_BP = 60


def _scheme(palette="default"):
    """(CATEGORY_COLORS, FAMILY_COLOR, HYPO_COLOR, categorize) — reuse the synteny scheme
    so genome maps and synteny figures colour genes identically. `palette` selects the colour
    set: 'default' or 'colorblind' (Paul Tol muted, colour-blind-safe). Safe fallback."""
    try:
        from synteny_figure import PALETTES, FAMILY_COLOR, HYPO_COLOR, categorize
        return PALETTES.get(palette, PALETTES["default"]), FAMILY_COLOR, HYPO_COLOR, categorize
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
        # store the name for EVERY gene (so the whole-genome locus GenBank keeps gp/locus names,
        # not /product="CDS"); the renderer decides which labels to DISPLAY via its density budget.
        genes.append({"start": s, "end": e, "strand": st, "role": role,
                      "label": label_of(m), "category": cat})
    return genes


MAP_TOOLS = ("dfv", "pub", "pygenomeviz", "matplotlib", "easyfig", "auto")


def draw(genes: list, anchor, out_base, title: str, log=print, track_name="genome",
         tool="dfv", genbank=None, labels=True, palette="default",
         functional_labels=False, module_brackets=False, marks=None):
    """Draw a linear genome map (PNG + SVG + PDF) coloured by functional category, your gene
    gold, the track labelled `track_name` (phage name; may be two lines name\\naccession).
    Genes are strand arrows (direction = strand); overlapping genes are stacked onto separate
    levels so nothing is hidden. `labels` toggles the gene-name labels. `tool`: 'dfv' (default
    — DNA Features Viewer, the cleanest publication renderer), 'pub' (the built-in matplotlib
    diagram, always available), 'pygenomeviz', or 'easyfig' (needs a `genbank` + an installed
    Easyfig; falls back). Any renderer that is unavailable falls back to 'pub'.
    `palette`: 'default' or 'colorblind' (Paul Tol muted). `functional_labels`: also tag the
    gene of interest + its overlap partner with their functional category. `module_brackets`:
    bracket contiguous same-category runs with the module name (dfv only).
    `marks`: optional [(nt_pos, label, colour)] drawn as labelled vertical ticks at genome
    coordinates — e.g. one per PREMATURE STOP on an interrupted/overprinted locus, which is
    the single most informative element of that figure and was previously inexpressible.
    Rendered by the 'dfv' and 'pub' renderers; the default None is a strict no-op, so every
    existing caller's output is unchanged. Returns the base path or None."""
    genes = [g for g in genes if g.get("start") is not None and g.get("end") is not None]
    if not genes:
        return None
    tool = (tool or "dfv").lower()
    if tool == "auto":
        tool = "dfv"
    if tool == "easyfig":
        try:
            return _draw_easyfig(genbank, out_base, title, log)
        except Exception as e:
            log(f"  (Easyfig unavailable: {e}; using DNA Features Viewer)")
            tool = "dfv"
    if tool == "pygenomeviz":
        try:
            return _draw_pgv(genes, out_base, title, track_name, log)
        except Exception as e:
            log(f"  (pyGenomeViz unavailable: {e}; using DNA Features Viewer)")
            tool = "dfv"
    if tool == "dfv":
        try:
            return _draw_dfv(genes, anchor, out_base, title, track_name, log, labels=labels,
                             palette=palette, functional_labels=functional_labels,
                             module_brackets=module_brackets, marks=marks)
        except Exception as e:
            log(f"  (DNA Features Viewer unavailable: {e}; using the built-in renderer)")
            tool = "pub"
    try:
        return _draw_pub(genes, anchor, out_base, title, track_name, log, labels=labels,
                         marks=marks)
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


def _dfv_tolerant_levels(features, tol):
    """Greedy row-packing for DNA Features Viewer that TOLERATES small overlaps: two genes
    share a row unless they overlap by more than `tol` bp. Phage genes very often share a
    base or two (overlapping start/stop codons); DFV's default stacks every such pair onto a
    new row, which turns a dense genome into a cluttered staircase. Collapsing sub-`tol`
    overlaps onto the baseline keeps the map clean, while a genuinely nested / overprinted
    gene (overlap >> tol, e.g. a gene inside an antisense RNA polymerase) still gets its own
    row so nothing is hidden. Returns {feature: row_level} (0 = baseline)."""
    feats = sorted(features, key=lambda f: (f.start, f.end))
    lane_end, levels = [], {}
    for f in feats:
        for i in range(len(lane_end)):
            if f.start >= lane_end[i] - tol:        # only a >tol overlap forces a new row
                lane_end[i] = max(lane_end[i], f.end)
                levels[f] = i
                break
        else:
            lane_end.append(f.end)
            levels[f] = len(lane_end) - 1
    return levels


def _module_runs(genes, CC):
    """Contiguous runs (>=2 genes) sharing a functional category — the phage's 'modules'
    (structural, replication, …). Returns [{cat, start, end, n}] in genome coords, skipping
    the anchor and hypothetical/unknown genes."""
    ordered = sorted([g for g in genes if g["role"] != "anchor" and g.get("category") in CC
                      and "hypothetical" not in (g.get("category") or "").lower()],
                     key=lambda g: g["start"])
    runs, cur = [], None
    for g in ordered:
        c = g["category"]
        if cur and cur["cat"] == c and g["start"] <= cur["end"] + (g["end"] - g["start"]) * 3:
            cur["end"] = max(cur["end"], g["end"])
            cur["n"] += 1
        else:
            cur = {"cat": c, "start": g["start"], "end": g["end"], "n": 1}
            runs.append(cur)
    return [r for r in runs if r["n"] >= 2]


def _draw_module_brackets(axes_list, runs):
    """Draw a labelled bracket above each functional-module run, on whichever wrapped line it
    falls (genoPlotR/Phamerator modular-organisation convention)."""
    for ax in axes_list:
        x0, x1 = sorted(ax.get_xlim())
        ymin, ymax = ax.get_ylim()
        h = (ymax - ymin) * 0.04
        drew = False
        for r in runs:
            rs, re = max(r["start"], x0), min(r["end"], x1)
            if re <= rs:
                continue
            yb = ymax + h
            ax.plot([rs, rs, re, re], [yb, yb + h, yb + h, yb],
                    color="#666", lw=0.8, clip_on=False, zorder=5)
            ax.text((rs + re) / 2, yb + h * 1.3, f"{r['cat'].split(' / ')[0]} module",
                    ha="center", va="bottom", fontsize=6.5, color="#555",
                    style="italic", clip_on=False, zorder=5)
            drew = True
        if drew:
            ax.set_ylim(ymin, ymax + h * 4)            # headroom for the brackets


def _draw_marks(axes_list, marks):
    """Draw each `(nt_pos, label, colour)` as a dashed vertical tick at its GENOME coordinate,
    labelled along the tick, on whichever wrapped line contains that coordinate. Built for the
    premature stops of an interrupted/overprinted locus: the stop is a single base, so a gene
    arrow cannot show it and only a positional tick can.

    Strict no-op for a None/empty list — that is what keeps every pre-existing caller's figure
    byte-identical. Individually tolerant of malformed entries (skips them) and it never moves
    the axes limits: `axvline` + a blended-transform label add no data extent, so the gene rows
    laid out above stay exactly where the renderer put them."""
    if not marks:
        return
    for ax in axes_list:
        x0, x1 = sorted(ax.get_xlim())
        for m in marks:
            try:
                pos = float(m[0])
                lab = str(m[1]) if len(m) > 1 and m[1] else ""
                col = (m[2] if len(m) > 2 and m[2] else "#c0392b")
            except Exception:
                continue                                    # malformed mark -> skip, never raise
            if not (x0 <= pos <= x1):
                continue                                    # falls on another wrapped line
            ax.axvline(pos, color=col, lw=1.1, ls=(0, (3, 2)), zorder=25)
            if lab:
                # x in DATA coords, y in AXES fraction -> the label always sits at the top of
                # the panel whatever y-range the renderer ended up with.
                ax.text(pos, 0.99, lab, transform=ax.get_xaxis_transform(), rotation=90,
                        ha="right", va="top", fontsize=6, color=col, zorder=26)


def _draw_dfv(genes, anchor, out_base, title, track_name, log, labels=True,
              palette="default", functional_labels=False, module_brackets=False,
              marks=None):
    """Publication genome map via **DNA Features Viewer** (Edinburgh Genome Foundry) — the
    default, cleanest renderer. Clean strand arrows; OVERLAPPING genes auto-stacked onto their
    own level (an overprint partner never hides the gene of interest); labels de-overlapped with
    leader lines; a real genome-coordinate axis. Genes coloured by functional category (gene of
    interest bold gold). `palette` 'default'/'colorblind'. Label policy is density-aware: the
    gene of interest + overlapping genes are ALWAYS labelled; other genes are labelled closest-
    to-the-anchor-first up to a width budget. `functional_labels` adds a category tag to the
    gene of interest + its overlap partner. Big genomes (>40 genes) WRAP onto multiple lines so
    every gene has room. `module_brackets` brackets contiguous same-category runs. PNG/SVG/PDF."""
    import math
    from dna_features_viewer import GraphicFeature, GraphicRecord
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    CC, FAM, HYPO, _ = _scheme(palette)
    wlo = min(g["start"] for g in genes)
    whi = max(g["end"] for g in genes)
    span = max(whi - wlo, 1)
    ngenes = len(genes)
    multiline = ngenes > 40
    fig_w = min(30.0, max(11.0, span / 1800.0))

    # density-aware label budget: anchor+overlap always; other genes closest-to-anchor first,
    # up to ~5 labels/inch (no cap when wrapped onto multiple lines, which gives the room).
    a_centre = (anchor[0] + anchor[1]) / 2 if anchor else (wlo + whi) / 2
    budget = (10 ** 9 if multiline else int(fig_w * 5)) if labels else 0
    flanks = sorted([g for g in genes if g["role"] not in ("anchor", "overlap") and g.get("label")],
                    key=lambda g: abs((g["start"] + g["end"]) / 2 - a_centre))
    keep = {id(g) for g in flanks[:max(0, budget)]}

    def _tag(g):
        c = g.get("category") or ""
        return f"\n[{c.split(' / ')[0]}]" if (functional_labels and c in CC) else ""

    def _lab(g):
        # The gene of interest is identified by its bold-gold colour + the legend, so it
        # gets NO inline text label: on a short anchor arrow DFV centres the wide
        # "gene of interest" box on the feature and the (always-raised, see below) gold
        # arrow then bisects it ("gen…rest"). functional_labels can still tag it.
        if g["role"] == "anchor":
            return (_tag(g).strip() or None)
        if g["role"] == "overlap":
            return (g.get("label") or "overlapping gene") + _tag(g)
        if labels and id(g) in keep and g.get("label"):
            return g["label"]
        return None

    feats = []
    for g in genes:
        is_a = g["role"] == "anchor"
        feats.append(GraphicFeature(
            start=g["start"], end=g["end"], strand=int(g.get("strand", 1)),
            color=(FAM if is_a else CC.get(g.get("category", ""), HYPO)),
            linecolor=("#1a1a1a" if is_a else "#33373d"),
            linewidth=(2.2 if is_a else 0.6),
            label=_lab(g),
            fontdict=({"fontsize": 9, "fontweight": "bold", "color": "#6b5300"} if is_a
                      else {"fontsize": (6 if multiline else 7), "color": "#1a2230"})))
    rec = GraphicRecord(sequence_length=span + 1, features=feats, first_index=wlo,
                        labels_spacing=14)
    # Row-packing for GENE ARROWS only: collapse few-bp start/stop overlaps onto the baseline
    # while stacking a real overprint onto its own row. DFV calls compute_features_levels 3×
    # (arrows + 2× label boxes); the 60 bp tolerance must NOT touch the label boxes (doing so
    # forced nearby labels onto one row — the overlap bug), so we route label / elevate-base
    # pseudo-features (they carry 'text'/'is_base' in .data) to DFV's exact packer.
    import dna_features_viewer as _dfvmod
    from dna_features_viewer.compute_features_levels import compute_features_levels as _cfl_orig
    _pg = _dfvmod.GraphicRecord.plot.__globals__
    _orig_cfl = _pg.get("compute_features_levels")

    def _patched_cfl(fs):
        if fs and any(("text" in getattr(f, "data", {})) or ("is_base" in getattr(f, "data", {}))
                      for f in fs):
            return _cfl_orig(fs)
        return _dfv_tolerant_levels(fs, OVERPRINT_MIN_OVERLAP_BP)

    _pg["compute_features_levels"] = _patched_cfl
    try:
        if multiline:
            n_lines = max(2, math.ceil(ngenes / 35))   # ~35 genes/line
            nucl = math.ceil((span + 1) / n_lines)
            fig, axes = rec.plot_on_multiple_lines(
                nucl_per_line=nucl, figure_width=min(22.0, max(12.0, nucl / 1800.0)),
                elevate_outline_annotations=True, annotate_inline=True,
                max_label_length=46, max_line_length=20)
            axes_list = list(axes) if hasattr(axes, "__iter__") else [axes]
        else:
            ax = rec.plot(figure_width=fig_w, elevate_outline_annotations=True,
                          annotate_inline=True, max_label_length=60, max_line_length=24)[0]
            axes_list = [ax]
            fig = ax.figure
    finally:
        if _orig_cfl is not None:
            _pg["compute_features_levels"] = _orig_cfl

    # Keep the gold gene-of-interest arrow VISIBLE in the raster: DFV draws the (wide) inline
    # "gene of interest" label box on top of the tiny anchor arrow, which fully occludes it in
    # the PNG/PDF (the arrow survives only in the SVG). Raise every gold family-colour arrow
    # above the label boxes so the paper's subject gene is never hidden.
    import matplotlib.colors as _mc
    _gold = _mc.to_hex(FAM)
    for _ax in axes_list:
        for _p in getattr(_ax, "patches", []):
            try:
                if _mc.to_hex(_p.get_facecolor()) == _gold:
                    _p.set_zorder(30)
            except Exception:
                pass

    if module_brackets:
        _draw_module_brackets(axes_list, _module_runs(genes, CC))
    # after the brackets (which resize the y-axis) so the tick labels land on the final top
    _draw_marks(axes_list, marks)

    nm = str(track_name).split("\n")
    header = nm[0] + (f"\n{nm[1]}" if len(nm) > 1 else "")
    handles = _legend_handles(genes, CC, FAM, HYPO)
    if multiline:
        # reserve a top strip (title, clear of line-1's module brackets) and a bottom strip
        # (caption + legend, clear of the last line's ruler); add inter-line spacing so a line's
        # brackets/labels don't touch the ruler of the line above. Grow the figure so these
        # reserved strips don't squeeze the gene rows.
        fw, fh = fig.get_size_inches()
        fh2 = fh + 1.9
        fig.set_size_inches(fw, fh2)
        fig.subplots_adjust(top=1 - 0.95 / fh2, bottom=1.30 / fh2, hspace=1.1)
        fig.suptitle(header, fontsize=12, fontweight="bold", y=1 - 0.32 / fh2)
        fig.text(0.5, 0.82 / fh2, title, ha="center", va="center", fontsize=8.5, color="#444")
        fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=7.5, frameon=False,
                   bbox_to_anchor=(0.5, 0.20 / fh2))
    else:
        ax = axes_list[0]
        ax.set_title(header, fontsize=12, fontweight="bold", pad=12)
        fw, fh = fig.get_size_inches()
        y0 = ax.get_position().y0                  # axes bottom (figure fraction)
        # caption + legend a fixed INCH gap below the axes so they always clear the ruler ticks
        fig.text(0.5, y0 - 0.55 / fh, title, ha="center", va="top", fontsize=8.5, color="#444")
        fig.legend(handles=handles, loc="upper center", ncol=5, fontsize=7.5, frameon=False,
                   bbox_to_anchor=(0.5, y0 - 0.85 / fh))
    out_base = Path(out_base)
    try:
        fig.savefig(f"{out_base}.png", dpi=300, bbox_inches="tight")
        fig.savefig(f"{out_base}.svg", bbox_inches="tight")
        fig.savefig(f"{out_base}.pdf", bbox_inches="tight")
    finally:
        plt.close(fig)
    log(f"  genome map -> {out_base.name}.png / .svg / .pdf")
    return out_base


def _draw_pub(genes, anchor, out_base, title, track_name, log, labels=True, marks=None):
    """Publication genome diagram (works for any genome): genes as strand arrows (direction
    = strand) coloured by functional category; the gene of interest is bold gold. Overlapping
    genes are PACKED onto separate lanes so nothing is hidden. Gene-name labels (toggle with
    `labels`) are stacked into rows using the ACTUAL rendered label widths so they never
    overlap, with leader lines. A genome-coordinate RULER (ticks + kb) runs the whole length.
    Figure size + margins scale to the genome and the content. PNG (300 dpi) + SVG + PDF."""
    import math
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
    assign, nlanes = _pack(genes)
    top_y = 0.0
    bot_y = top_y - (nlanes - 1) * lane_h
    lbase = top_y + H + 0.25
    # provisional figure (final WIDTH; height set after labels are measured). Fix the x
    # margins now so measured label widths stay valid after the later height resize.
    fig, ax = plt.subplots(figsize=(width, 6.0))
    fig.subplots_adjust(left=0.04, right=0.99)
    ax.set_xlim(xlo, xhi)
    # arrows on their packed lanes
    for g in sorted(genes, key=lambda x: x["start"]):
        st = 1 if int(g.get("strand", 1)) >= 0 else -1
        is_a = g["role"] == "anchor"
        _gene_arrow(ax, g["start"], g["end"], st, top_y - assign[id(g)] * lane_h,
                    H * (1.3 if is_a else 1.0), head,
                    FAM if is_a else CC.get(g.get("category", ""), HYPO),
                    "#1a1a1a" if is_a else "#333333", 1.8 if is_a else 0.5)
    # labels: create, MEASURE the real rendered width, stack into non-overlapping rows
    label_rows = 0
    if labels:
        items = []
        for g in sorted([g for g in genes if g.get("label")
                         and g["role"] in ("anchor", "overlap", "flank")],
                        key=lambda x: (x["start"] + x["end"]) / 2):
            cx = (g["start"] + g["end"]) / 2
            is_a = g["role"] == "anchor"
            t = ax.text(cx, lbase + rowstep, str(g["label"])[:30], ha="center", va="bottom",
                        fontsize=(8.5 if is_a else 6.6), fontweight=("bold" if is_a else "normal"),
                        color=("#806000" if is_a else "#1a2230"), zorder=4)
            items.append((g, cx, is_a, t))
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        inv = ax.transData.inverted()
        occ = []
        for g, cx, is_a, t in items:
            bb = t.get_window_extent(renderer=rend)
            x0 = inv.transform((bb.x0, bb.y0))[0]
            x1 = inv.transform((bb.x1, bb.y0))[0]
            hw = abs(x1 - x0) / 2 + span * 0.004
            r = 0
            while r < len(occ) and (cx - hw) < occ[r]:
                r += 1
            if r == len(occ):
                occ.append(-1e18)
            occ[r] = cx + hw
            yt = lbase + rowstep * (r + 1)
            t.set_position((cx, yt))
            gy = top_y - assign[id(g)] * lane_h + H + 0.05
            ax.plot([cx, cx], [gy, yt], color="#bbb", lw=0.4, zorder=1)
        label_rows = len(occ)
    # genome-coordinate ruler (full length) below the lowest lane
    tick = _nice_bar(span)
    yruler = bot_y - (H + 0.85)
    ax.plot([lo, hi], [yruler, yruler], color="#444", lw=1.2, zorder=2)
    x = math.ceil(lo / tick) * tick
    while x <= hi + 1:
        ax.plot([x, x], [yruler, yruler - 0.13], color="#444", lw=1.0, zorder=2)
        ax.text(x, yruler - 0.2, (f"{x // 1000:g} kb" if x >= 1000 else f"{int(x)} bp"),
                ha="center", va="top", fontsize=6.8, color="#333")
        x += tick
    # finalise size now that label rows are known; keep x-margins fixed (widths still valid)
    height = max(4.2, 3.0 + 0.5 * nlanes + 0.42 * label_rows)
    fig.set_size_inches(width, height)
    ytop = (lbase + rowstep * (label_rows + 1)) if labels else (top_y + H + 0.4)
    ax.set_ylim(yruler - 0.55, ytop + 0.3)
    _draw_marks([ax], marks)          # after the final ylim, so tick labels sit at the top
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    nm = str(track_name).split("\n")
    fig.subplots_adjust(top=1 - 1.25 / height, bottom=1.0 / height, left=0.04, right=0.99)
    fig.suptitle(nm[0], fontsize=12, fontweight="bold", y=1 - 0.38 / height)
    if len(nm) > 1:
        fig.text(0.5, 1 - 0.80 / height, nm[1], ha="center", fontsize=9, color="#555")
    fig.text(0.5, 0.62 / height, title, ha="center", fontsize=8.5, color="#444")
    fig.legend(handles=_legend_handles(genes, CC, FAM, HYPO), loc="lower center",
               bbox_to_anchor=(0.5, 0.05 / height), ncol=5, fontsize=7.5, frameon=False)
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
        # g["start"]/g["end"] are 1-based INCLUSIVE; the sub-sequence's first base (1-based wlo)
        # is index 0. Biopython FeatureLocation is 0-based half-open, so a 1-based-inclusive
        # [s,e] maps to [s-wlo, e-wlo+1] — the +1 on the end is required or every CDS prints 1 bp
        # short at its 3' end.
        feats.append(SeqFeature(FeatureLocation(max(0, g["start"] - wlo), max(1, g["end"] - wlo + 1),
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


def _overprint_partners(genes):
    """The genes that genuinely ANTISENSE-OVERPRINT the gene of interest: opposite strand AND
    overlapping it by at least OVERPRINT_MIN_OVERLAP_BP.

    Both conditions are load-bearing, and `role == 'overlap'` (which is >=1 bp, either strand)
    satisfies neither:
      * SAME-STRAND overlaps are not overprinting at all — they are ordinary translational
        coupling (ATGA / TAATG), the commonest gene arrangement in a phage operon.
      * A few-bp opposite-strand overlap of two convergent genes is a gene boundary, not a
        nested gene; only an overlap that is a real fraction of a gene can host one.
    Returns [] when the anchor is missing or carries no coordinates, so a caller cannot
    conjure an overprint out of a gene list that never described one.
    """
    anchor = next((g for g in genes if g.get("role") == "anchor"), None)
    if not anchor:
        return []
    try:
        a_s, a_e = int(anchor["start"]), int(anchor["end"])
        a_st = int(anchor.get("strand", 1) or 0)
    except Exception:
        return []
    if not a_st:                       # unstranded anchor: "antisense" is undefined
        return []
    out = []
    for g in genes:
        if g is anchor or g.get("role") == "anchor":
            continue
        try:
            s, e, st = int(g["start"]), int(g["end"]), int(g.get("strand", 1) or 0)
        except Exception:
            continue
        if not st or (st > 0) == (a_st > 0):
            continue                                       # same strand (or unstranded)
        if min(e, a_e) - max(s, a_s) + 1 >= OVERPRINT_MIN_OVERLAP_BP:
            out.append(g)
    return out


def _legend_handles(genes, CC, FAM, HYPO):
    """Legend patches for the categories actually present, each annotated with a COUNT
    (e.g. 'structural (3)') so the functional composition is readable without counting arrows."""
    from collections import Counter
    from matplotlib.patches import Patch
    counts = Counter(g.get("category") for g in genes if g["role"] != "anchor")
    present = [c for c in CC if counts.get(c)]
    n_hypo = sum(v for c, v in counts.items() if c not in CC)
    handles = [Patch(facecolor=FAM, edgecolor="#1a1a1a", linewidth=1.2, label="gene of interest")]
    handles += [Patch(facecolor=CC[c], edgecolor="#33373d", label=f"{c} ({counts[c]})") for c in present]
    if n_hypo:
        handles += [Patch(facecolor=HYPO, edgecolor="#33373d", label=f"hypothetical / other ({n_hypo})")]
    # A gene that really is nested antisense inside the gene of interest was never named in
    # the legend, so a reader saw an arrow straddling the gold one with nothing to say what
    # that means — for this project it is the headline biology. The swatch is deliberately
    # UNFILLED: an overprinting partner is coloured by its own function (transcription,
    # replication, …) like every other gene, so claiming a colour here would contradict the
    # category entries above. No count suffix (unlike the category entries): the role is a
    # RELATIONSHIP to the gene of interest, and the arrows carrying it are already the ones
    # straddling the gold one.
    # The test is _overprint_partners, NOT `role == 'overlap'`: role is >=1 bp on EITHER
    # strand, so keying off it stamped "overprinting partner (antisense)" on a 4 bp
    # same-strand ATGA overlap — an ordinary operon arrangement this file already treats as
    # noise (it collapses sub-OVERPRINT_MIN_OVERLAP_BP overlaps onto one row).
    if _overprint_partners(genes):
        handles += [Patch(facecolor="none", edgecolor="#33373d", linewidth=1.2, linestyle="--",
                          label="overprinting partner (antisense)")]
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
        seg.add_feature(g["start"] - wlo, g["end"] - wlo + 1, int(g.get("strand", 1)),
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
