"""
pipeline/alignment.py — MSA (MAFFT/Clustal Omega) + trimAl + quality metrics.

Wraps external alignment tools and provides alignment quality statistics
and ASCII previews suitable for display in the Shiny app.
"""
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

from .utils import find_tool, run_cmd

try:
    from Bio import AlignIO
    from Bio.Align import MultipleSeqAlignment
    _BIO_AVAILABLE = True
except ImportError:
    _BIO_AVAILABLE = False


# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------

def accuracy_flags(n_seqs: int, max_linsi: int = 500) -> "tuple[list, str]":
    """MAFFT flags for the most accurate strategy that stays computationally
    tractable for the given sequence count.

    L-INS-i (``--localpair --maxiterate 1000``) is the gold-standard MAFFT
    strategy for aligning homologous domains embedded in variable-length context
    — the situation here. It is O(N^2) in pairwise comparisons, so it is used up
    to ``max_linsi`` sequences and the pipeline falls back to MAFFT's ``--auto``
    (which selects FFT-NS-2 / PartTree) for larger sets to stay practical.

    Returns ``(extra_flags, human_label)``.
    """
    if 2 <= n_seqs <= max_linsi:
        return (["--maxiterate", "1000", "--localpair"], f"L-INS-i (accurate, n={n_seqs})")
    return ([], f"auto/FFT-NS-2 (n={n_seqs} > {max_linsi})")


def run_mafft(
    faa_path: Path,
    out_path: Path,
    cpu: int = 4,
    extra_flags: Optional[list] = None,
) -> Path:
    """Run MAFFT and return the path to the aligned output.

    Parameters
    ----------
    faa_path : Path
        Input un-aligned protein FASTA.
    out_path : Path
        Destination path for the aligned FASTA.
    cpu : int
        Number of threads (``--thread``).
    extra_flags : list, optional
        Explicit MAFFT strategy flags, e.g. ``["--localpair", "--maxiterate",
        "1000"]`` (L-INS-i). When given, these REPLACE ``--auto`` so the requested
        accurate strategy is honoured (passing both ``--auto`` and ``--localpair``
        lets ``--auto`` silently override the accurate choice). When omitted,
        ``--auto`` is used.

    Returns
    -------
    Path
        Path to the aligned file, or empty Path on failure.
    """
    faa_path = Path(faa_path)
    out_path = Path(out_path)

    if not faa_path.exists():
        print(f"ERROR: Input not found: {faa_path}", file=sys.stderr)
        return Path()

    mafft_bin = find_tool("mafft")
    if mafft_bin is None:
        print("ERROR: mafft not found on PATH.", file=sys.stderr)
        return Path()

    out_path.parent.mkdir(parents=True, exist_ok=True)

    strategy = list(extra_flags) if extra_flags else ["--auto"]
    cmd = [mafft_bin, "--thread", str(cpu), *strategy, str(faa_path)]

    result = run_cmd(cmd)
    if result.returncode != 0:
        print(f"ERROR: mafft failed:\n{result.stderr}", file=sys.stderr)
        return Path()

    out_path.write_text(result.stdout)

    if out_path.stat().st_size == 0:
        print("WARNING: mafft produced empty output.", file=sys.stderr)
        return Path()

    return out_path


def run_clustalo(
    faa_path: Path,
    out_path: Path,
    cpu: int = 4,
) -> Path:
    """Run Clustal Omega and return the path to the aligned output.

    Parameters
    ----------
    faa_path : Path
        Input un-aligned protein FASTA.
    out_path : Path
        Destination path for the aligned FASTA.
    cpu : int
        Number of threads (``--threads``).

    Returns
    -------
    Path
        Path to the aligned file, or empty Path on failure.
    """
    faa_path = Path(faa_path)
    out_path = Path(out_path)

    if not faa_path.exists():
        print(f"ERROR: Input not found: {faa_path}", file=sys.stderr)
        return Path()

    clustalo_bin = find_tool("clustalo")
    if clustalo_bin is None:
        print("ERROR: clustalo not found on PATH.", file=sys.stderr)
        return Path()

    out_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        clustalo_bin,
        "-i", str(faa_path),
        "-o", str(out_path),
        "--threads", str(cpu),
        "--force",
    ]

    result = run_cmd(cmd)
    if result.returncode != 0:
        print(f"ERROR: clustalo failed:\n{result.stderr}", file=sys.stderr)
        return Path()

    if not out_path.exists() or out_path.stat().st_size == 0:
        print("WARNING: clustalo produced empty output.", file=sys.stderr)
        return Path()

    return out_path


def run_trimal(
    aln_path: Path,
    out_path: Path,
    method: str = "automated1",
) -> Path:
    """Run trimAl to trim a multiple sequence alignment.

    Parameters
    ----------
    aln_path : Path
        Input aligned FASTA.
    out_path : Path
        Output trimmed FASTA.
    method : str
        trimAl method flag without the leading dash.
        Common values: ``"automated1"``, ``"gappyout"``, ``"strict"``.

    Returns
    -------
    Path
        Path to trimmed file, or empty Path on failure.
    """
    aln_path = Path(aln_path)
    out_path = Path(out_path)

    if not aln_path.exists():
        print(f"ERROR: Alignment not found: {aln_path}", file=sys.stderr)
        return Path()

    trimal_bin = find_tool("trimal")
    if trimal_bin is None:
        print("ERROR: trimal not found on PATH.", file=sys.stderr)
        return Path()

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # trimal expects the flag as -automated1, -gappyout, etc.
    method_flag = method if method.startswith("-") else f"-{method}"
    cmd = [
        trimal_bin,
        "-in", str(aln_path),
        "-out", str(out_path),
        method_flag,
    ]

    result = run_cmd(cmd)
    if result.returncode != 0:
        print(f"ERROR: trimal failed:\n{result.stderr}", file=sys.stderr)
        return Path()

    if not out_path.exists() or out_path.stat().st_size == 0:
        print("WARNING: trimal produced empty output.", file=sys.stderr)
        return Path()

    return out_path


# ---------------------------------------------------------------------------
# Quality metrics
# ---------------------------------------------------------------------------

def alignment_quality(aln_path: Path) -> dict:
    """Compute quality statistics for a multiple sequence alignment.

    Parameters
    ----------
    aln_path : Path
        Aligned FASTA file.

    Returns
    -------
    dict with keys:
        aln_length, n_sequences, gap_pct, conserved_columns,
        avg_pairwise_id, flagged_sequences
    """
    aln_path = Path(aln_path)

    empty_result = {
        "aln_length": 0,
        "n_sequences": 0,
        "gap_pct": 0.0,
        "conserved_columns": 0,
        "avg_pairwise_id": 0.0,
        "flagged_sequences": [],
    }

    if not aln_path.exists():
        print(f"ERROR: Alignment not found: {aln_path}", file=sys.stderr)
        return empty_result

    if not _BIO_AVAILABLE:
        print("ERROR: Biopython not installed.", file=sys.stderr)
        return empty_result

    try:
        aln = AlignIO.read(str(aln_path), "fasta")
    except Exception as exc:
        print(f"ERROR: Cannot read alignment {aln_path}: {exc}", file=sys.stderr)
        return empty_result

    n_seqs = len(aln)
    aln_len = aln.get_alignment_length()
    if n_seqs == 0 or aln_len == 0:
        return empty_result

    # Gap percentage across entire alignment
    total_chars = n_seqs * aln_len
    total_gaps = sum(str(rec.seq).count("-") for rec in aln)
    gap_pct = round(100.0 * total_gaps / total_chars, 2)

    # Conserved columns (single dominant residue ≥ 80% non-gap)
    conserved = 0
    for col_idx in range(aln_len):
        col = [aln[row_idx, col_idx] for row_idx in range(n_seqs)]
        non_gap = [c for c in col if c != "-"]
        if not non_gap:
            continue
        most_common_frac = max(non_gap.count(c) for c in set(non_gap)) / len(non_gap)
        if most_common_frac >= 0.80:
            conserved += 1

    # Average pairwise identity (sample up to 200 pairs for speed)
    import itertools, random
    pairs = list(itertools.combinations(range(n_seqs), 2))
    if len(pairs) > 200:
        pairs = random.sample(pairs, 200)

    identities = []
    for i, j in pairs:
        seq_i = str(aln[i].seq)
        seq_j = str(aln[j].seq)
        matches = sum(
            a == b for a, b in zip(seq_i, seq_j) if a != "-" and b != "-"
        )
        compared = sum(1 for a, b in zip(seq_i, seq_j) if a != "-" and b != "-")
        if compared > 0:
            identities.append(matches / compared)

    avg_pairwise_id = round(100.0 * (sum(identities) / len(identities)), 2) if identities else 0.0

    # Flag sequences with >80% gaps
    flagged = []
    for rec in aln:
        gap_frac = str(rec.seq).count("-") / aln_len
        if gap_frac > 0.80:
            flagged.append(rec.id)

    return {
        "aln_length": aln_len,
        "n_sequences": n_seqs,
        "gap_pct": gap_pct,
        "conserved_columns": conserved,
        "avg_pairwise_id": avg_pairwise_id,
        "flagged_sequences": flagged,
    }


def alignment_figure(
    aln_path: Path,
    out_dir: Path,
    max_seqs: int = 60,
    max_cols: int = 300,
    fmt: str = "png",
    dpi: int = 300,
) -> bytes:
    """
    Export a publication-ready coloured multiple-sequence alignment image.

    Amino acids are coloured by physicochemical class (ClustalX scheme).
    Residues shown at full opacity; gaps are white.

    Parameters
    ----------
    aln_path : Path
        Aligned FASTA (output of MAFFT / trimAl).
    out_dir : Path
        Directory where ``alignment_figure.<fmt>`` is saved.
    max_seqs : int
        Maximum sequences to display (top rows).
    max_cols : int
        Maximum alignment columns to display.
    fmt : str
        ``"png"`` (300 dpi) · ``"svg"`` · ``"pdf"``.
    dpi : int
        PNG resolution (ignored for SVG/PDF).

    Returns
    -------
    bytes  in the requested format, or ``b""`` on failure.
    """
    aln_path = Path(aln_path)
    out_dir  = Path(out_dir)

    if not aln_path.exists():
        print(f"ERROR: Alignment not found: {aln_path}", file=sys.stderr)
        return b""

    if not _BIO_AVAILABLE:
        print("ERROR: Biopython not installed.", file=sys.stderr)
        return b""

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import numpy as np
    except ImportError:
        print("ERROR: matplotlib/numpy not installed.", file=sys.stderr)
        return b""

    # ── ClustalX amino acid colour scheme ─────────────────────────────────
    # Matches the scheme used in Jalview, MUSCLE output, and most journals.
    _AA_COLORS: dict[str, str] = {
        # Hydrophobic
        "A": "#80a0f0", "V": "#80a0f0", "I": "#80a0f0",
        "L": "#80a0f0", "M": "#80a0f0", "F": "#80a0f0",
        "W": "#80a0f0", "P": "#ffff00",
        # Positively charged
        "K": "#f01505", "R": "#f01505", "H": "#15c015",
        # Negatively charged
        "D": "#c048c0", "E": "#c048c0",
        # Polar uncharged
        "S": "#15c015", "T": "#15c015", "N": "#15c015", "Q": "#15c015",
        # Cysteine / aromatic
        "C": "#f08080", "Y": "#15a8a8", "G": "#f09048",
        # Stop / unknown
        "*": "#ffffff", "X": "#cccccc", "B": "#cccccc", "Z": "#cccccc",
        "-": "#ffffff",  # gap = white
    }
    _DEFAULT_COL = "#eeeeee"

    try:
        from Bio import AlignIO
        aln  = AlignIO.read(str(aln_path), "fasta")
    except Exception as exc:
        print(f"ERROR: Cannot read alignment: {exc}", file=sys.stderr)
        return b""

    records  = list(aln)[:max_seqs]
    aln_len  = min(aln.get_alignment_length(), max_cols)
    n_seqs   = len(records)

    if n_seqs == 0 or aln_len == 0:
        return b""

    # ── Build colour matrix ─────────────────────────────────────────────────
    import matplotlib.colors as mcolors

    color_mat = np.zeros((n_seqs, aln_len, 4), dtype=float)   # RGBA
    for i, rec in enumerate(records):
        for j, aa in enumerate(str(rec.seq)[:aln_len]):
            hex_c = _AA_COLORS.get(aa.upper(), _DEFAULT_COL)
            r, g, b = mcolors.to_rgb(hex_c)
            color_mat[i, j] = [r, g, b, 1.0 if aa != "-" else 0.0]

    # ── Per-column conservation (CLC-style abundance track): fraction of
    #    sequences sharing the most common residue in each column. ────────────
    from collections import Counter
    conservation = np.zeros(aln_len)
    consensus = []                                   # most-common residue per column
    for j in range(aln_len):
        col = [str(rec.seq)[j] if j < len(rec.seq) else "-" for rec in records]
        non_gap = [c for c in col if c != "-"]
        if non_gap:
            res, cnt = Counter(non_gap).most_common(1)[0]
            conservation[j] = cnt / n_seqs
            consensus.append(res.upper())
        else:
            consensus.append("-")

    # ── Figure layout: alignment (top) + conservation bar chart (bottom) ─────
    cell_w  = 0.11   # inch per column
    cell_h  = 0.16   # inch per row
    fig_w   = max(8.0,  min(aln_len  * cell_w + 3.2, 26.0))
    fig_h   = max(3.8,  min(n_seqs   * cell_h + 3.4, 22.0))

    plt.rcParams.update({
        "font.family": "monospace",
        "font.size":   5,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })

    fig = plt.figure(figsize=(fig_w, fig_h + 0.5))
    gs  = fig.add_gridspec(3, 1, height_ratios=[max(n_seqs, 6), 1.2, 4], hspace=0.06)
    ax  = fig.add_subplot(gs[0])                 # alignment
    axk = fig.add_subplot(gs[1], sharex=ax)      # consensus row
    axc = fig.add_subplot(gs[2], sharex=ax)      # conservation / abundance bars

    ax.imshow(color_mat, aspect="auto", interpolation="none",
              origin="upper", extent=[0, aln_len, n_seqs, 0])

    # Amino-acid LETTERS overlaid on the cells (CLC-style). Drawn when the grid is
    # legible; vector outputs (SVG/PDF) stay zoomable in Inkscape regardless.
    draw_letters = aln_len <= 200 and n_seqs <= 70
    if draw_letters:
        lfs = max(3.0, min(7.5, 110.0 / aln_len))
        for i, rec in enumerate(records):
            s = str(rec.seq)
            for j in range(aln_len):
                aa = s[j] if j < len(s) else "-"
                if aa != "-":
                    ax.text(j + 0.5, i + 0.5, aa, ha="center", va="center",
                            fontsize=lfs, color="#1a1a1a")

    # Sequence ID labels on the left
    ax.set_yticks([i + 0.5 for i in range(n_seqs)])
    ax.set_yticklabels([r.id[:38] for r in records], fontsize=5.5)
    ax.set_xlim(0, aln_len)
    plt.setp(ax.get_xticklabels(), visible=False)   # x-axis lives on the bottom panel
    ax.set_title(
        f"Multiple Sequence Alignment  ({n_seqs} sequences × {aln_len} columns"
        + ("" if draw_letters else "; residues shown as colour blocks — zoom the SVG for letters")
        + ")",
        fontsize=8, fontweight="bold", pad=5,
    )
    ax.tick_params(axis="both", which="both", length=2, width=0.5)

    # ── Consensus row (most common residue per column), CLC-style ────────────
    cons_mat = np.zeros((1, aln_len, 4), dtype=float)
    for j, aa in enumerate(consensus):
        r, g, b = mcolors.to_rgb(_AA_COLORS.get(aa, _DEFAULT_COL))
        cons_mat[0, j] = [r, g, b, 1.0 if aa != "-" else 0.0]
    axk.imshow(cons_mat, aspect="auto", interpolation="none", origin="upper",
               extent=[0, aln_len, 1, 0])
    if draw_letters:
        for j, aa in enumerate(consensus):
            if aa != "-":
                axk.text(j + 0.5, 0.5, aa, ha="center", va="center",
                         fontsize=lfs, fontweight="bold", color="#1a1a1a")
    axk.set_yticks([0.5]); axk.set_yticklabels(["consensus"], fontsize=6)
    axk.set_xlim(0, aln_len)
    plt.setp(axk.get_xticklabels(), visible=False)
    axk.tick_params(axis="both", which="both", length=0)

    # Conservation / abundance bar chart, one bar per column, shaded by level.
    cmap = plt.cm.YlGn
    axc.bar(np.arange(aln_len) + 0.5, conservation, width=1.0,
            color=[cmap(0.25 + 0.75 * float(c)) for c in conservation], edgecolor="none")
    axc.set_ylim(0, 1.0)
    axc.set_ylabel("Conservation\n(per column)", fontsize=6)
    axc.set_yticks([0, 0.5, 1.0])
    axc.tick_params(labelsize=5)
    tick_positions = list(range(0, aln_len, 50)) + [aln_len]
    axc.set_xticks(tick_positions)
    axc.set_xticklabels([str(p + 1) for p in tick_positions], fontsize=6)
    axc.set_xlabel("Alignment column", fontsize=7)
    for sp in ("top", "right"):
        axc.spines[sp].set_visible(False)

    # Colour legend
    legend_aa_groups = [
        ("Hydrophobic (A/V/I/L/M/F/W)", "#80a0f0"),
        ("Proline (P)",                  "#ffff00"),
        ("Pos. charged (K/R)",           "#f01505"),
        ("His (H) / Polar (S/T/N/Q)",    "#15c015"),
        ("Neg. charged (D/E)",           "#c048c0"),
        ("Cys (C)",                      "#f08080"),
        ("Tyr (Y)",                      "#15a8a8"),
        ("Gly (G)",                      "#f09048"),
        ("Gap",                          "#ffffff"),
    ]
    handles = [
        mpatches.Patch(facecolor=col, edgecolor="#888", linewidth=0.5, label=lbl)
        for lbl, col in legend_aa_groups
    ]
    ax.legend(
        handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.01),
        ncol=3, fontsize=5, framealpha=0.9, edgecolor="#cccccc",
        handlelength=1.0, borderpad=0.4, columnspacing=0.6,
    )

    import io as _io
    buf = _io.BytesIO()
    fmt_lower = fmt.lower().lstrip(".")
    if fmt_lower not in ("png", "svg", "pdf"):
        fmt_lower = "png"
    save_kw: dict = {"format": fmt_lower, "bbox_inches": "tight"}
    if fmt_lower == "png":
        save_kw["dpi"] = dpi
    fig.savefig(buf, **save_kw)
    plt.close(fig)
    buf.seek(0)
    data = buf.read()

    # Save to disk
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"alignment_figure.{fmt_lower}"
    out_path.write_bytes(data)

    return data


def alignment_preview(
    aln_path: Path,
    max_seqs: int = 20,
    max_cols: int = 80,
) -> str:
    """Return an ASCII text grid preview of an alignment.

    Parameters
    ----------
    aln_path : Path
        Aligned FASTA file.
    max_seqs : int
        Maximum number of sequences to display.
    max_cols : int
        Maximum number of alignment columns to display.

    Returns
    -------
    str
        Formatted text grid, empty string on failure.
    """
    aln_path = Path(aln_path)
    if not aln_path.exists():
        return ""

    if not _BIO_AVAILABLE:
        return "ERROR: Biopython not installed."

    try:
        aln = AlignIO.read(str(aln_path), "fasta")
    except Exception as exc:
        return f"ERROR: {exc}"

    if not aln:
        return ""

    records = list(aln)[:max_seqs]
    id_width = min(max(len(r.id) for r in records), 20)

    lines = []
    header = " " * (id_width + 2) + "".join(
        str((i // 10) % 10) for i in range(min(max_cols, aln.get_alignment_length()))
    )
    lines.append(header)

    for rec in records:
        label = rec.id[:id_width].ljust(id_width)
        seq = str(rec.seq)[:max_cols]
        lines.append(f"{label}  {seq}")

    return "\n".join(lines)
