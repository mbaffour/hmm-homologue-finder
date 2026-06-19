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
        Additional MAFFT flags, e.g. ``["--localpair", "--maxiterate", "1000"]``.

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

    cmd = [mafft_bin, "--auto", "--thread", str(cpu)]
    if extra_flags:
        cmd.extend(extra_flags)
    cmd.append(str(faa_path))

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

    # ── Figure layout ───────────────────────────────────────────────────────
    cell_w  = 0.10   # inch per column — gives ~30mm per 300 cols
    cell_h  = 0.14   # inch per row
    fig_w   = max(8.0,  min(aln_len  * cell_w + 3.0, 24.0))
    fig_h   = max(3.0,  min(n_seqs   * cell_h + 1.8, 20.0))

    plt.rcParams.update({
        "font.family": "monospace",
        "font.size":   5,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(color_mat, aspect="auto", interpolation="none",
              origin="upper", extent=[0, aln_len, n_seqs, 0])

    # Sequence ID labels on left
    ax.set_yticks([i + 0.5 for i in range(n_seqs)])
    ax.set_yticklabels([r.id[:28] for r in records], fontsize=5.5)
    ax.set_ylabel("", fontsize=7)

    # Column position tick marks every 50 columns
    tick_positions = list(range(0, aln_len, 50)) + [aln_len]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([str(p + 1) for p in tick_positions], fontsize=6)
    ax.set_xlabel("Alignment column", fontsize=7)

    ax.set_title(
        f"Multiple Sequence Alignment  ({n_seqs} sequences × {aln_len} columns)",
        fontsize=8, fontweight="bold", pad=5,
    )
    ax.tick_params(axis="both", which="both", length=2, width=0.5)

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

    plt.tight_layout(pad=0.8)

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
