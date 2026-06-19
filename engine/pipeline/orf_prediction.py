"""
pipeline/orf_prediction.py — ORF prediction: Prodigal or 6-frame translation.

Provides two prediction backends:
  1. Prodigal  — best for assembled contigs; uses meta mode for short / mixed inputs.
  2. Six-frame — delegates to scripts/04_translate_sixframe.py for raw nucleotides.
"""
import shutil
import subprocess
import sys
from pathlib import Path

from .utils import find_tool, run_cmd


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def predict_orfs_prodigal(
    fna_path: Path,
    out_dir: Path,
    cpu: int = 4,
) -> Path:
    """Run Prodigal in metagenomic mode and return the translated .faa output.

    Uses ``-p meta`` which is safe for both short contigs and longer assemblies.
    Prodigal is single-threaded; the *cpu* parameter is accepted for API
    compatibility but not forwarded to Prodigal itself.

    Parameters
    ----------
    fna_path : Path
        Input nucleotide FASTA.
    out_dir : Path
        Directory for output files (created if missing).
    cpu : int
        Ignored for Prodigal; kept for API symmetry.

    Returns
    -------
    Path
        Path to the output .faa file, or an empty Path on failure.
    """
    fna_path = Path(fna_path)
    out_dir = Path(out_dir)

    if not fna_path.exists():
        print(f"ERROR: Input file not found: {fna_path}", file=sys.stderr)
        return Path()

    prodigal_bin = find_tool("prodigal")
    if prodigal_bin is None:
        print(
            "WARNING: prodigal not found on PATH. "
            "Install Prodigal or use method='sixframe'.",
            file=sys.stderr,
        )
        return Path()

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = fna_path.stem.replace(".fna", "").replace(".fa", "")
    faa_out = out_dir / f"{stem}_prodigal.faa"
    gff_out = out_dir / f"{stem}_prodigal.gff"

    cmd = [
        prodigal_bin,
        "-i", str(fna_path),
        "-a", str(faa_out),
        "-o", str(gff_out),
        "-f", "gff",
        "-p", "meta",
        "-q",
    ]

    result = run_cmd(cmd)
    if result.returncode != 0:
        print(
            f"ERROR: prodigal failed (exit {result.returncode}):\n{result.stderr}",
            file=sys.stderr,
        )
        return Path()

    if not faa_out.exists() or faa_out.stat().st_size == 0:
        print("WARNING: prodigal produced an empty output.", file=sys.stderr)
        return Path()

    return faa_out


def predict_orfs_sixframe(
    fna_path: Path,
    out_dir: Path,
    scripts_dir: Path,
    min_aa: int = 30,
) -> Path:
    """Six-frame translate a nucleotide FASTA into all possible ORFs.

    Tries the external ``04_translate_sixframe.py`` script first; if it is
    not available, falls back to an inline Biopython implementation that
    requires no additional dependencies beyond what the app already uses.

    Parameters
    ----------
    fna_path : Path
        Input nucleotide FASTA.
    out_dir : Path
        Directory for output files (created if missing).
    scripts_dir : Path
        Directory that *may* contain ``04_translate_sixframe.py``
        (used for the script fallback only).
    min_aa : int
        Minimum ORF length in amino acids.

    Returns
    -------
    Path
        Path to the output .faa file, or an empty Path on failure.
    """
    fna_path = Path(fna_path)
    out_dir = Path(out_dir)
    scripts_dir = Path(scripts_dir)

    if not fna_path.exists():
        print(f"ERROR: Input file not found: {fna_path}", file=sys.stderr)
        return Path()

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = fna_path.stem
    faa_out = out_dir / f"{stem}_sixframe.faa"

    # ── Try external script first ────────────────────────────────────
    translate_script = scripts_dir / "04_translate_sixframe.py"
    if translate_script.exists():
        cmd = [
            sys.executable,
            str(translate_script),
            str(fna_path),
            str(faa_out),
            "--min-aa", str(min_aa),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0 and faa_out.exists() and faa_out.stat().st_size > 0:
            return faa_out
        print(
            f"WARNING: 04_translate_sixframe.py failed — using built-in fallback.\n"
            f"  stderr: {result.stderr[:200]}",
            file=sys.stderr,
        )

    # ── Inline Biopython 6-frame translation (always available) ──────
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
    except ImportError:
        print("ERROR: Biopython not installed — cannot do 6-frame translation.", file=sys.stderr)
        return Path()

    import warnings
    from Bio import BiopythonWarning

    n_written = 0
    try:
        with faa_out.open("w") as fh:
            for rec in SeqIO.parse(str(fna_path), "fasta"):
                seq = rec.seq
                for strand, nuc in [(1, seq), (-1, seq.reverse_complement())]:
                    for frame in range(3):
                        frame_seq = nuc[frame:]
                        usable = len(frame_seq) - (len(frame_seq) % 3)
                        if usable <= 0:
                            continue
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", BiopythonWarning)
                            trans = str(frame_seq[:usable].translate())
                        for i, peptide in enumerate(trans.split("*")):
                            if len(peptide) >= min_aa:
                                orf_id = f"{rec.id}_s{strand}_f{frame}_o{i}"
                                fh.write(f">{orf_id}\n{peptide}\n")
                                n_written += 1
    except Exception as exc:
        print(f"ERROR: Inline 6-frame translation failed: {exc}", file=sys.stderr)
        return Path()

    if n_written == 0:
        print("WARNING: 6-frame translation produced no ORFs meeting the length cutoff.", file=sys.stderr)
        return Path()

    return faa_out


def choose_and_predict(
    fna_path: Path,
    out_dir: Path,
    method: str,
    scripts_dir: Path,
    **kwargs,
) -> Path:
    """Dispatch ORF prediction to the requested backend.

    Parameters
    ----------
    fna_path : Path
        Input nucleotide FASTA.
    out_dir : Path
        Output directory.
    method : str
        ``"prodigal"`` or ``"sixframe"``.
    scripts_dir : Path
        Parent scripts directory (used only for ``"sixframe"``).
    **kwargs
        Extra keyword arguments forwarded to the chosen backend
        (e.g. ``cpu``, ``min_aa``).

    Returns
    -------
    Path
        Path to output .faa, or empty Path on failure.
    """
    method = method.lower().strip()

    if method == "prodigal":
        return predict_orfs_prodigal(
            fna_path=fna_path,
            out_dir=out_dir,
            cpu=kwargs.get("cpu", 4),
        )
    elif method in ("sixframe", "six_frame", "6frame"):
        return predict_orfs_sixframe(
            fna_path=fna_path,
            out_dir=out_dir,
            scripts_dir=scripts_dir,
            min_aa=kwargs.get("min_aa", 30),
        )
    else:
        print(
            f"ERROR: Unknown ORF prediction method '{method}'. "
            "Choose 'prodigal' or 'sixframe'.",
            file=sys.stderr,
        )
        return Path()
