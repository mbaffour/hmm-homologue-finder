"""
pipeline/motifs.py — MEME/FIMO motif discovery.

Optional: enabled only when meme is on PATH.

Runs MEME on a protein FASTA to discover de-novo motifs, then uses FIMO
to scan those motifs against a target FASTA.  Both tools must be part of
the MEME Suite (https://meme-suite.org/).
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

from .utils import find_tool, run_cmd
from typing import Optional

import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run_command(cmd: list[str], timeout: int = 7200) -> tuple[int, str, str]:
    """Run a shell command with augmented PATH; return (returncode, stdout, stderr)."""
    try:
        result = run_cmd(cmd, timeout=timeout)
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return 1, "", "Command timed out"
    except FileNotFoundError:
        return 127, "", f"Binary not found: {cmd[0]}"
    except Exception as exc:
        return 1, "", str(exc)


def _meme_available() -> bool:
    return find_tool("meme") is not None


def _fimo_available() -> bool:
    return find_tool("fimo") is not None


# ---------------------------------------------------------------------------
# MEME
# ---------------------------------------------------------------------------

def run_meme(
    faa_path: Path,
    out_dir: Path,
    n_motifs: int = 5,
    min_width: int = 6,
    max_width: int = 50,
    cpu: int = 4,
) -> dict:
    """Run MEME de-novo motif discovery on a protein FASTA.

    Parameters
    ----------
    faa_path : Path
        Input protein FASTA (all-uppercase, no gaps).
    out_dir : Path
        Output directory (passed to ``-oc``).
    n_motifs : int
        Number of motifs to find (``-nmotifs``).
    min_width : int
        Minimum motif width (``-minw``).
    max_width : int
        Maximum motif width (``-maxw``).
    cpu : int
        Parallel processes (``-p``).

    Returns
    -------
    dict
        {meme_dir: Path, n_motifs_found: int, motifs: list[dict],
         success: bool, error: str}
        Each motif dict: {id, consensus, evalue, width, nsites}.
    """
    faa_path = Path(faa_path)
    out_dir  = Path(out_dir)

    result: dict = {
        "meme_dir":       out_dir,
        "meme_txt":       None,
        "n_motifs_found": 0,
        "motifs":         [],
        "success":        False,
        "error":          "",
    }

    if not _meme_available():
        result["error"] = "meme not found on PATH"
        print(f"WARNING: {result['error']}", file=sys.stderr)
        return result

    if not faa_path.exists():
        result["error"] = f"Input FASTA not found: {faa_path}"
        return result

    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "meme",
        str(faa_path),
        "-oc",     str(out_dir),
        "-protein",
        "-nmotifs", str(n_motifs),
        "-minw",    str(min_width),
        "-maxw",    str(max_width),
        "-p",       str(cpu),
        "-mod",     "zoops",  # zero or one occurrence per sequence
    ]

    print(f"INFO: Running: {' '.join(cmd)}", file=sys.stderr)
    rc, stdout, stderr = _run_command(cmd)

    if rc != 0:
        result["error"] = f"meme failed (rc={rc}): {stderr[-1000:]}"
        print(f"ERROR: {result['error']}", file=sys.stderr)
        return result

    meme_txt = out_dir / "meme.txt"
    if not meme_txt.exists():
        result["error"] = "meme finished but meme.txt not found"
        return result

    motifs = parse_meme_txt(meme_txt)
    result.update(
        {
            "meme_txt":       meme_txt,
            "n_motifs_found": len(motifs),
            "motifs":         motifs,
            "success":        True,
            "error":          "",
        }
    )
    return result


# ---------------------------------------------------------------------------
# FIMO
# ---------------------------------------------------------------------------

def run_fimo(
    meme_txt: Path,
    faa_path: Path,
    out_dir: Path,
    pthresh: float = 1e-4,
) -> pd.DataFrame:
    """Run FIMO to scan MEME motifs against a protein FASTA.

    Parameters
    ----------
    meme_txt : Path
        ``meme.txt`` output file from a MEME run.
    faa_path : Path
        Target protein FASTA.
    out_dir : Path
        Output directory (passed to ``--oc``).
    pthresh : float
        P-value threshold (``--thresh``).

    Returns
    -------
    pd.DataFrame
        Columns: motif_id, sequence_name, start, stop, strand,
                 score, pvalue, matched_sequence.
        Empty DataFrame on failure.
    """
    meme_txt = Path(meme_txt)
    faa_path = Path(faa_path)
    out_dir  = Path(out_dir)

    _EMPTY_COLS = [
        "motif_id", "sequence_name", "start", "stop",
        "strand", "score", "pvalue", "matched_sequence",
    ]

    if not _fimo_available():
        print("WARNING: fimo not found on PATH; skipping motif scanning.", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_COLS)

    if not meme_txt.exists():
        print(f"ERROR: meme.txt not found: {meme_txt}", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_COLS)

    if not faa_path.exists():
        print(f"ERROR: Target FASTA not found: {faa_path}", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_COLS)

    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "fimo",
        "--thresh",  str(pthresh),
        "--oc",      str(out_dir),
        "--verbosity", "1",
        str(meme_txt),
        str(faa_path),
    ]

    print(f"INFO: Running: {' '.join(cmd)}", file=sys.stderr)
    rc, stdout, stderr = _run_command(cmd)

    if rc != 0:
        print(f"ERROR: fimo failed (rc={rc}): {stderr[-1000:]}", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_COLS)

    # FIMO writes fimo.tsv (or fimo.txt in older versions)
    fimo_tsv = out_dir / "fimo.tsv"
    fimo_txt = out_dir / "fimo.txt"

    if fimo_tsv.exists():
        return _parse_fimo_tsv(fimo_tsv)
    if fimo_txt.exists():
        return _parse_fimo_txt(fimo_txt)

    print("WARNING: fimo finished but neither fimo.tsv nor fimo.txt found.", file=sys.stderr)
    return pd.DataFrame(columns=_EMPTY_COLS)


def _parse_fimo_tsv(tsv_path: Path) -> pd.DataFrame:
    """Parse FIMO TSV output (MEME Suite ≥ 5.x format)."""
    _COLS = [
        "motif_id", "sequence_name", "start", "stop",
        "strand", "score", "pvalue", "matched_sequence",
    ]
    try:
        # Skip comment lines starting with '#'
        rows = []
        for line in tsv_path.read_text(errors="replace").splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            rows.append(
                {
                    "motif_id":        parts[0],
                    "sequence_name":   parts[2],
                    "start":           _safe_int(parts[3]),
                    "stop":            _safe_int(parts[4]),
                    "strand":          parts[5],
                    "score":           _safe_float(parts[6]),
                    "pvalue":          _safe_float(parts[7]),
                    "matched_sequence": parts[9] if len(parts) > 9 else "",
                }
            )
        if not rows:
            return pd.DataFrame(columns=_COLS)
        df = pd.DataFrame(rows)
        return df[_COLS]
    except Exception as exc:
        print(f"ERROR: Cannot parse FIMO TSV: {exc}", file=sys.stderr)
        return pd.DataFrame(columns=_COLS)


def _parse_fimo_txt(txt_path: Path) -> pd.DataFrame:
    """Parse legacy FIMO text output (MEME Suite < 5.x)."""
    _COLS = [
        "motif_id", "sequence_name", "start", "stop",
        "strand", "score", "pvalue", "matched_sequence",
    ]
    try:
        rows = []
        for line in txt_path.read_text(errors="replace").splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            rows.append(
                {
                    "motif_id":        parts[0],
                    "sequence_name":   parts[1],
                    "start":           _safe_int(parts[2]),
                    "stop":            _safe_int(parts[3]),
                    "strand":          parts[4],
                    "score":           _safe_float(parts[5]),
                    "pvalue":          _safe_float(parts[6]),
                    "matched_sequence": parts[7] if len(parts) > 7 else "",
                }
            )
        if not rows:
            return pd.DataFrame(columns=_COLS)
        df = pd.DataFrame(rows)
        return df[_COLS]
    except Exception as exc:
        print(f"ERROR: Cannot parse FIMO text file: {exc}", file=sys.stderr)
        return pd.DataFrame(columns=_COLS)


# ---------------------------------------------------------------------------
# Parse MEME text output
# ---------------------------------------------------------------------------

def parse_meme_txt(meme_txt_path: Path) -> list[dict]:
    """Parse MEME text output for motif summaries.

    Parameters
    ----------
    meme_txt_path : Path
        Path to ``meme.txt``.

    Returns
    -------
    list[dict]
        Each dict: {id, consensus, evalue, width, nsites}.
        Empty list on failure or if no motifs found.
    """
    meme_txt_path = Path(meme_txt_path)
    if not meme_txt_path.exists():
        return []

    try:
        text = meme_txt_path.read_text(errors="replace")
    except Exception as exc:
        print(f"ERROR: Cannot read {meme_txt_path}: {exc}", file=sys.stderr)
        return []

    motifs: list[dict] = []

    # Patterns for the MEME text format
    # Block header: "MOTIF  1  MEME-1"  or  "MOTIF 1 MEME-1"
    motif_block_re = re.compile(
        r"MOTIF\s+(\S+).*?letter-probability matrix.*?w=\s*(\d+).*?nsites=\s*(\d+).*?E=\s*(\S+)",
        re.DOTALL,
    )

    # Simpler per-motif summary section
    summary_re = re.compile(
        r"^MOTIF\s+(\S+)",
        re.MULTILINE,
    )

    # --- Try the letter-probability matrix blocks first ---
    for m in motif_block_re.finditer(text):
        motif_id  = m.group(1)
        width     = _safe_int(m.group(2))
        nsites    = _safe_int(m.group(3))
        evalue    = _safe_float(m.group(4))

        # Extract consensus from IUPAC consensus line immediately after block
        consensus = _extract_consensus_after(text, m.end())

        motifs.append(
            {
                "id":        f"MEME-{motif_id}",
                "consensus": consensus,
                "evalue":    evalue,
                "width":     width,
                "nsites":    nsites,
            }
        )

    if motifs:
        return motifs

    # --- Fallback: minimal summary-line parsing ---
    # "MOTIF 1  width = 15  sites = 42  llr = 382  E-value = 1.2e-30"
    simple_re = re.compile(
        r"MOTIF\s+(\S+)\s+.*?width\s*=\s*(\d+)\s+sites\s*=\s*(\d+)\s+.*?E-value\s*=\s*(\S+)",
        re.DOTALL | re.IGNORECASE,
    )
    for m in simple_re.finditer(text):
        motif_id  = m.group(1)
        width     = _safe_int(m.group(2))
        nsites    = _safe_int(m.group(3))
        evalue    = _safe_float(m.group(4))
        motifs.append(
            {
                "id":        f"MEME-{motif_id}",
                "consensus": "",
                "evalue":    evalue,
                "width":     width,
                "nsites":    nsites,
            }
        )

    return motifs


def _extract_consensus_after(text: str, pos: int, max_chars: int = 500) -> str:
    """Look for a consensus/IUPAC line in the text after position ``pos``."""
    window = text[pos : pos + max_chars]

    # MEME sometimes writes consensus as a line of amino-acid characters
    # immediately following the probability matrix header
    cons_re = re.compile(r"\n([A-Z]{4,})\n")
    m = cons_re.search(window)
    if m:
        return m.group(1)

    return ""


# ---------------------------------------------------------------------------
# Small numeric helpers
# ---------------------------------------------------------------------------

def _safe_int(value: str, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except (ValueError, TypeError):
        return default


def _safe_float(value: str, default: float = 0.0) -> float:
    try:
        return float(str(value).strip())
    except (ValueError, TypeError):
        return default
