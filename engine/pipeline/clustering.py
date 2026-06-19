"""
pipeline/clustering.py — Sequence clustering (CD-HIT or MMseqs2).

Optional: enabled only when cd-hit or mmseqs is on PATH.  Both tools produce
a membership DataFrame with columns: protein_id, cluster_id, is_representative.

Dispatch order:  cd-hit  →  mmseqs2  →  empty stub result
"""
from __future__ import annotations

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

_EMPTY_MEMBERSHIP_COLS = ["protein_id", "cluster_id", "is_representative"]


def _empty_result() -> dict:
    return {
        "cluster_file":   None,
        "rep_faa":        None,
        "n_clusters":     0,
        "membership_df":  pd.DataFrame(columns=_EMPTY_MEMBERSHIP_COLS),
        "error":          "",
    }


def _run_command(cmd: list[str], timeout: int = 3600) -> tuple[int, str, str]:
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


# ---------------------------------------------------------------------------
# CD-HIT
# ---------------------------------------------------------------------------

def cluster_cdhit(
    faa_path: Path,
    out_dir: Path,
    identity: float = 0.40,
    coverage: float = 0.80,
    threads: int = 4,
) -> dict:
    """Cluster sequences with CD-HIT.

    Parameters
    ----------
    faa_path : Path
        Input protein FASTA.
    out_dir : Path
        Output directory.
    identity : float
        Sequence identity threshold (0–1).  CD-HIT requires ≥0.40 for ``-n 2``.
    coverage : float
        Alignment coverage threshold (shorter sequence, ``-aL``).
    threads : int
        CPU threads (``-T``).

    Returns
    -------
    dict
        {cluster_file: Path, rep_faa: Path, n_clusters: int,
         membership_df: pd.DataFrame, error: str}
    """
    faa_path = Path(faa_path)
    out_dir  = Path(out_dir)
    result   = _empty_result()

    if not find_tool("cd-hit"):
        result["error"] = "cd-hit not found on PATH"
        return result

    if not faa_path.exists():
        result["error"] = f"Input FASTA not found: {faa_path}"
        return result

    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / "cdhit_clusters"
    clstr_file = Path(str(out_prefix) + ".clstr")

    # CD-HIT word size depends on identity threshold
    if   identity >= 0.70: n_word = 5
    elif identity >= 0.60: n_word = 4
    elif identity >= 0.50: n_word = 3
    else:                  n_word = 2

    cmd = [
        "cd-hit",
        "-i", str(faa_path),
        "-o", str(out_prefix),
        "-c", str(identity),
        "-aL", str(coverage),
        "-T", str(threads),
        "-M", "4000",
        "-n", str(n_word),
        "-d", "0",       # full sequence name in .clstr
    ]

    print(f"INFO: Running: {' '.join(cmd)}", file=sys.stderr)
    rc, stdout, stderr = _run_command(cmd)

    if rc != 0:
        result["error"] = f"cd-hit failed (rc={rc}): {stderr[-1000:]}"
        print(f"ERROR: {result['error']}", file=sys.stderr)
        return result

    if not clstr_file.exists():
        result["error"] = "cd-hit finished but .clstr file not found"
        return result

    membership_df = _parse_cdhit_clstr(clstr_file)
    n_clusters    = membership_df["cluster_id"].nunique() if not membership_df.empty else 0

    result.update(
        {
            "cluster_file":  clstr_file,
            "rep_faa":       out_prefix,   # CD-HIT writes rep sequences to out_prefix (no extension)
            "n_clusters":    n_clusters,
            "membership_df": membership_df,
            "error":         "",
        }
    )
    return result


def _parse_cdhit_clstr(clstr_file: Path) -> pd.DataFrame:
    """Parse a CD-HIT ``.clstr`` file into a membership DataFrame.

    Returns
    -------
    pd.DataFrame
        Columns: protein_id, cluster_id, is_representative.
    """
    rows: list[dict] = []
    current_cluster: Optional[int] = None

    try:
        text = clstr_file.read_text(errors="replace")
    except Exception as exc:
        print(f"ERROR: Cannot read {clstr_file}: {exc}", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_MEMBERSHIP_COLS)

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">Cluster"):
            try:
                current_cluster = int(line.split()[1])
            except (IndexError, ValueError):
                current_cluster = len(rows)
            continue

        if current_cluster is None:
            continue

        # Example line: 0  123aa, >WP_001234567.1... *
        is_rep = line.endswith("*")
        # Extract protein ID from ">ID..."
        if ">" not in line:
            continue
        try:
            pid_part = line.split(">")[1].split("...")[0].strip()
        except IndexError:
            continue

        rows.append(
            {
                "protein_id":       pid_part,
                "cluster_id":       current_cluster,
                "is_representative": is_rep,
            }
        )

    if not rows:
        return pd.DataFrame(columns=_EMPTY_MEMBERSHIP_COLS)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# MMseqs2
# ---------------------------------------------------------------------------

def cluster_mmseqs(
    faa_path: Path,
    out_dir: Path,
    identity: float = 0.40,
    coverage: float = 0.80,
    threads: int = 4,
) -> dict:
    """Cluster sequences with MMseqs2 ``easy-cluster``.

    Parameters
    ----------
    faa_path : Path
        Input protein FASTA.
    out_dir : Path
        Output directory; temp files go in ``out_dir/tmp``.
    identity : float
        Minimum sequence identity (``--min-seq-id``).
    coverage : float
        Minimum coverage (``-c``).
    threads : int
        CPU threads.

    Returns
    -------
    dict
        {cluster_file: Path, rep_faa: Path, n_clusters: int,
         membership_df: pd.DataFrame, error: str}
    """
    faa_path = Path(faa_path)
    out_dir  = Path(out_dir)
    result   = _empty_result()

    binary = find_tool("mmseqs")
    if not binary:
        result["error"] = "mmseqs not found on PATH"
        return result

    if not faa_path.exists():
        result["error"] = f"Input FASTA not found: {faa_path}"
        return result

    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir    = out_dir / "mmseqs_tmp"
    out_prefix = out_dir / "mmseqs_clusters"

    tmp_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "mmseqs", "easy-cluster",
        str(faa_path),
        str(out_prefix),
        str(tmp_dir),
        "--min-seq-id", str(identity),
        "-c", str(coverage),
        "--threads", str(threads),
        "--cluster-mode", "0",
        "-v", "1",
    ]

    print(f"INFO: Running: {' '.join(cmd)}", file=sys.stderr)
    rc, stdout, stderr = _run_command(cmd)

    if rc != 0:
        result["error"] = f"mmseqs failed (rc={rc}): {stderr[-1000:]}"
        print(f"ERROR: {result['error']}", file=sys.stderr)
        return result

    # MMseqs2 easy-cluster writes: PREFIX_cluster.tsv, PREFIX_rep_seq.fasta
    cluster_tsv = Path(str(out_prefix) + "_cluster.tsv")
    rep_faa     = Path(str(out_prefix) + "_rep_seq.fasta")

    if not cluster_tsv.exists():
        result["error"] = "mmseqs finished but cluster TSV not found"
        return result

    membership_df = _parse_mmseqs_cluster_tsv(cluster_tsv)
    n_clusters    = membership_df["cluster_id"].nunique() if not membership_df.empty else 0

    result.update(
        {
            "cluster_file":  cluster_tsv,
            "rep_faa":       rep_faa if rep_faa.exists() else None,
            "n_clusters":    n_clusters,
            "membership_df": membership_df,
            "error":         "",
        }
    )
    return result


def _parse_mmseqs_cluster_tsv(tsv_path: Path) -> pd.DataFrame:
    """Parse MMseqs2 cluster TSV (rep_id\\tmember_id).

    Returns
    -------
    pd.DataFrame
        Columns: protein_id, cluster_id (integer), is_representative.
    """
    try:
        df = pd.read_csv(tsv_path, sep="\t", header=None, names=["rep_id", "member_id"])
    except Exception as exc:
        print(f"ERROR: Cannot parse MMseqs2 cluster TSV: {exc}", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_MEMBERSHIP_COLS)

    if df.empty:
        return pd.DataFrame(columns=_EMPTY_MEMBERSHIP_COLS)

    # Assign integer cluster IDs based on representative ordering
    rep_ids      = df["rep_id"].unique()
    rep_to_int   = {r: i for i, r in enumerate(rep_ids)}

    rows: list[dict] = []
    for _, row in df.iterrows():
        rep_id    = row["rep_id"]
        member_id = row["member_id"]
        rows.append(
            {
                "protein_id":       member_id,
                "cluster_id":       rep_to_int.get(rep_id, 0),
                "is_representative": member_id == rep_id,
            }
        )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

def cluster_dispatch(faa_path: Path, out_dir: Path, **kwargs) -> dict:
    """Try CD-HIT first, then MMseqs2, else return an empty stub.

    Parameters
    ----------
    faa_path : Path
        Input protein FASTA.
    out_dir : Path
        Output directory.
    **kwargs
        Passed to whichever clustering function is selected.
        Common keys: identity, coverage, threads.

    Returns
    -------
    dict
        Same structure as :func:`cluster_cdhit` / :func:`cluster_mmseqs`.
        ``error`` field contains the reason if neither tool is available.
    """
    faa_path = Path(faa_path)
    out_dir  = Path(out_dir)

    if find_tool("cd-hit"):
        print("INFO: Using cd-hit for clustering.", file=sys.stderr)
        result = cluster_cdhit(faa_path, out_dir, **kwargs)
        if result["error"] == "" and result["n_clusters"] > 0:
            return result
        print(
            f"WARNING: cd-hit clustering failed or returned 0 clusters: {result['error']}",
            file=sys.stderr,
        )

    if find_tool("mmseqs"):
        print("INFO: Using mmseqs for clustering.", file=sys.stderr)
        result = cluster_mmseqs(faa_path, out_dir, **kwargs)
        if result["error"] == "" and result["n_clusters"] > 0:
            return result
        print(
            f"WARNING: mmseqs clustering failed or returned 0 clusters: {result['error']}",
            file=sys.stderr,
        )

    print(
        "WARNING: Neither cd-hit nor mmseqs found on PATH; skipping clustering.",
        file=sys.stderr,
    )
    empty = _empty_result()
    empty["error"] = "No clustering tool available (cd-hit or mmseqs required)"
    return empty


# ---------------------------------------------------------------------------
# Cluster summary
# ---------------------------------------------------------------------------

def cluster_summary(membership_df: pd.DataFrame) -> pd.DataFrame:
    """Summarise clustering results per cluster.

    Parameters
    ----------
    membership_df : pd.DataFrame
        Columns: protein_id, cluster_id, is_representative.

    Returns
    -------
    pd.DataFrame
        Columns: cluster_id, size, representative_id.
        Sorted by size descending.  Empty DataFrame on failure.
    """
    if membership_df is None or membership_df.empty:
        return pd.DataFrame(columns=["cluster_id", "size", "representative_id"])

    required = {"protein_id", "cluster_id", "is_representative"}
    missing  = required - set(membership_df.columns)
    if missing:
        print(
            f"WARNING: cluster_summary: missing columns {missing}",
            file=sys.stderr,
        )
        return pd.DataFrame(columns=["cluster_id", "size", "representative_id"])

    # Size per cluster
    size_df = (
        membership_df
        .groupby("cluster_id")["protein_id"]
        .count()
        .reset_index()
        .rename(columns={"protein_id": "size"})
    )

    # Representative per cluster
    rep_df = (
        membership_df[membership_df["is_representative"] == True]  # noqa: E712
        .drop_duplicates(subset=["cluster_id"])
        [["cluster_id", "protein_id"]]
        .rename(columns={"protein_id": "representative_id"})
    )

    summary = size_df.merge(rep_df, on="cluster_id", how="left")
    summary = summary.sort_values("size", ascending=False).reset_index(drop=True)

    return summary
