"""
pipeline/annotation.py — TM topology + signal peptide + domain architecture.

Optional: enabled when phobius/tmhmm + Pfam domtblout available.

Provides:
  - Phobius wrapper (TM + signal peptide)
  - TMHMM wrapper (TM topology)
  - Pfam domain architecture from hmmsearch domtblout
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from .utils import find_tool, run_cmd
from typing import Optional

import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run_command(
    cmd: list[str],
    stdin_text: Optional[str] = None,
    timeout: int = 3600,
) -> tuple[int, str, str]:
    """Run a shell command; return (returncode, stdout, stderr)."""
    try:
        proc = subprocess.run(
            cmd,
            input=stdin_text,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        return proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired:
        return 1, "", "Command timed out"
    except FileNotFoundError:
        return 127, "", f"Binary not found: {cmd[0]}"
    except Exception as exc:
        return 1, "", str(exc)


_EMPTY_PHOBIUS_COLS = ["protein_id", "tm_count", "signal_peptide",
                       "predicted_class", "topology_string"]
_EMPTY_TMHMM_COLS  = ["protein_id", "tm_count", "topology_string"]


# ---------------------------------------------------------------------------
# Transmembrane / signal peptide prediction
# ---------------------------------------------------------------------------

def run_phobius(faa_path: Path, out_dir: Path) -> pd.DataFrame:
    """Run Phobius in short format to predict TM topology and signal peptides.

    Requires ``phobius.pl`` (or ``phobius``) on PATH.

    Parameters
    ----------
    faa_path : Path
        Input protein FASTA.
    out_dir : Path
        Output directory; ``phobius_output.txt`` is written here.

    Returns
    -------
    pd.DataFrame
        Columns: protein_id, tm_count, signal_peptide (Y/N), predicted_class,
        topology_string.
        predicted_class ∈ {TM_protein, SP_protein, SP+TM, Cytoplasmic}.
        Empty DataFrame if Phobius is unavailable.
    """
    faa_path = Path(faa_path)
    out_dir  = Path(out_dir)

    phobius_bin = find_tool("phobius.pl") or find_tool("phobius")
    if not phobius_bin:
        print("WARNING: phobius.pl not found on PATH; skipping TM/SP prediction.", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_PHOBIUS_COLS)

    if not faa_path.exists():
        print(f"ERROR: Input FASTA not found: {faa_path}", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_PHOBIUS_COLS)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "phobius_output.txt"

    cmd = [phobius_bin, "-short", str(faa_path)]
    print(f"INFO: Running: {' '.join(cmd)}", file=sys.stderr)
    rc, stdout, stderr = _run_command(cmd)

    if rc != 0:
        print(f"ERROR: Phobius failed (rc={rc}): {stderr[-500:]}", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_PHOBIUS_COLS)

    out_file.write_text(stdout)
    return _parse_phobius_short(stdout)


def _parse_phobius_short(text: str) -> pd.DataFrame:
    """Parse Phobius ``-short`` output.

    Format (space-separated):
        SEQID  TM  SP  PREDICTION

    where TM = number of TM segments, SP = Y/0, PREDICTION = text.
    """
    records = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("SEQNAME") or line.startswith("SEQID") or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) < 3:
            continue

        protein_id = parts[0]
        tm_count   = _safe_int(parts[1]) if len(parts) > 1 else 0
        sp         = parts[2].upper() if len(parts) > 2 else "0"
        prediction = parts[3] if len(parts) > 3 else ""

        signal_pep = "Y" if sp == "Y" else "N"

        if tm_count > 0 and signal_pep == "Y":
            cls = "SP+TM"
        elif tm_count > 0:
            cls = "TM_protein"
        elif signal_pep == "Y":
            cls = "SP_protein"
        else:
            cls = "Cytoplasmic"

        records.append({
            "protein_id":      protein_id,
            "tm_count":        tm_count,
            "signal_peptide":  signal_pep,
            "predicted_class": cls,
            "topology_string": prediction,
        })

    return pd.DataFrame(records) if records else pd.DataFrame(
        columns=_EMPTY_PHOBIUS_COLS
    )


def run_tmhmm(faa_path: Path, out_dir: Path) -> pd.DataFrame:
    """Run TMHMM2 in short output mode to predict TM topology.

    Parameters
    ----------
    faa_path : Path
        Input protein FASTA.
    out_dir : Path
        Output directory; ``tmhmm_out.txt`` is written here.

    Returns
    -------
    pd.DataFrame
        Columns: protein_id, tm_count, topology_string.
        Empty DataFrame if TMHMM is unavailable.
    """
    faa_path = Path(faa_path)
    out_dir  = Path(out_dir)

    tmhmm_bin = find_tool("tmhmm2") or find_tool("tmhmm")
    if not tmhmm_bin:
        print("WARNING: tmhmm not found on PATH; skipping TM prediction.", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_TMHMM_COLS)

    if not faa_path.exists():
        print(f"ERROR: Input FASTA not found: {faa_path}", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_TMHMM_COLS)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "tmhmm_out.txt"

    cmd = [tmhmm_bin, "--short", str(faa_path)]
    print(f"INFO: Running: {' '.join(cmd)}", file=sys.stderr)
    rc, stdout, stderr = _run_command(cmd)

    if rc != 0:
        print(f"ERROR: TMHMM failed (rc={rc}): {stderr[-500:]}", file=sys.stderr)
        return pd.DataFrame(columns=_EMPTY_TMHMM_COLS)

    out_file.write_text(stdout)
    return _parse_tmhmm_short(stdout)


def _parse_tmhmm_short(text: str) -> pd.DataFrame:
    """Parse TMHMM2 short output (tab- or space-separated).

    New format (tab-separated per-field):
        ID<TAB>len=N<TAB>ExpAA=N<TAB>First60=N<TAB>PredHel=N<TAB>Topology=...

    Old format (space-separated):
        # protein_id Length=N ExpAA=N First60=N PredHel=N Topology=...
    """
    records = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        # Skip pure comment lines that don't contain data
        if line.startswith("#") and "PredHel=" not in line:
            continue

        # Strip leading '#' if present (old format comment-data lines)
        data_line = line.lstrip("# ").strip()
        # Split on tabs or whitespace
        if "\t" in data_line:
            parts = data_line.split("\t")
        else:
            parts = re.split(r"\s+", data_line)

        if len(parts) < 2:
            continue

        protein_id = parts[0]
        tm_count   = 0
        topology   = ""

        for part in parts[1:]:
            if part.startswith("PredHel="):
                tm_count = _safe_int(part.split("=", 1)[1])
            elif part.startswith("Topology="):
                topology = part.split("=", 1)[1]

        records.append({
            "protein_id":      protein_id,
            "tm_count":        tm_count,
            "topology_string": topology,
        })

    return pd.DataFrame(records) if records else pd.DataFrame(
        columns=_EMPTY_TMHMM_COLS
    )


def annotate_from_tm(hits_df: pd.DataFrame, tm_df: pd.DataFrame) -> pd.DataFrame:
    """Merge TM/SP predictions into the main hits table.

    Handles both Phobius output (has ``signal_peptide``, ``predicted_class``)
    and TMHMM output (has ``topology_string``).

    Parameters
    ----------
    hits_df : pd.DataFrame
        Main hits table; must contain ``protein_id``.
    tm_df : pd.DataFrame
        Output of :func:`run_phobius` or :func:`run_tmhmm`.

    Returns
    -------
    pd.DataFrame
        hits_df with columns added/updated:
        tm_topology, signal_peptide, predicted_localization.
    """
    if hits_df is None or hits_df.empty:
        return hits_df

    if tm_df is None or tm_df.empty:
        for col in ("tm_topology", "signal_peptide", "predicted_localization"):
            if col not in hits_df.columns:
                hits_df = hits_df.copy()
                hits_df[col] = pd.NA
        return hits_df

    df = hits_df.copy()

    # Determine which columns to carry over from tm_df
    if "predicted_class" in tm_df.columns:
        merge_cols = [c for c in ["protein_id", "tm_count", "signal_peptide",
                                   "predicted_class", "topology_string"]
                      if c in tm_df.columns]
    else:
        merge_cols = [c for c in ["protein_id", "tm_count", "topology_string"]
                      if c in tm_df.columns]

    tm_sub = tm_df[merge_cols].drop_duplicates(subset=["protein_id"]).copy()
    tm_sub = tm_sub.rename(columns={
        "topology_string": "tm_topology",
        "predicted_class": "predicted_localization",
    })

    # Drop columns that will be replaced to avoid _x/_y suffixes
    for col in tm_sub.columns:
        if col != "protein_id" and col in df.columns:
            df.drop(columns=[col], inplace=True)

    df = df.merge(tm_sub, on="protein_id", how="left")

    # Ensure all three canonical columns exist
    if "tm_topology" not in df.columns:
        df["tm_topology"] = pd.NA
    if "signal_peptide" not in df.columns:
        df["signal_peptide"] = pd.NA
    if "predicted_localization" not in df.columns:
        df["predicted_localization"] = pd.NA

    return df


# ---------------------------------------------------------------------------
# Domain architecture
# ---------------------------------------------------------------------------

def domain_architecture(domtblout_path: Path, hmm_length: int = 0) -> pd.DataFrame:
    """Build domain architecture strings from a Pfam --domtblout result.

    For each protein_id, sorts hits by ali_from and concatenates domain names.

    Parameters
    ----------
    domtblout_path : Path
        Path to a Pfam hmmsearch ``--domtblout`` output file.
    hmm_length : int
        Length of the query HMM (kept for API consistency; not used internally).

    Returns
    -------
    pd.DataFrame
        Columns: protein_id, domain_architecture, n_domains.
        Architecture string: domains joined by ``|``, sorted by ali_from.
        E.g. ``"RRM_1|RdRp_4|Methyltransf_2"``.
    """
    if not domtblout_path or not Path(domtblout_path).exists():
        return pd.DataFrame(columns=["protein_id", "domain_architecture", "n_domains"])

    records = []
    try:
        with open(domtblout_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 23:
                    continue
                try:
                    records.append({
                        "protein_id":  parts[0],
                        "domain_name": parts[3].split(".")[0],  # strip clan suffix
                        "ali_from":    int(parts[17]),
                        "ali_to":      int(parts[18]),
                        "evalue":      float(parts[11]),
                    })
                except (IndexError, ValueError):
                    continue
    except Exception as exc:
        print(f"[annotation] domtblout parse error: {exc}", file=sys.stderr)
        return pd.DataFrame(columns=["protein_id", "domain_architecture", "n_domains"])

    if not records:
        return pd.DataFrame(columns=["protein_id", "domain_architecture", "n_domains"])

    df = pd.DataFrame(records)

    # Keep best e-value per (protein, domain) pair to avoid duplicates
    df = df.sort_values("evalue").drop_duplicates(subset=["protein_id", "domain_name"])
    df = df.sort_values(["protein_id", "ali_from"])

    try:
        result = (
            df.groupby("protein_id", sort=False)
            .apply(lambda g: "|".join(g["domain_name"].tolist()), include_groups=False)
            .reset_index()
        )
    except TypeError:
        # Older pandas without include_groups
        result = (
            df.groupby("protein_id", sort=False)
            .apply(lambda g: "|".join(g["domain_name"].tolist()))
            .reset_index()
        )

    result.columns = ["protein_id", "domain_architecture"]
    result["n_domains"] = result["domain_architecture"].apply(lambda x: len(x.split("|")))
    return result


# ---------------------------------------------------------------------------
# Numeric helper
# ---------------------------------------------------------------------------

def _safe_int(value: str, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except (ValueError, TypeError):
        return default
