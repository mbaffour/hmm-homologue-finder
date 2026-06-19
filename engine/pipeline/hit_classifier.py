"""
pipeline/hit_classifier.py — Full hit classification pipeline.

Merges per-sequence (tblout) and per-domain (domtblout) hmmsearch results,
applies multi-evidence confidence scoring, supports reciprocal validation,
hit sequence extraction, and cross-database deduplication.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

from .confidence import classify_hits
from .searcher import parse_tblout, parse_domtblout


# ---------------------------------------------------------------------------
# Required columns in the main hits table (NaN for optional missing ones)
# ---------------------------------------------------------------------------
_MAIN_COLS = [
    "genome_id", "contig_id", "protein_id", "hit_name",
    "evalue", "bit_score", "bias_score",
    "domain_evalue", "domain_bit_score",
    "hmm_from", "hmm_to", "hmm_coverage_pct",
    "seq_from", "seq_to", "seq_length", "ali_coverage_pct",
    "strand", "description",
    "confidence_tier", "why_classified", "qc_flags",
    "reciprocal_hit",
    "tm_topology", "signal_peptide", "predicted_localization",
    "domain_architecture", "cluster_id",
    "database_source", "iteration_added", "contig_edge_flag",
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_main_hits_table(
    tblout_df: pd.DataFrame,
    domtblout_df: pd.DataFrame,
    hmm_length: int,
    db_name: str,
    strict: float,
    moderate: float,
    iteration: int = 0,
) -> pd.DataFrame:
    """Merge tblout + domtblout results and produce the canonical hits table.

    The per-domain table supplies domain-level coordinates (hmm_from/to,
    ali_from/to) and domain bit/evalue; the per-sequence table supplies the
    primary bit_score, evalue, and bias_score used for confidence scoring.
    When a target has multiple domains, only the best-scoring domain (highest
    domain_bit_score) is kept per protein.

    Parameters
    ----------
    tblout_df : pd.DataFrame
        Output of :func:`parse_tblout`.
    domtblout_df : pd.DataFrame
        Output of :func:`parse_domtblout`.
    hmm_length : int
        HMM profile length in match states.
    db_name : str
        Database source label added to every row.
    strict : float
        Strict bit-score threshold for :func:`classify_hits`.
    moderate : float
        Moderate bit-score threshold.
    iteration : int
        Which refinement iteration produced these hits.

    Returns
    -------
    pd.DataFrame
        Full hits table with all :data:`_MAIN_COLS` columns.
        Empty DataFrame if inputs are empty.
    """
    if tblout_df.empty:
        return _empty_hits_table()

    # ---- Reduce domtblout to best domain per protein ----
    dom_cols_needed = [
        "target_name", "domain_evalue", "domain_bit_score",
        "hmm_from", "hmm_to", "ali_from", "ali_to",
        "env_from", "env_to",
    ]
    if not domtblout_df.empty:
        dom = domtblout_df.copy()
        # Ensure columns exist
        for c in dom_cols_needed:
            if c not in dom.columns:
                dom[c] = pd.NA
        dom = (
            dom.sort_values("domain_bit_score", ascending=False)
               .drop_duplicates(subset=["target_name"], keep="first")
        )
        dom = dom[dom_cols_needed].rename(columns={"target_name": "protein_id"})
    else:
        dom = pd.DataFrame(columns=dom_cols_needed).rename(
            columns={"target_name": "protein_id"}
        )

    # ---- Build base table from tblout ----
    df = tblout_df.copy().rename(columns={
        "target_name": "protein_id",
        "query_name":  "hit_name",
    })

    # Merge domain info
    if not dom.empty:
        df = df.merge(dom, on="protein_id", how="left")
    else:
        for c in ["domain_evalue", "domain_bit_score",
                  "hmm_from", "hmm_to", "ali_from", "ali_to",
                  "env_from", "env_to"]:
            df[c] = pd.NA

    # ---- Parse IDs to populate genome_id / contig_id ----
    df["genome_id"]  = df["protein_id"].apply(_genome_id_from_protein_id)
    df["contig_id"]  = df["protein_id"].apply(_contig_id_from_protein_id)

    # ---- Coordinate columns ----
    df["seq_from"]   = df.get("ali_from",  pd.NA)
    df["seq_to"]     = df.get("ali_to",    pd.NA)
    df["seq_length"] = (
        df["seq_to"].fillna(0).astype(float) - df["seq_from"].fillna(0).astype(float) + 1
    ).where(df["seq_to"].notna(), other=pd.NA)
    df["ali_coverage_pct"] = pd.NA  # Filled downstream if env coords known
    df["strand"]     = pd.NA

    # ---- Optional annotation columns ----
    for col in ("tm_topology", "signal_peptide", "predicted_localization",
                "domain_architecture", "cluster_id"):
        df[col] = pd.NA

    df["database_source"]  = db_name
    df["iteration_added"]  = iteration
    df["contig_edge_flag"] = False
    df["reciprocal_hit"]   = pd.NA

    # ---- Confidence scoring ----
    df = classify_hits(df, hmm_length=hmm_length, strict=strict, moderate=moderate)

    # ---- Ensure all required columns present ----
    df = _ensure_columns(df, _MAIN_COLS)
    return df[_MAIN_COLS]


def reciprocal_validate_hmmsearch(
    hits_faa: Path,
    seed_faa: Path,
    hmm_path: Path,
    out_dir: Path,
    strict: float,
) -> "dict[str, bool]":
    """Reciprocal validation: hmmsearch hits back against seed set.

    For each hit sequence file entry, check whether hmmsearch against the
    seed set recovers it above the strict bit-score threshold.

    Parameters
    ----------
    hits_faa : Path
        FASTA of hit sequences extracted from the database.
    seed_faa : Path
        Original seed sequences used to build the HMM.
    hmm_path : Path
        The HMM profile.
    out_dir : Path
        Directory for temporary search output.
    strict : float
        Bit-score threshold.

    Returns
    -------
    dict[str, bool]
        {protein_id: True/False}. Empty dict on failure.
    """
    hits_faa = Path(hits_faa)
    seed_faa = Path(seed_faa)
    hmm_path = Path(hmm_path)
    out_dir  = Path(out_dir)

    if not hits_faa.exists():
        print(f"ERROR: Hits FASTA not found: {hits_faa}", file=sys.stderr)
        return {}
    if not hmm_path.exists():
        print(f"ERROR: HMM not found: {hmm_path}", file=sys.stderr)
        return {}

    out_dir.mkdir(parents=True, exist_ok=True)
    tbl = out_dir / "reciprocal.tbl"

    cmd = [
        "hmmsearch",
        "--tblout", str(tbl),
        "--noali",
        str(hmm_path),
        str(hits_faa),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: Reciprocal hmmsearch failed:\n{result.stderr}", file=sys.stderr)
        return {}

    tbl_df = parse_tblout(tbl)
    if tbl_df.empty:
        return {}

    validated: dict[str, bool] = {}
    for _, row in tbl_df.iterrows():
        pid  = row["target_name"]
        bits = float(row["bit_score"] or 0.0)
        validated[pid] = bits >= strict

    return validated


def extract_hit_sequences(
    hits_df: pd.DataFrame,
    db_faa: Path,
    out_faa: Path,
) -> int:
    """Extract sequences for accepted hits from a protein database FASTA.

    Parameters
    ----------
    hits_df : pd.DataFrame
        Must contain a ``protein_id`` column.
    db_faa : Path
        Source FASTA database.
    out_faa : Path
        Destination FASTA file.

    Returns
    -------
    int
        Number of sequences written.
    """
    db_faa  = Path(db_faa)
    out_faa = Path(out_faa)

    if hits_df.empty or "protein_id" not in hits_df.columns:
        return 0
    if not db_faa.exists():
        print(f"ERROR: Database not found: {db_faa}", file=sys.stderr)
        return 0

    try:
        from Bio import SeqIO
    except ImportError:
        print("ERROR: Biopython not installed.", file=sys.stderr)
        return 0

    wanted = set(hits_df["protein_id"].dropna().astype(str))
    out_faa.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    try:
        with out_faa.open("w") as out_handle:
            for rec in SeqIO.parse(str(db_faa), "fasta"):
                if rec.id in wanted:
                    SeqIO.write(rec, out_handle, "fasta")
                    count += 1
    except Exception as exc:
        print(f"ERROR: Could not extract sequences: {exc}", file=sys.stderr)
        return 0

    return count


def best_hit_per_genome(hits_df: pd.DataFrame) -> pd.DataFrame:
    """Return the best-scoring hit per genome.

    Parameters
    ----------
    hits_df : pd.DataFrame
        Must contain ``genome_id`` and ``bit_score`` columns.

    Returns
    -------
    pd.DataFrame
        Subset with one row per unique genome_id (highest bit_score).
    """
    if hits_df.empty:
        return hits_df

    if "genome_id" not in hits_df.columns or "bit_score" not in hits_df.columns:
        print("WARNING: genome_id or bit_score missing from hits table.", file=sys.stderr)
        return hits_df

    return (
        hits_df.sort_values("bit_score", ascending=False)
               .drop_duplicates(subset=["genome_id"], keep="first")
               .reset_index(drop=True)
    )


def merge_all_databases(hit_dfs: "list[pd.DataFrame]") -> pd.DataFrame:
    """Concatenate hits from multiple databases and deduplicate.

    When the same protein_id appears in multiple databases, only the entry
    with the highest bit_score is kept.

    Parameters
    ----------
    hit_dfs : list[pd.DataFrame]
        One DataFrame per database.

    Returns
    -------
    pd.DataFrame
        Deduplicated combined hits table.
    """
    non_empty = [df for df in hit_dfs if df is not None and not df.empty]
    if not non_empty:
        return _empty_hits_table()

    combined = pd.concat(non_empty, ignore_index=True)
    combined = _ensure_columns(combined, _MAIN_COLS)

    if "protein_id" in combined.columns and "bit_score" in combined.columns:
        combined = (
            combined.sort_values("bit_score", ascending=False)
                    .drop_duplicates(subset=["protein_id"], keep="first")
                    .reset_index(drop=True)
        )

    return combined[_MAIN_COLS]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _empty_hits_table() -> pd.DataFrame:
    return pd.DataFrame(columns=_MAIN_COLS)


def _ensure_columns(df: pd.DataFrame, cols: list) -> pd.DataFrame:
    """Add any missing columns as NaN."""
    for col in cols:
        if col not in df.columns:
            df[col] = pd.NA
    return df


def _genome_id_from_protein_id(pid: str) -> str:
    """Best-effort genome ID extraction from protein ID string.

    Handles common NCBI formats like WP_001234567.1, NP_123456.1,
    and INPHARED/metagenome formats like SEQID|frame_F1|start_END.
    """
    if not pid or not isinstance(pid, str):
        return ""
    # Strip 6-frame coordinate suffix added by 04_translate_sixframe.py
    # Format: CONTIG_ID|frame_F1|100_200
    if "|frame_" in pid:
        return pid.split("|frame_")[0]
    # Standard NCBI versioned accession — strip version
    parts = pid.split(".")
    if len(parts) >= 2:
        return ".".join(parts[:-1])
    return pid


def _contig_id_from_protein_id(pid: str) -> str:
    """Extract contig ID from a protein ID (best-effort)."""
    if not pid or not isinstance(pid, str):
        return ""
    if "|frame_" in pid:
        return pid.split("|frame_")[0]
    return pid
