"""
pipeline/iterative.py — Iterative HMM refinement (jackhmmer-style).

Manages the logic for expanding an HMM through successive rounds of searching,
checking for convergence, and growing the seed set with approved hits.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def iteration_candidates(
    hits_df: pd.DataFrame,
    seeds_faa: Path,
    strict: float,
) -> pd.DataFrame:
    """Return high-confidence hits that are not already in the seed set.

    Candidates are used as candidates for adding to the next iteration's
    seed alignment. The caller reviews them before committing.

    Parameters
    ----------
    hits_df : pd.DataFrame
        Full hits table from the current iteration.
    seeds_faa : Path
        Current seed FASTA file (used to exclude already-present IDs).
    strict : float
        Minimum bit score to accept as a candidate.

    Returns
    -------
    pd.DataFrame
        Subset of hits_df with columns:
        protein_id, bit_score, evalue, hmm_coverage_pct, confidence_tier,
        database_source, description. Sorted by bit_score descending.
    """
    seeds_faa = Path(seeds_faa)

    keep_cols = [
        "protein_id", "bit_score", "evalue", "hmm_coverage_pct",
        "confidence_tier", "database_source", "description",
    ]

    if hits_df.empty:
        return pd.DataFrame(columns=keep_cols)

    # Load existing seed IDs
    seed_ids: set[str] = set()
    if seeds_faa.exists():
        try:
            from Bio import SeqIO
            for rec in SeqIO.parse(str(seeds_faa), "fasta"):
                seed_ids.add(rec.id)
        except Exception as exc:
            print(f"WARNING: Could not read seed FASTA {seeds_faa}: {exc}", file=sys.stderr)

    # Filter to high-confidence hits above strict threshold
    mask = (
        hits_df.get("confidence_tier", pd.Series(dtype=str)) == "high_confidence"
    ) & (
        hits_df.get("bit_score", pd.Series(dtype=float)).fillna(0.0) >= strict
    )
    candidates = hits_df[mask].copy()

    # Exclude sequences already in seeds
    if "protein_id" in candidates.columns and seed_ids:
        candidates = candidates[~candidates["protein_id"].isin(seed_ids)]

    # Ensure expected columns
    for col in keep_cols:
        if col not in candidates.columns:
            candidates[col] = pd.NA

    return (
        candidates[keep_cols]
        .sort_values("bit_score", ascending=False)
        .reset_index(drop=True)
    )


def convergence_check(
    prev_count: int,
    curr_count: int,
    prev_leng: int,
    curr_leng: int,
) -> bool:
    """Determine whether iterative refinement has converged.

    Convergence requires BOTH conditions:
      - Hit count change < 5 % relative to previous count.
      - HMM LENG change < 3 match-state positions.

    Parameters
    ----------
    prev_count : int
        Hit count at the previous iteration.
    curr_count : int
        Hit count at the current iteration.
    prev_leng : int
        HMM length (LENG) at the previous iteration.
    curr_leng : int
        HMM length at the current iteration.

    Returns
    -------
    bool
        True if converged.
    """
    if prev_count == 0 and curr_count == 0:
        return True

    # Avoid division by zero
    count_change = abs(curr_count - prev_count) / max(prev_count, 1)
    leng_change  = abs(curr_leng - prev_leng)

    count_converged = count_change < 0.05
    leng_converged  = leng_change < 3

    return count_converged and leng_converged


def convergence_data(history: "list[dict]") -> dict:
    """Prepare convergence history for a Plotly line chart.

    Parameters
    ----------
    history : list[dict]
        Each dict has keys: iteration, hit_count, hmm_leng, diversity.

    Returns
    -------
    dict
        {iterations, hit_counts, leng_values, diversity_values}
        suitable for plotly go.Scatter / go.Figure construction.
    """
    if not history:
        return {
            "iterations":       [],
            "hit_counts":       [],
            "leng_values":      [],
            "diversity_values": [],
        }

    iterations       = [d.get("iteration",  i) for i, d in enumerate(history)]
    hit_counts       = [d.get("hit_count",  0) for d in history]
    leng_values      = [d.get("hmm_leng",   0) for d in history]
    diversity_values = [d.get("diversity",  0) for d in history]

    return {
        "iterations":       iterations,
        "hit_counts":       hit_counts,
        "leng_values":      leng_values,
        "diversity_values": diversity_values,
    }


def append_to_seeds(
    existing_faa: Path,
    new_seqs_faa: Path,
    approved_ids: "list[str]",
    out_faa: Path,
) -> int:
    """Concatenate approved sequences from new_seqs_faa into the seed set.

    Parameters
    ----------
    existing_faa : Path
        Current seed FASTA.
    new_seqs_faa : Path
        FASTA containing candidate sequences to append.
    approved_ids : list[str]
        Subset of protein IDs from new_seqs_faa to include.
    out_faa : Path
        Output path for the expanded seed FASTA.

    Returns
    -------
    int
        Number of sequences added (not counting pre-existing seeds).
    """
    existing_faa = Path(existing_faa)
    new_seqs_faa = Path(new_seqs_faa)
    out_faa      = Path(out_faa)

    try:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
    except ImportError:
        print("ERROR: Biopython not installed.", file=sys.stderr)
        return 0

    # Load existing seeds
    existing_records: list[SeqRecord] = []
    existing_ids: set[str] = set()
    if existing_faa.exists():
        try:
            existing_records = list(SeqIO.parse(str(existing_faa), "fasta"))
            existing_ids = {r.id for r in existing_records}
        except Exception as exc:
            print(f"WARNING: Could not parse existing seeds {existing_faa}: {exc}",
                  file=sys.stderr)

    # Load approved new sequences
    approved_set = set(approved_ids)
    new_records: list[SeqRecord] = []
    if new_seqs_faa.exists() and approved_set:
        try:
            for rec in SeqIO.parse(str(new_seqs_faa), "fasta"):
                if rec.id in approved_set and rec.id not in existing_ids:
                    new_records.append(rec)
        except Exception as exc:
            print(f"WARNING: Could not parse candidate FASTA {new_seqs_faa}: {exc}",
                  file=sys.stderr)

    # Write combined output
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    try:
        all_records = existing_records + new_records
        SeqIO.write(all_records, str(out_faa), "fasta")
    except Exception as exc:
        print(f"ERROR: Could not write {out_faa}: {exc}", file=sys.stderr)
        return 0

    return len(new_records)
