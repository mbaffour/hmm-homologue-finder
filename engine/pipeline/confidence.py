"""
pipeline/confidence.py — Multi-evidence hit confidence scoring.

Confidence tiers:
  high_confidence  — bit >= strict AND hmm_cov >= 0.60 AND (reciprocal OR domain_match)
  putative         — bit >= strict AND hmm_cov >= 0.30
  divergent        — bit >= moderate AND (hmm_cov < 0.30 OR domain match present)
  likely_fp        — high bias, low complexity, or fails all evidence

QC flags (pipe-separated string):
  high_bias        — bias_score > 5.0
  short_alignment  — hmm_cov < 0.50
  low_complexity   — sequence contains >40% single amino acid
  contig_edge      — protein < 80 aa truncated at sequence edge, OR hit within
                     30 aa of real contig boundary (requires contig_length column)
"""
from __future__ import annotations

import sys
from typing import Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def score_hit(
    row: pd.Series,
    hmm_length: int,
    strict: float,
    moderate: float,
    hmm_cov_floor: float = 0.30,
) -> Tuple[str, str]:
    """Classify a single search hit and return (confidence_tier, why_classified).

    Parameters
    ----------
    row : pd.Series
        Must contain: bit_score, bias_score, hmm_from, hmm_to.
        Optional: reciprocal_hit (bool), domain_match (bool),
                  taxonomy_outlier (bool).
    hmm_length : int
        Total length of the HMM in match states.
    strict : float
        Strict bit-score threshold (high-confidence boundary).
    moderate : float
        Moderate bit-score threshold (putative / divergent boundary).
    hmm_cov_floor : float
        Minimum HMM coverage fraction to count as "enough" alignment.

    Returns
    -------
    tuple[str, str]
        (confidence_tier, why_classified)
    """
    if hmm_length <= 0:
        return "likely_fp", "HMM length is zero or negative"

    import pandas as _pd

    def _safe_float(val, default=0.0):
        return default if (val is None or (_pd.isna(val) if not isinstance(val, str) else False)) else float(val)

    def _safe_int(val, default=0):
        return default if (val is None or (_pd.isna(val) if not isinstance(val, str) else False)) else int(val)

    def _safe_bool(val, default=False):
        return default if (val is None or (_pd.isna(val) if not isinstance(val, str) else False)) else bool(val)

    bit        = _safe_float(row.get("bit_score",  0.0))
    bias       = _safe_float(row.get("bias_score", 0.0))
    domain_ok  = _safe_bool(row.get("domain_match",    False))

    # Reciprocal validation: three states
    #   True  — ran and confirmed (boosts confidence)
    #   False AND column exists and non-null — ran and FAILED (demotes)
    #   None / NA  — not yet run; classify on bit+coverage alone
    _recip_raw = row.get("reciprocal_hit")
    _recip_is_na = _recip_raw is None or (
        not isinstance(_recip_raw, (bool, int)) and _pd.isna(_recip_raw)
    )
    reciprocal_run    = not _recip_is_na          # was the step actually executed?
    reciprocal        = bool(_recip_raw) if reciprocal_run else False

    # Prefer pre-computed hmm_coverage_pct (populated from domtblout).
    # Fall back to hmm_from/hmm_to if available.
    # If none of these are populated, assume full coverage (1.0) so that
    # per-sequence tblout-only searches are not penalised for missing domain data.
    if not _pd.isna(row.get("hmm_coverage_pct")) and row.get("hmm_coverage_pct") is not None:
        hmm_cov = _safe_float(row.get("hmm_coverage_pct")) / 100.0
    else:
        hmm_from = _safe_int(row.get("hmm_from", 0))
        hmm_to   = _safe_int(row.get("hmm_to",   0))
        if hmm_from > 0 and hmm_to > 0 and hmm_length > 0:
            hmm_cov = (hmm_to - hmm_from + 1) / hmm_length
        else:
            # No domain-level coordinates available; assume full coverage
            hmm_cov = 1.0

    reasons: list[str] = [
        f"bit={bit:.1f}",
        f"hmm_cov={hmm_cov:.2f}",
    ]
    if reciprocal_run:
        reasons.append("reciprocal=confirmed" if reciprocal else "reciprocal=failed")
    if domain_ok:
        reasons.append("domain_match=yes")
    if bias > 5.0:
        reasons.append(f"bias={bias:.1f}(high)")

    # ---- Tier assignment ----
    # Likely FP: high bias overwhelming the score
    if bias > bit * 0.80 and bit < strict:
        return "likely_fp", "; ".join(reasons + ["bias_dominates_score"])

    if bit >= strict:
        if hmm_cov >= 0.60:
            # If reciprocal was run and FAILED, demote to putative
            if reciprocal_run and not reciprocal and not domain_ok:
                return "putative", "; ".join(reasons + ["reciprocal_failed"])
            # Strong bit + good coverage = high_confidence
            # (reciprocal confirms it; domain_match also confirms; neither required)
            return "high_confidence", "; ".join(reasons)
        if hmm_cov >= hmm_cov_floor:
            return "putative", "; ".join(reasons + ["coverage_below_60pct"])
        return "divergent", "; ".join(reasons + ["low_hmm_coverage"])

    if bit >= moderate:
        return "divergent", "; ".join(reasons + ["below_strict_threshold"])

    return "likely_fp", "; ".join(reasons + ["below_moderate_threshold"])


def classify_hits(
    df: pd.DataFrame,
    hmm_length: int,
    strict: float = 45.0,
    moderate: float = 30.0,
    hmm_cov_floor: float = 0.30,
) -> pd.DataFrame:
    """Apply :func:`score_hit` to every row; add classification columns.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns: bit_score, bias_score, hmm_from, hmm_to.
    hmm_length : int
        HMM profile length in match states.
    strict : float
        Strict bit-score threshold.
    moderate : float
        Moderate bit-score threshold.
    hmm_cov_floor : float
        Minimum HMM coverage for putative tier.

    Returns
    -------
    pd.DataFrame
        Input df with added columns:
        hmm_coverage_pct, confidence_tier, why_classified, qc_flags.
    """
    if df.empty:
        for col in ("hmm_coverage_pct", "confidence_tier", "why_classified", "qc_flags"):
            df[col] = pd.Series(dtype="object")
        return df

    df = df.copy()

    # Ensure required columns exist with sensible defaults.
    # reciprocal_hit deliberately left as pd.NA (not False) so that
    # score_hit() can distinguish "not run" from "ran and failed".
    for col, default in [
        ("bit_score",      0.0),
        ("bias_score",     0.0),
        ("hmm_from",       0),
        ("hmm_to",         0),
        ("domain_match",   False),
    ]:
        if col not in df.columns:
            df[col] = default

    # Compute HMM coverage percentage.
    # Only populate from hmm_from/hmm_to when they contain real data (> 0).
    # When both are 0 (tblout-only search — no domain coords), fall back to
    # 100% so that full-length protein hits are not penalised for missing data.
    if hmm_length > 0:
        hmm_from_s = df["hmm_from"].fillna(0).astype(int)
        hmm_to_s   = df["hmm_to"].fillna(0).astype(int)
        has_domain_coords = (hmm_from_s > 0) | (hmm_to_s > 0)
        raw_cov = ((hmm_to_s - hmm_from_s + 1) / hmm_length * 100.0).round(2)
        df["hmm_coverage_pct"] = raw_cov.where(has_domain_coords, other=100.0)
    else:
        df["hmm_coverage_pct"] = 100.0

    # Score each hit
    tiers: list[str] = []
    whys:  list[str] = []
    for _, row in df.iterrows():
        tier, why = score_hit(row, hmm_length, strict, moderate, hmm_cov_floor)
        tiers.append(tier)
        whys.append(why)

    df["confidence_tier"] = tiers
    df["why_classified"]  = whys

    # QC flags
    df["qc_flags"] = add_qc_flags(df)

    return df


def add_qc_flags(
    df: pd.DataFrame,
    bias_threshold: float = 5.0,
    cov_threshold: float = 0.50,
) -> pd.Series:
    """Compute pipe-separated QC flag strings for each row.

    Parameters
    ----------
    df : pd.DataFrame
        Hit table; expects bias_score, hmm_coverage_pct (0-100 scale),
        and optionally seq_length, seq_from, seq_to, sequence columns.
    bias_threshold : float
        bias_score above this value → ``high_bias`` flag.
    cov_threshold : float
        HMM coverage fraction below this value → ``short_alignment`` flag.

    Returns
    -------
    pd.Series
        String Series of pipe-separated flags, one entry per row.
    """
    if df.empty:
        return pd.Series(dtype="str")

    df = df.copy()

    flags_list: list[str] = []

    # Compute HMM coverage fraction (0–1).
    # When hmm_coverage_pct is present and non-null, use it directly.
    # When hmm_from/hmm_to are present and > 0, compute from them.
    # Otherwise default to 1.0 (full coverage assumed — tblout-only search).
    if "hmm_coverage_pct" in df.columns and df["hmm_coverage_pct"].notna().any():
        hmm_cov_frac = df["hmm_coverage_pct"].fillna(100.0) / 100.0
    elif "hmm_from" in df.columns and "hmm_to" in df.columns:
        # hmm_length not available here; skip fractional coverage, default to 1.0
        hmm_cov_frac = pd.Series(1.0, index=df.index)
    else:
        hmm_cov_frac = pd.Series(1.0, index=df.index)

    bias   = df["bias_score"].fillna(0.0) if "bias_score" in df.columns else pd.Series(0.0, index=df.index)
    seq_len = df.get("seq_length", pd.Series(0, index=df.index)).fillna(0).astype(int)
    seq_from = df.get("seq_from",  pd.Series(0, index=df.index)).fillna(0).astype(int)
    seq_to   = df.get("seq_to",    pd.Series(0, index=df.index)).fillna(0).astype(int)

    for idx in df.index:
        row_flags: list[str] = []

        if bias.at[idx] > bias_threshold:
            row_flags.append("high_bias")

        if hmm_cov_frac.at[idx] < cov_threshold:
            row_flags.append("short_alignment")

        # Low complexity: single AA > 40% of aligned region
        if "sequence" in df.columns:
            seq_str = str(df.at[idx, "sequence"] or "")
            if seq_str:
                char_counts = {c: seq_str.count(c) for c in set(seq_str) if c not in ("-", "*")}
                if char_counts:
                    max_frac = max(char_counts.values()) / len(seq_str)
                    if max_frac > 0.40:
                        row_flags.append("low_complexity")

        # Contig edge: flag only when the protein looks genuinely truncated.
        # A protein is a potential contig edge hit when it is short (< 80 aa)
        # AND the alignment butts against one end of the sequence (within 5 aa
        # of position 1 or of seq_length).  This avoids false positives on
        # full-length proteins where the HMM simply happens to start near
        # position 1 of the target (which is nearly every protein).
        # The flag is also respected when an explicit "contig_length" column is
        # present (populated from genome/nucleotide coordinate data).
        s_len  = seq_len.at[idx]
        s_from = seq_from.at[idx]
        s_to   = seq_to.at[idx]
        _contig_len = int(df.at[idx, "contig_length"]) if "contig_length" in df.columns else 0
        if _contig_len > 0:
            # Nucleotide-derived coordinates: real contig edge check
            if s_from <= 30 or (_contig_len - s_to) <= 30:
                row_flags.append("contig_edge")
        elif s_len > 0 and s_len < 80 and (s_from <= 5 or (s_len - s_to) <= 5):
            # Short protein truncated at sequence edge → likely metagenome fragment
            row_flags.append("contig_edge")

        flags_list.append("|".join(row_flags) if row_flags else "")

    return pd.Series(flags_list, index=df.index)
