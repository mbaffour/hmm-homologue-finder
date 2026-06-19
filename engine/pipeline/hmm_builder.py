"""
pipeline/hmm_builder.py — hmmbuild wrapper + HMM build report + logo data.

Provides functions to build profile HMMs from alignments, parse .hmm files,
extract logo data from match emission probabilities, and validate HMMs by
self-searching against the seed sequences.
"""
import re
import subprocess
import sys
from pathlib import Path
from typing import Optional

from .utils import find_tool, run_cmd


# HMMER amino acid ordering in .hmm emission probability lines
_HMMER_AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_hmmbuild(
    aln_path: Path,
    hmm_path: Path,
    hmm_name: str = "novel_phage_gene",
) -> dict:
    """Build a profile HMM from a multiple sequence alignment.

    Parameters
    ----------
    aln_path : Path
        Trimmed, aligned FASTA input.
    hmm_path : Path
        Destination .hmm file.
    hmm_name : str
        Name embedded in the HMM (``--name``).

    Returns
    -------
    dict
        {leng, alph, nseq, cksum, hmmbuild_version} or empty dict on failure.
    """
    aln_path = Path(aln_path)
    hmm_path = Path(hmm_path)

    if not aln_path.exists():
        print(f"ERROR: Alignment not found: {aln_path}", file=sys.stderr)
        return {}

    hmm_path.parent.mkdir(parents=True, exist_ok=True)

    hmmbuild_bin = find_tool("hmmbuild") or "hmmbuild"
    cmd = [
        hmmbuild_bin,
        "-n", hmm_name,
        str(hmm_path),
        str(aln_path),
    ]

    result = run_cmd(cmd)
    if result.returncode != 0:
        print(f"ERROR: hmmbuild failed:\n{result.stderr}", file=sys.stderr)
        return {}

    if not hmm_path.exists():
        print("ERROR: hmmbuild did not produce output file.", file=sys.stderr)
        return {}

    # Parse stdout for key stats
    parsed = _parse_hmmbuild_stdout(result.stdout)

    # Also parse the .hmm file for cross-verification
    file_parsed = parse_hmm_file(hmm_path)

    return {
        "leng": file_parsed.get("LENG", parsed.get("leng", 0)),
        "alph": file_parsed.get("ALPH", parsed.get("alph", "amino")),
        "nseq": file_parsed.get("NSEQ", parsed.get("nseq", 0)),
        "cksum": file_parsed.get("CKSUM", parsed.get("cksum", "")),
        "hmmbuild_version": parsed.get("version", ""),
        "name": file_parsed.get("NAME", hmm_name),
    }


def parse_hmm_file(hmm_path: Path) -> dict:
    """Read a .hmm text file and extract header metadata.

    Parameters
    ----------
    hmm_path : Path
        Path to a HMMER3 .hmm file.

    Returns
    -------
    dict
        {NAME, LENG, ALPH, NSEQ, DATE, CKSUM} — values are cast to int
        where appropriate. Empty dict on failure.
    """
    hmm_path = Path(hmm_path)
    if not hmm_path.exists():
        print(f"ERROR: HMM file not found: {hmm_path}", file=sys.stderr)
        return {}

    result: dict = {}
    int_fields = {"LENG", "NSEQ"}
    target_fields = {"NAME", "LENG", "ALPH", "NSEQ", "DATE", "CKSUM"}

    try:
        with hmm_path.open() as fh:
            for line in fh:
                parts = line.split(None, 1)
                if len(parts) < 2:
                    continue
                key = parts[0].strip()
                if key in target_fields:
                    val = parts[1].strip()
                    result[key] = int(val) if key in int_fields else val
                # Stop reading once we hit the model body
                if line.startswith("HMM "):
                    break
    except Exception as exc:
        print(f"ERROR: Cannot parse {hmm_path}: {exc}", file=sys.stderr)
        return {}

    return result


def logo_data(hmm_path: Path, sample_every: int = 1) -> "list[dict]":
    """Extract match emission probabilities from a .hmm file for logo rendering.

    The HMMER3 .hmm body format has 3 lines per node:
      Line 1 — match emissions (20 aa log-odds)
      Line 2 — insert emissions
      Line 3 — state transitions

    Parameters
    ----------
    hmm_path : Path
        Path to the .hmm file.
    sample_every : int
        Sample every N-th position (1 = all positions).

    Returns
    -------
    list[dict]
        One dict per position: {pos, top_aa (list of {aa, prob} sorted desc, top 5)}.
        Empty list on failure.
    """
    hmm_path = Path(hmm_path)
    if not hmm_path.exists():
        print(f"ERROR: HMM file not found: {hmm_path}", file=sys.stderr)
        return []

    positions: list = []

    try:
        text = hmm_path.read_text()
    except Exception as exc:
        print(f"ERROR: Cannot read {hmm_path}: {exc}", file=sys.stderr)
        return []

    # Locate the HMM body — starts after the "HMM " header line and the
    # two annotation lines (aa list and transitions list)
    in_body = False
    skip_header_lines = 0
    lines = text.splitlines()
    line_idx = 0

    while line_idx < len(lines):
        line = lines[line_idx]
        if re.match(r"^HMM\s+", line):
            in_body = True
            skip_header_lines = 2  # skip the aa-order and transition-name lines
            line_idx += 1
            continue
        if in_body and skip_header_lines > 0:
            skip_header_lines -= 1
            line_idx += 1
            continue
        if in_body:
            # Try to match a position line: starts with integer node number
            # Format: "     1   <20 floats>   <extra>"
            stripped = line.strip()
            if stripped.startswith("//"):
                break
            match = re.match(r"^\s*(\d+)\s+([\d.*\s]+)", stripped)
            if match:
                pos = int(match.group(1))
                tokens = stripped.split()
                # tokens[0] = position, tokens[1..20] = match emissions (nats, * = -inf)
                raw_emissions = tokens[1:21]
                if len(raw_emissions) == 20:
                    probs = _emissions_to_probs(raw_emissions)
                    aa_probs = [
                        {"aa": aa, "prob": p}
                        for aa, p in zip(_HMMER_AA_ORDER, probs)
                    ]
                    aa_probs.sort(key=lambda x: x["prob"], reverse=True)
                    if pos % sample_every == 0:
                        positions.append({
                            "pos": pos,
                            "top_aa": aa_probs[:5],
                        })
                # Skip insert emissions (line 2) and transitions (line 3)
                line_idx += 3
                continue
        line_idx += 1

    return positions


def self_search_recovery(
    hmm_path: Path,
    seed_faa: Path,
    strict_bits: float = 45.0,
) -> dict:
    """Run hmmsearch of the HMM against its own seed sequences.

    This measures how well the HMM recovers the sequences used to build it.

    Parameters
    ----------
    hmm_path : Path
        The HMM profile.
    seed_faa : Path
        Seed protein FASTA.
    strict_bits : float
        Bit score threshold for counting a sequence as "recovered".

    Returns
    -------
    dict
        {recovered, total, recovery_rate, min_score, max_score}
        or default zeros on failure.
    """
    hmm_path = Path(hmm_path)
    seed_faa = Path(seed_faa)

    empty = {
        "recovered": 0,
        "total": 0,
        "recovery_rate": 0.0,
        "min_score": 0.0,
        "max_score": 0.0,
    }

    if not hmm_path.exists():
        print(f"ERROR: HMM not found: {hmm_path}", file=sys.stderr)
        return empty
    if not seed_faa.exists():
        print(f"ERROR: Seed FASTA not found: {seed_faa}", file=sys.stderr)
        return empty

    import tempfile, os
    with tempfile.NamedTemporaryFile(suffix=".tbl", delete=False) as tmp:
        tbl_path = Path(tmp.name)

    hmmsearch_bin = find_tool("hmmsearch") or "hmmsearch"
    cmd = [
        hmmsearch_bin,
        "--tblout", str(tbl_path),
        "--noali",
        str(hmm_path),
        str(seed_faa),
    ]

    result = run_cmd(cmd)
    if result.returncode != 0:
        print(f"ERROR: hmmsearch failed:\n{result.stderr}", file=sys.stderr)
        tbl_path.unlink(missing_ok=True)
        return empty

    # Parse tblout
    scores: list = []
    try:
        with tbl_path.open() as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    bit = float(parts[5])
                    scores.append(bit)
                except ValueError:
                    continue
    finally:
        tbl_path.unlink(missing_ok=True)

    # Count seed sequences
    total = 0
    try:
        from Bio import SeqIO
        total = sum(1 for _ in SeqIO.parse(str(seed_faa), "fasta"))
    except Exception:
        pass

    recovered = sum(1 for s in scores if s >= strict_bits)

    return {
        "recovered": recovered,
        "total": total,
        "recovery_rate": round(recovered / total, 4) if total > 0 else 0.0,
        "min_score": round(min(scores), 2) if scores else 0.0,
        "max_score": round(max(scores), 2) if scores else 0.0,
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_hmmbuild_stdout(stdout: str) -> dict:
    """Extract key stats from hmmbuild stdout text."""
    result: dict = {}

    version_match = re.search(r"HMMER\s+([\d.]+)", stdout)
    if version_match:
        result["version"] = version_match.group(1)

    for label, key in [
        ("LENG", "leng"),
        ("ALPH", "alph"),
        ("NSEQ", "nseq"),
        ("CKSUM", "cksum"),
    ]:
        m = re.search(rf"{label}\s+(\S+)", stdout)
        if m:
            val = m.group(1)
            result[key] = int(val) if label in ("LENG", "NSEQ") else val

    return result


def _emissions_to_probs(raw: list) -> list:
    """Convert HMMER log-odds emission strings to probabilities.

    HMMER stores emissions as -ln(p), where ``*`` represents -infinity (prob 0).
    Convert back to linear probabilities and normalise.
    """
    import math
    neg_log_probs = []
    for tok in raw:
        if tok == "*":
            neg_log_probs.append(float("inf"))
        else:
            try:
                neg_log_probs.append(float(tok))
            except ValueError:
                neg_log_probs.append(float("inf"))

    # Convert -ln(p) → p
    probs = []
    for val in neg_log_probs:
        if val == float("inf"):
            probs.append(0.0)
        else:
            probs.append(math.exp(-val))

    # Normalise to sum to 1
    total = sum(probs)
    if total > 0:
        probs = [p / total for p in probs]

    return probs
