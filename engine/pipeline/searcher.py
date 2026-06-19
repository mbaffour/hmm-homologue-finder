"""
pipeline/searcher.py — hmmsearch wrapper; handles protein and nucleotide DBs.

Runs hmmsearch against protein or nucleotide databases (translating the latter
via 6-frame translation first) and parses tblout/domtblout into DataFrames.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

from .utils import find_tool, run_cmd


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_hmmsearch_protein(
    hmm_path: Path,
    db_faa: Path,
    out_dir: Path,
    db_name: str,
    evalue: float = 1e-5,
    cpu: int = 0,
) -> dict:
    """Run hmmsearch against a protein database.

    Parameters
    ----------
    hmm_path : Path
        Profile HMM.
    db_faa : Path
        Target protein FASTA database.
    out_dir : Path
        Output directory for result files.
    db_name : str
        Short name for this database (used in output file names).
    evalue : float
        E-value inclusion threshold.
    cpu : int
        Number of threads. 0 = let HMMER decide.

    Returns
    -------
    dict
        {tblout (Path), domtblout (Path), out (Path), hit_count (int),
         strict_count (int, bit >= 45)}
    """
    hmm_path = Path(hmm_path)
    db_faa = Path(db_faa)
    out_dir = Path(out_dir)

    empty = {"tblout": Path(), "domtblout": Path(), "out": Path(),
             "hit_count": 0, "strict_count": 0}

    if not hmm_path.exists():
        print(f"ERROR: HMM not found: {hmm_path}", file=sys.stderr)
        return empty
    if not db_faa.exists():
        print(f"ERROR: Database not found: {db_faa}", file=sys.stderr)
        return empty

    out_dir.mkdir(parents=True, exist_ok=True)
    safe_name = _safe_name(db_name)
    tblout    = out_dir / f"{safe_name}.tblout"
    domtblout = out_dir / f"{safe_name}.domtblout"
    out_file  = out_dir / f"{safe_name}.out"

    hmmsearch_bin = find_tool("hmmsearch") or "hmmsearch"
    cmd = [
        hmmsearch_bin,
        "--tblout",    str(tblout),
        "--domtblout", str(domtblout),
        "-E",          str(evalue),
    ]
    if cpu > 0:
        cmd.extend(["--cpu", str(cpu)])
    cmd.extend([str(hmm_path), str(db_faa)])

    result = run_cmd(cmd)
    if result.returncode != 0:
        print(f"ERROR: hmmsearch failed:\n{result.stderr}", file=sys.stderr)
        return empty

    out_file.write_text(result.stdout)

    hits_df = parse_tblout(tblout)
    hit_count = len(hits_df)
    strict_count = int((hits_df["bit_score"] >= 45.0).sum()) if hit_count > 0 else 0

    return {
        "tblout": tblout,
        "domtblout": domtblout,
        "out": out_file,
        "hit_count": hit_count,
        "strict_count": strict_count,
    }


def run_hmmsearch_nucleotide(
    hmm_path: Path,
    db_fna: Path,
    out_dir: Path,
    db_name: str,
    scripts_dir: Path,
    evalue: float = 1e-5,
    min_aa: int = 30,
    cpu: int = 0,
) -> dict:
    """Translate a nucleotide DB via 6-frame then run hmmsearch.

    Parameters
    ----------
    hmm_path : Path
        Profile HMM.
    db_fna : Path
        Target nucleotide FASTA database.
    out_dir : Path
        Output directory.
    db_name : str
        Short database name.
    scripts_dir : Path
        Directory containing ``04_translate_sixframe.py``.
    evalue : float
        E-value threshold.
    min_aa : int
        Minimum ORF length passed to the translation script.
    cpu : int
        Number of HMMER threads.

    Returns
    -------
    dict
        Same structure as :func:`run_hmmsearch_protein`.
    """
    from .orf_prediction import predict_orfs_sixframe

    hmm_path   = Path(hmm_path)
    db_fna     = Path(db_fna)
    out_dir    = Path(out_dir)
    scripts_dir = Path(scripts_dir)

    empty = {"tblout": Path(), "domtblout": Path(), "out": Path(),
             "hit_count": 0, "strict_count": 0}

    if not db_fna.exists():
        print(f"ERROR: Nucleotide DB not found: {db_fna}", file=sys.stderr)
        return empty

    trans_dir = out_dir / "translated"
    faa_path = predict_orfs_sixframe(db_fna, trans_dir, scripts_dir, min_aa=min_aa)
    if not faa_path or not faa_path.exists():
        print("ERROR: Translation step failed.", file=sys.stderr)
        return empty

    return run_hmmsearch_protein(
        hmm_path=hmm_path,
        db_faa=faa_path,
        out_dir=out_dir,
        db_name=db_name,
        evalue=evalue,
        cpu=cpu,
    )


def parse_tblout(tblout_path: Path) -> pd.DataFrame:
    """Parse an hmmsearch ``--tblout`` file into a DataFrame.

    Columns returned: target_name, query_name, evalue, bit_score, bias_score,
    description.

    Parameters
    ----------
    tblout_path : Path
        Path to the tblout file.

    Returns
    -------
    pd.DataFrame
        Empty DataFrame (with correct columns) if file missing or no hits.
    """
    cols = ["target_name", "query_name", "evalue", "bit_score", "bias_score",
            "description"]

    tblout_path = Path(tblout_path)
    if not tblout_path.exists():
        return pd.DataFrame(columns=cols)

    rows: list = []
    try:
        with tblout_path.open() as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                # tblout format (space-delimited, description may contain spaces):
                # target_name accession query_name accession evalue bitscore bias
                # ... (18 fixed fields) ... description
                parts = line.split()
                if len(parts) < 18:
                    continue
                target_name = parts[0]
                query_name  = parts[2]
                evalue      = _safe_float(parts[4])
                bit_score   = _safe_float(parts[5])
                bias_score  = _safe_float(parts[6])
                description = " ".join(parts[18:]) if len(parts) > 18 else ""
                rows.append({
                    "target_name": target_name,
                    "query_name":  query_name,
                    "evalue":      evalue,
                    "bit_score":   bit_score,
                    "bias_score":  bias_score,
                    "description": description,
                })
    except Exception as exc:
        print(f"ERROR: Cannot parse tblout {tblout_path}: {exc}", file=sys.stderr)
        return pd.DataFrame(columns=cols)

    if not rows:
        return pd.DataFrame(columns=cols)

    return pd.DataFrame(rows, columns=cols)


def parse_domtblout(domtblout_path: Path) -> pd.DataFrame:
    """Parse an hmmsearch ``--domtblout`` file into a DataFrame.

    Columns: target_name, query_name, evalue, bit_score, bias_score,
    domain_evalue, domain_bit_score, hmm_from, hmm_to, ali_from, ali_to,
    env_from, env_to, description.

    Parameters
    ----------
    domtblout_path : Path

    Returns
    -------
    pd.DataFrame
        Empty DataFrame with correct columns on failure.
    """
    cols = [
        "target_name", "query_name", "evalue", "bit_score", "bias_score",
        "domain_evalue", "domain_bit_score",
        "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
        "description",
    ]

    domtblout_path = Path(domtblout_path)
    if not domtblout_path.exists():
        return pd.DataFrame(columns=cols)

    rows: list = []
    try:
        with domtblout_path.open() as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                # domtblout has 22 fixed fields then description
                if len(parts) < 22:
                    continue
                try:
                    row = {
                        "target_name":      parts[0],
                        "query_name":       parts[3],
                        "evalue":           _safe_float(parts[6]),
                        "bit_score":        _safe_float(parts[7]),
                        "bias_score":       _safe_float(parts[8]),
                        "domain_evalue":    _safe_float(parts[11]),
                        "domain_bit_score": _safe_float(parts[12]),
                        "hmm_from":         _safe_int(parts[15]),
                        "hmm_to":           _safe_int(parts[16]),
                        "ali_from":         _safe_int(parts[17]),
                        "ali_to":           _safe_int(parts[18]),
                        "env_from":         _safe_int(parts[19]),
                        "env_to":           _safe_int(parts[20]),
                        "description":      " ".join(parts[22:]) if len(parts) > 22 else "",
                    }
                    rows.append(row)
                except (IndexError, ValueError) as exc:
                    print(f"WARNING: Skipping malformed domtblout line: {exc}", file=sys.stderr)
                    continue
    except Exception as exc:
        print(f"ERROR: Cannot parse domtblout {domtblout_path}: {exc}", file=sys.stderr)
        return pd.DataFrame(columns=cols)

    if not rows:
        return pd.DataFrame(columns=cols)

    return pd.DataFrame(rows, columns=cols)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val: str) -> float:
    try:
        return float(val)
    except (ValueError, TypeError):
        return float("nan")


def _safe_int(val: str) -> int:
    try:
        return int(val)
    except (ValueError, TypeError):
        return 0


def _safe_name(name: str) -> str:
    """Convert a database name to a safe filename stem."""
    import re
    return re.sub(r"[^\w\-]", "_", name)[:80]
