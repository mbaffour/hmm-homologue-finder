"""
core/logger.py — Structured audit trail + tool availability checker.

Every subprocess call is logged to logs/audit_trail.jsonl with:
  timestamp, step, command, tool_version, input_hashes, output_hashes,
  exit_code, duration_sec
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Optional tools the app can use; keyed by the executable name
# ---------------------------------------------------------------------------
OPTIONAL_TOOLS: dict[str, str] = {
    "diamond":    "DIAMOND reciprocal BLAST",
    "cd-hit":     "Sequence clustering (CD-HIT)",
    "cd-hit-est": "NT clustering (CD-HIT-EST)",
    "mmseqs":     "Sequence clustering (MMseqs2)",
    "meme":       "Motif discovery (MEME)",
    "fimo":       "Motif scanning (FIMO)",
    "phobius.pl": "TM/signal peptide (Phobius)",
    "tmhmm":      "TM topology (TMHMM)",
    "foldseek":   "Structural similarity (Foldseek)",
    "iqtree2":    "Phylogenetics (IQ-TREE 2)",
    "iqtree":     "Phylogenetics (IQ-TREE)",
    "prodigal":   "ORF prediction (Prodigal)",
    "mafft":      "Multiple alignment (MAFFT)",
    "clustalo":   "Multiple alignment (Clustal Omega)",
    "trimal":     "Alignment trimming (trimAl)",
    "hmmbuild":   "HMM building (HMMER)",
    "hmmsearch":  "HMM search (HMMER)",
}

REQUIRED_TOOLS: list[str] = ["hmmbuild", "hmmsearch", "mafft", "trimal"]


def _extra_search_paths() -> list[str]:
    """Return extra directories to search for tools beyond the current PATH."""
    extras = []
    # Active conda environment (CONDA_PREFIX or CONDA_DEFAULT_ENV)
    conda_prefix = os.environ.get("CONDA_PREFIX") or os.environ.get("CONDA_DEFAULT_ENV")
    if conda_prefix and os.path.isdir(os.path.join(conda_prefix, "bin")):
        extras.append(os.path.join(conda_prefix, "bin"))
    # Common miniforge/miniconda/anaconda locations
    home = Path.home()
    user_scripts = home / "Library" / "Python" / f"{sys.version_info.major}.{sys.version_info.minor}" / "bin"
    if user_scripts.exists():
        extras.append(str(user_scripts))
    local_bin = home / ".local" / "bin"
    if local_bin.exists():
        extras.append(str(local_bin))
    for base in [home / "miniforge3", home / "miniconda3", home / "anaconda3",
                 Path("/opt/anaconda3"), Path("/opt/miniconda3")]:
        for sub in ["bin", "envs/hmm_env/bin"]:
            d = base / sub
            if d.exists():
                extras.append(str(d))
    return extras


def check_tools(proj_dir: Path) -> dict[str, dict]:
    """
    Discover which tools are on PATH (plus common conda env locations).
    Writes reports/tools.json.
    Returns a dict keyed by tool name with keys: available, path, version.
    """
    result: dict[str, dict] = {}
    all_tools = {**{t: OPTIONAL_TOOLS[t] for t in OPTIONAL_TOOLS}}
    for t in REQUIRED_TOOLS:
        all_tools.setdefault(t, t)

    # Build an augmented PATH that includes conda env bin dirs
    extra_paths = _extra_search_paths()
    augmented_path = os.pathsep.join(
        [os.environ.get("PATH", "")] + extra_paths
    )

    for exe, desc in all_tools.items():
        path = shutil.which(exe, path=augmented_path)
        version = _get_version(exe, path) if path else None
        result[exe] = {
            "available": path is not None,
            "path": path,
            "version": version,
            "description": desc,
            "required": exe in REQUIRED_TOOLS,
        }

    reports_dir = proj_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    out = reports_dir / "tools.json"
    out.write_text(json.dumps(result, indent=2))
    return result


def _get_version(exe: str, path: Optional[str]) -> Optional[str]:
    if not path:
        return None
    flags = {
        "hmmbuild":  ["-h"],
        "hmmsearch": ["-h"],
        "mafft":     ["--version"],
        "trimal":    ["--version"],
        "iqtree2":   ["--version"],
        "iqtree":    ["--version"],
        "iqtree3":   ["--version"],
        "prodigal":  ["-v"],
        "clustalo":  ["--version"],
        "diamond":   ["version"],
        "meme":      ["--version"],
        "fimo":      ["--version"],
        "cd-hit":    ["-h"],
        "cd-hit-est":["-h"],
        "mmseqs":    ["version"],
    }
    # Patterns to extract the version string from tool output
    import re
    version_patterns = {
        "hmmbuild":  r"HMMER\s+(\S+)",
        "hmmsearch": r"HMMER\s+(\S+)",
        "cd-hit":    r"CD-HIT version\s+(\S+)",
        "cd-hit-est":r"CD-HIT version\s+(\S+)",
    }
    try:
        flag = flags.get(exe, ["--version"])
        proc = subprocess.run(
            [path] + flag,
            capture_output=True, text=True, timeout=10,
        )
        combined = (proc.stdout + proc.stderr).strip()
        pat = version_patterns.get(exe)
        if pat:
            m = re.search(pat, combined)
            if m:
                return m.group(1)
        lines = combined.splitlines()
        return lines[0] if lines else "unknown"
    except Exception:
        return "unknown"


# ---------------------------------------------------------------------------
# Audit trail
# ---------------------------------------------------------------------------

class AuditLogger:
    """
    Append-only structured logger. Each call to record() writes one JSON
    line to logs/audit_trail.jsonl.
    """

    def __init__(self, proj_dir: Path):
        self.proj_dir = proj_dir
        logs_dir = proj_dir / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        self.trail_file = logs_dir / "audit_trail.jsonl"
        self.run_log = logs_dir / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    def record(
        self,
        step: str,
        command: list[str],
        exit_code: int,
        duration_sec: float,
        input_files: Optional[dict[str, str]] = None,
        output_files: Optional[dict[str, str]] = None,
        extra: Optional[dict] = None,
    ) -> None:
        entry = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "step": step,
            "command": " ".join(str(c) for c in command),
            "exit_code": exit_code,
            "duration_sec": round(duration_sec, 2),
            "input_files": input_files or {},
            "output_files": output_files or {},
        }
        if extra:
            entry.update(extra)
        with open(self.trail_file, "a") as fh:
            fh.write(json.dumps(entry) + "\n")

    def log(self, message: str) -> None:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {message}"
        with open(self.run_log, "a") as fh:
            fh.write(line + "\n")

    def read_trail(self) -> list[dict]:
        if not self.trail_file.exists():
            return []
        records = []
        for line in self.trail_file.read_text().splitlines():
            line = line.strip()
            if line:
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
        return records


def file_hash(path: Path) -> str:
    """SHA-256 of file contents, truncated to 12 hex chars."""
    if not path.exists():
        return "missing"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()[:12]
