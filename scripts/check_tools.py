#!/usr/bin/env python3
"""
check_tools.py — preflight: are all required tools installed?
=============================================================
Verifies every external program and Python package the pipeline needs is
available before a (long) run starts. Importable by hmm_finder.py, and also
runnable on its own:

    python3 check_tools.py            # report only
    python3 check_tools.py --install  # try to auto-install anything missing

Auto-install runs the tool folder's setup.sh, which creates/updates the
`hmm-discovery` conda environment from the deployable repo's environment.yml.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# External command-line tools the workflow calls.
REQUIRED_TOOLS = [
    "hmmsearch", "hmmbuild", "hmmscan",   # HMMER
    "mafft", "trimal",                    # alignment
    "prodigal", "seqkit",                 # ORF calling / FASTA handling
    "cd-hit",                             # clustering
    "iqtree",                             # phylogeny
    "meme", "fimo",                       # motifs
    "clinker",                            # synteny figures
    "curl",                               # downloads
]
# Python packages.
REQUIRED_PYTHON = ["Bio", "pandas"]


def check() -> tuple[list[str], list[str]]:
    """Return (missing_tools, missing_python)."""
    missing_tools = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    missing_python = []
    for mod in REQUIRED_PYTHON:
        try:
            __import__(mod)
        except Exception:
            missing_python.append(mod)
    return missing_tools, missing_python


def report(missing_tools: list[str], missing_python: list[str]) -> None:
    print("Tool check:")
    for t in REQUIRED_TOOLS:
        ok = shutil.which(t) is not None
        print(f"  [{'OK ' if ok else 'MISSING'}] {t}")
    for m in REQUIRED_PYTHON:
        try:
            __import__(m)
            ok = True
        except Exception:
            ok = False
        print(f"  [{'OK ' if ok else 'MISSING'}] python:{m}")
    if missing_tools or missing_python:
        print(f"\n  Missing tools:  {missing_tools or 'none'}")
        print(f"  Missing python: {missing_python or 'none'}")
    else:
        print("\n  All required software is installed.")


def ensure(install: bool = False) -> bool:
    """Check; optionally auto-install. Returns True if everything is present."""
    missing_tools, missing_python = check()
    report(missing_tools, missing_python)
    if not missing_tools and not missing_python:
        return True
    if not install:
        print("\nRun with --install to set up the environment automatically,")
        print("or run  ./setup.sh  in the tool folder.")
        return False
    setup = Path(__file__).resolve().parent.parent / "setup.sh"
    if not setup.exists():
        print(f"setup.sh not found at {setup}")
        return False
    print("\nInstalling missing software via setup.sh …")
    subprocess.run(["bash", str(setup)], check=False)
    missing_tools, missing_python = check()
    report(missing_tools, missing_python)
    return not (missing_tools or missing_python)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--install", action="store_true", help="auto-install anything missing")
    args = ap.parse_args()
    ok = ensure(install=args.install)
    sys.exit(0 if ok else 1)
