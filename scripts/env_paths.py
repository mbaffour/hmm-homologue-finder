#!/usr/bin/env python3
"""
env_paths.py — put the `hmm-discovery` conda env's bin/ on PATH.

Single source of truth shared by the standalone helper scripts so the bundled
tools (hmmsearch, prodigal, mafft, cd-hit, iqtree, clinker, …) resolve even when
a script is launched by a process that never ran `conda activate`.
"""
from __future__ import annotations

import os
from pathlib import Path


def ensure_env_on_path(env_name: str = "hmm-discovery"):
    """Prepend the conda env's bin dir to PATH; return the dir used or None.

    Search order: the active CONDA_PREFIX first (if set), then the env under the
    common Miniforge/Miniconda install roots in the user's home directory.
    """
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(Path(os.environ["CONDA_PREFIX"]) / "bin")
    candidates += [Path.home() / _n / "envs" / env_name / "bin"
                   for _n in ("miniforge3", "mambaforge", "miniconda3", "anaconda3")]
    for b in candidates:
        if b.is_dir():
            os.environ["PATH"] = f"{b}{os.pathsep}{os.environ.get('PATH', '')}"
            return str(b)
    return None
