"""
pipeline/utils.py — Shared utilities for pipeline modules.

Key exports:
  find_tool(name)   → full path or None, using augmented PATH
  get_env()         → os.environ copy with augmented PATH for subprocess
  run_cmd(cmd)      → subprocess.CompletedProcess with augmented env
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional


def _extra_paths() -> list[str]:
    """Return extra bin directories beyond the current PATH."""
    extras: list[str] = []
    # Active conda env
    for var in ("CONDA_PREFIX", "CONDA_DEFAULT_ENV"):
        val = os.environ.get(var)
        if val:
            candidate = os.path.join(val, "bin")
            if os.path.isdir(candidate):
                extras.append(candidate)
                break
    # Common miniforge / miniconda / anaconda locations
    home = Path.home()
    user_scripts = home / "Library" / "Python" / f"{sys.version_info.major}.{sys.version_info.minor}" / "bin"
    if user_scripts.exists():
        extras.append(str(user_scripts))
    local_bin = home / ".local" / "bin"
    if local_bin.exists():
        extras.append(str(local_bin))
    for base in [
        home / "miniforge3",
        home / "miniconda3",
        home / "anaconda3",
        Path("/opt/anaconda3"),
        Path("/opt/miniconda3"),
        Path("/usr/local/miniconda3"),
    ]:
        # Base conda bin
        d = base / "bin"
        if d.exists():
            extras.append(str(d))
        # Every named environment — catches meme-tools, hmm_env, etc.
        envs_dir = base / "envs"
        if envs_dir.is_dir():
            for env_bin in sorted(envs_dir.glob("*/bin")):
                if env_bin.is_dir():
                    extras.append(str(env_bin))
    return extras


def _augmented_path() -> str:
    """Return PATH string augmented with conda env bin dirs."""
    return os.pathsep.join([os.environ.get("PATH", "")] + _extra_paths())


def ensure_tools_on_path() -> None:
    """Mutate ``os.environ['PATH']`` in place so conda-env tools are findable
    by *in-process* libraries (e.g. toyplot looking for ``gs``, pygenomeviz,
    matplotlib backends). Safe to call multiple times — only appends dirs
    that are not already present.

    Call once at application startup so every in-process subprocess and
    library can locate the bioinformatics binaries without each call site
    having to pass an augmented env.
    """
    current = os.environ.get("PATH", "")
    current_dirs = current.split(os.pathsep)
    additions = [d for d in _extra_paths() if d and d not in current_dirs]
    if additions:
        os.environ["PATH"] = os.pathsep.join(current_dirs + additions)


def find_tool(name: str) -> Optional[str]:
    """Return the full path to *name* (searching augmented PATH) or None."""
    return shutil.which(name, path=_augmented_path())


def get_env() -> dict:
    """Return a copy of os.environ with PATH augmented for conda envs."""
    env = dict(os.environ)
    env["PATH"] = _augmented_path()
    return env


def run_cmd(
    cmd: list,
    cwd: Optional[Path] = None,
    timeout: int = 3600,
    text: bool = True,
    capture_output: bool = True,
) -> subprocess.CompletedProcess:
    """Run *cmd* using the augmented environment.

    Automatically resolves the executable to a full path so that subprocess
    can find tools installed in conda environments that are not on the
    system PATH.
    """
    env = get_env()
    # Resolve executable to full path when possible
    if cmd and not os.path.isabs(cmd[0]):
        resolved = find_tool(cmd[0])
        if resolved:
            cmd = [resolved] + list(cmd[1:])
    return subprocess.run(
        cmd,
        env=env,
        cwd=str(cwd) if cwd else None,
        capture_output=capture_output,
        text=text,
        timeout=timeout,
    )
