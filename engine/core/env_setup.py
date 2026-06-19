"""
core/env_setup.py — First-run environment checker and guided installer.

Called by the app on startup (or by the Setup tab) to:
  1. Check which required / optional tools are available
  2. Check Python package versions
  3. Provide install instructions and a runnable setup command
  4. Return a structured report consumable by the UI

No imports outside the standard library + the project's own utils.
"""
from __future__ import annotations

import importlib
import os
import platform
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional


# ---- Tool and package requirements -----------------------------------------

REQUIRED_TOOLS: list[dict] = [
    {"name": "hmmbuild",  "pkg": "hmmer",  "channel": "bioconda",
     "desc": "Build profile HMMs from MSA", "auto_install": True},
    {"name": "hmmsearch", "pkg": "hmmer",  "channel": "bioconda",
     "desc": "Search protein / nucleotide databases", "auto_install": True},
    {"name": "mafft",     "pkg": "mafft",  "channel": "bioconda",
     "desc": "Multiple sequence alignment", "auto_install": True},
    {"name": "trimal",    "pkg": "trimal", "channel": "bioconda",
     "desc": "Alignment column trimming", "auto_install": True},
]

OPTIONAL_TOOLS: list[dict] = [
    {"name": "iqtree2",   "pkg": "iqtree2",  "channel": "bioconda",
     "desc": "Phylogenetic tree inference", "full_run": True, "auto_install": True},
    {"name": "iqtree",    "pkg": "iqtree",   "channel": "bioconda",
     "desc": "Phylogenetic tree inference (older)", "full_run": True, "auto_install": True},
    {"name": "clustalo",  "pkg": "clustalo", "channel": "bioconda",
     "desc": "Alternative multiple sequence aligner"},
    {"name": "prodigal",  "pkg": "prodigal", "channel": "bioconda",
     "desc": "Gene prediction for nucleotide input", "full_run": True, "auto_install": True},
    {"name": "cd-hit",    "pkg": "cd-hit",   "channel": "bioconda",
     "desc": "Sequence clustering (fast)", "full_run": True, "auto_install": True},
    {"name": "mmseqs",    "pkg": "mmseqs2",  "channel": "bioconda",
     "desc": "Sequence clustering (ultra-fast)", "full_run": True, "auto_install": True},
    {"name": "diamond",   "pkg": "diamond",  "channel": "bioconda",
     "desc": "Reciprocal BLAST confirmation", "auto_install": True},
    {"name": "meme",      "pkg": "meme",     "channel": "bioconda",
     "desc": "Motif discovery (MEME suite)", "full_run": True, "auto_install": True},
    {"name": "fimo",      "pkg": "meme",     "channel": "bioconda",
     "desc": "Motif scanning (part of MEME suite)", "full_run": True, "auto_install": True},
    {"name": "clinker",   "pkg": "clinker", "channel": "pip",
     "desc": "Interactive synteny gene-cluster comparison", "full_run": True, "auto_install": True},
    {"name": "gs",        "pkg": "ghostscript",   "channel": "conda-forge",
     "desc": "PNG export for phylogenetic tree figures", "full_run": True, "auto_install": True},
    {"name": "phobius.pl","pkg": "phobius",  "channel": "bioconda",
     "desc": "TM topology + signal peptide", "auto_install": False},
    {"name": "tmhmm",     "pkg": "tmhmm",    "channel": "bioconda",
     "desc": "TM topology prediction", "auto_install": False},
    {"name": "foldseek",  "pkg": "foldseek", "channel": "bioconda",
     "desc": "Structural similarity search", "auto_install": True},
]

FULL_RUN_TOOL_GROUPS: list[dict] = [
    {
        "name": "iqtree",
        "desc": "Phylogenetic tree inference",
        "any_of": ["iqtree2", "iqtree"],
    },
]

REQUIRED_PYTHON: list[dict] = [
    {"pkg": "shiny",        "import": "shiny",       "min_version": "1.5.0"},
    {"pkg": "shinyswatch",  "import": "shinyswatch",  "min_version": "0.10.0"},
    {"pkg": "pandas",       "import": "pandas",       "min_version": "2.0.0"},
    {"pkg": "plotly",       "import": "plotly",       "min_version": "5.0.0"},
    {"pkg": "biopython",    "import": "Bio",          "min_version": "1.81"},
    {"pkg": "jinja2",       "import": "jinja2",       "min_version": "3.0.0"},
    {"pkg": "matplotlib",   "import": "matplotlib",   "min_version": "3.7.0"},
    {"pkg": "numpy",        "import": "numpy",        "min_version": "1.24.0"},
    {"pkg": "aiohttp",      "import": "aiohttp",      "min_version": "3.9.0"},
    {"pkg": "aiofiles",     "import": "aiofiles",     "min_version": "23.0.0"},
    {"pkg": "openpyxl",     "import": "openpyxl",     "min_version": "3.1.0"},
    {"pkg": "urllib3<2",    "import": "urllib3",      "min_version": "1.26.0",
     "max_version": "1.99.99"},
    {"pkg": "scipy",        "import": "scipy",        "min_version": "1.10.0",
     "full_run": True},
    {"pkg": "toytree<3.0.11", "import": "toytree",    "min_version": "3.0.0",
     "max_version": "3.0.10",
     "full_run": True},
    {"pkg": "toyplot",      "import": "toyplot",      "min_version": "1.0.0",
     "full_run": True},
    {"pkg": "pygenomeviz",  "import": "pygenomeviz",  "min_version": "0.4.0",
     "full_run": True},
    {"pkg": "kaleido",      "import": "kaleido",      "min_version": "0.2.1",
     "full_run": True},
]


# ---- Augmented PATH (mirrors core/logger.py) --------------------------------

def _augmented_path() -> str:
    extras: list[str] = []
    for var in ("CONDA_PREFIX", "CONDA_DEFAULT_ENV"):
        val = os.environ.get(var)
        if val:
            c = os.path.join(val, "bin")
            if os.path.isdir(c):
                extras.append(c)
                break
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
    return os.pathsep.join([os.environ.get("PATH", "")] + extras)


def _find_tool(name: str) -> Optional[str]:
    return shutil.which(name, path=_augmented_path())


# ---- Conda detection --------------------------------------------------------

def _find_conda() -> Optional[str]:
    """Return path to mamba or conda, preferring mamba."""
    aug = _augmented_path()
    for cmd in ("mamba", "conda"):
        p = shutil.which(cmd, path=aug)
        if p:
            return p
    return None


def _tool_by_name(name: str) -> dict:
    for tool in REQUIRED_TOOLS + OPTIONAL_TOOLS:
        if tool["name"] == name:
            return tool
    return {"name": name, "pkg": name, "channel": "conda-forge"}


def _missing_full_run_tools(report: dict) -> list[dict]:
    """Return missing, auto-installable tools for a complete run."""
    available = {
        t["name"]: t["available"]
        for t in report.get("required_tools", []) + report.get("optional_tools", [])
    }
    missing: list[dict] = []

    for tool in report.get("required_tools", []):
        if not tool.get("available"):
            missing.append(tool)

    for tool in report.get("optional_tools", []):
        if tool.get("full_run") and tool.get("auto_install", True) and not tool.get("available"):
            if any(tool["name"] in group["any_of"] for group in FULL_RUN_TOOL_GROUPS):
                continue
            missing.append(tool)

    for group in FULL_RUN_TOOL_GROUPS:
        if not any(available.get(name) for name in group["any_of"]):
            first = _tool_by_name(group["any_of"][0])
            missing.append({
                **first,
                "name": group["name"],
                "desc": group["desc"],
                "alternatives": group["any_of"],
            })

    seen: set[str] = set()
    unique: list[dict] = []
    for tool in missing:
        key = tool.get("name", tool.get("pkg", ""))
        if key and key not in seen:
            seen.add(key)
            unique.append(tool)
    return unique


def _missing_full_run_python(report: dict) -> list[dict]:
    return [
        p for p in report.get("python_packages", [])
        if p.get("full_run") and not p.get("ok")
    ]


# ---- Version parsing --------------------------------------------------------

def _parse_version(ver_str: str) -> tuple:
    """Convert '3.4.1' → (3, 4, 1) for comparison."""
    parts = []
    for p in str(ver_str).split(".")[:3]:
        try:
            parts.append(int(p))
        except ValueError:
            parts.append(0)
    return tuple(parts)


# ---- Public API -------------------------------------------------------------

def check_environment() -> dict:
    """
    Check the full environment and return a structured report.

    Returns
    -------
    dict with keys:
        python_version : str
        platform       : str
        conda_available: bool
        conda_path     : str | None
        required_tools : list[dict]  — {name, available, path, desc, required}
        optional_tools : list[dict]
        python_packages: list[dict]  — {pkg, installed, version, ok}
        all_required_ok: bool
        all_python_ok  : bool
        install_cmd    : str  — one-liner to fix everything
        setup_script   : str  — path to setup_environment.sh
    """
    report: dict = {
        "python_version": sys.version.split()[0],
        "platform":       platform.platform(),
        "conda_available": False,
        "conda_path":     None,
        "required_tools": [],
        "optional_tools": [],
        "python_packages": [],
        "all_required_ok": False,
        "all_python_ok":   False,
        "all_full_run_ok": False,
        "missing_full_run_tools": [],
        "missing_full_run_python": [],
        "install_cmd":    "",
        "setup_script":   str(Path(__file__).parent.parent / "setup_environment.sh"),
    }

    # Conda
    conda = _find_conda()
    report["conda_available"] = conda is not None
    report["conda_path"] = conda

    # Required tools
    all_req_ok = True
    for tool in REQUIRED_TOOLS:
        path = _find_tool(tool["name"])
        entry = {
            "name": tool["name"],
            "available": path is not None,
            "path": path,
            "desc": tool["desc"],
            "pkg": tool["pkg"],
            "channel": tool["channel"],
            "required": True,
            "auto_install": tool.get("auto_install", True),
            "full_run": True,
        }
        report["required_tools"].append(entry)
        if not path:
            all_req_ok = False

    # Optional tools
    for tool in OPTIONAL_TOOLS:
        path = _find_tool(tool["name"])
        report["optional_tools"].append({
            "name": tool["name"],
            "available": path is not None,
            "path": path,
            "desc": tool["desc"],
            "pkg": tool["pkg"],
            "channel": tool["channel"],
            "required": False,
            "auto_install": tool.get("auto_install", True),
            "full_run": tool.get("full_run", False),
        })

    report["all_required_ok"] = all_req_ok

    # Python packages
    all_py_ok = True
    for pkgdef in REQUIRED_PYTHON:
        try:
            mod = importlib.import_module(pkgdef["import"])
            ver = getattr(mod, "__version__", "unknown")
            min_ver = pkgdef["min_version"]
            max_ver = pkgdef.get("max_version")
            ok = True
            if ver != "unknown":
                ok = _parse_version(ver) >= _parse_version(min_ver)
                if max_ver:
                    ok = ok and _parse_version(ver) <= _parse_version(max_ver)
            report["python_packages"].append({
                "pkg": pkgdef["pkg"],
                "installed": True,
                "version": ver,
                "ok": ok,
                "min_version": min_ver,
                "max_version": max_ver,
                "full_run": pkgdef.get("full_run", False),
            })
            if not ok:
                all_py_ok = False
        except Exception as exc:
            installed = importlib.util.find_spec(pkgdef["import"]) is not None
            report["python_packages"].append({
                "pkg": pkgdef["pkg"],
                "installed": installed,
                "version": None,
                "ok": False,
                "min_version": pkgdef["min_version"],
                "max_version": pkgdef.get("max_version"),
                "full_run": pkgdef.get("full_run", False),
                "error": str(exc),
            })
            all_py_ok = False

    report["all_python_ok"] = all_py_ok

    missing_full_tools = _missing_full_run_tools(report)
    missing_full_py = _missing_full_run_python(report)
    report["missing_full_run_tools"] = missing_full_tools
    report["missing_full_run_python"] = missing_full_py
    report["all_full_run_ok"] = (
        report["all_required_ok"]
        and report["all_python_ok"]
        and not missing_full_tools
        and not missing_full_py
    )

    # Generate install command
    missing_tools = [t["pkg"] for t in missing_full_tools if t.get("channel") != "pip"]
    missing_pip_tools = [t["pkg"] for t in missing_full_tools if t.get("channel") == "pip"]
    missing_py = [p["pkg"] for p in report["python_packages"] if not p["ok"]]

    parts = []
    if missing_tools:
        parts.append(
            f"conda install -c bioconda -c conda-forge {' '.join(set(missing_tools))}"
        )
    if missing_pip_tools or missing_py:
        quoted = " ".join(shlex.quote(pkg) for pkg in sorted(set(missing_pip_tools + missing_py)))
        parts.append(f"pip install {quoted}")
    report["install_cmd"] = " && ".join(parts) if parts else "# All tools already installed ✓"

    return report


def install_missing(
    report: dict,
    env_name: str = "hmm_env",
    dry_run: bool = False,
    include_optional: bool = True,
) -> list[str]:
    """
    Install missing required (and optionally optional) tools and Python packages.

    Parameters
    ----------
    report : dict
        Output of :func:`check_environment`.
    env_name : str
        Conda environment to install into.
    dry_run : bool
        If True, return commands without executing them.
    include_optional : bool
        When True (default), also install missing optional tools.

    Returns
    -------
    list[str]
        Log lines describing actions taken (or to be taken in dry-run mode).
    """
    log: list[str] = []
    conda = report.get("conda_path")

    # Collect missing tools — required always, optional when requested
    all_tools = list(report.get("required_tools", []))
    if include_optional:
        all_tools += list(report.get("optional_tools", []))

    missing_conda = [
        t for t in all_tools
        if not t["available"] and t.get("auto_install", True) and t.get("channel") != "pip"
    ]
    missing_pip_tools = [
        t["pkg"] for t in all_tools
        if not t["available"] and t.get("auto_install", True) and t.get("channel") == "pip"
    ]

    if missing_conda and conda:
        seen_pkgs: set[str] = set()
        for tool in missing_conda:
            pkg = tool["pkg"]
            if pkg in seen_pkgs:
                continue
            seen_pkgs.add(pkg)
            channel = tool.get("channel", "bioconda")
            cmd = [conda, "install", "-n", env_name, "-y", "-c", channel, "-c", "conda-forge", pkg]
            log.append(f"Running: {' '.join(cmd)}")
            if not dry_run:
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
                    if result.returncode == 0:
                        log.append(f"Conda install OK: {pkg}")
                    else:
                        log.append(f"Conda install failed for {pkg}:\n{result.stderr[-800:]}")
                except Exception as exc:
                    log.append(f"Conda install error for {pkg}: {exc}")
    elif missing_conda:
        log.append(
            "Cannot install tools automatically — conda not found. "
            "Run setup_environment.sh to install."
        )

    # Install missing Python packages
    missing_py = [p["pkg"] for p in report.get("python_packages", []) if not p["ok"]]
    pip_pkgs = sorted(set(missing_pip_tools + missing_py))
    if pip_pkgs:
        cmd = [sys.executable, "-m", "pip", "install", "--quiet"] + pip_pkgs
        log.append(f"Running: {' '.join(cmd)}")
        if not dry_run:
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                if result.returncode == 0:
                    log.append("pip install: OK")
                else:
                    log.append(f"pip install failed:\n{result.stderr[-500:]}")
            except Exception as exc:
                log.append(f"pip install error: {exc}")

    if not log:
        log.append("All tools already installed — nothing to do.")

    return log


def first_run_needed(app_dir: Optional[Path] = None) -> bool:
    """Return True if this appears to need full-run environment setup."""
    report = check_environment()
    return not report.get("all_full_run_ok", False)


async def auto_install_async(
    log_callback,
    env_name: str = "hmm_env",
    include_optional: bool = False,
) -> bool:
    """
    Install missing tools asynchronously, streaming log lines via *log_callback*.

    Parameters
    ----------
    log_callback : callable
        Called with a single str log line as each line arrives.
    env_name : str
        Conda environment to install into (default: hmm_env).
    include_optional : bool
        Also install optional tools. False by default at startup to keep it fast.

    Returns
    -------
    bool
        True if all required tools are now available, False otherwise.
    """
    import asyncio

    report = check_environment()
    conda = report.get("conda_path")

    # Collect missing tools
    missing_req = [
        t for t in report.get("required_tools", [])
        if not t["available"] and t.get("auto_install", True)
    ]
    if include_optional:
        missing_req += [
            t for t in report.get("optional_tools", [])
            if not t["available"] and t.get("auto_install", True)
        ]
    else:
        missing_req += [
            t for t in report.get("missing_full_run_tools", [])
            if t.get("auto_install", True)
        ]
    missing_py = [p["pkg"] for p in report.get("python_packages", []) if not p["ok"]]

    if not missing_req and not missing_py:
        log_callback("✅ Full-run environment already installed — nothing to do.")
        return True

    # ── Conda tools ─────────────────────────────────────────────────────────
    missing_conda = [t for t in missing_req if t.get("channel") != "pip"]
    missing_pip_tools = [t["pkg"] for t in missing_req if t.get("channel") == "pip"]

    if missing_conda:
        if not conda:
            log_callback(
                "❌ conda / mamba not found — cannot auto-install.\n"
                "   Run:  bash setup_environment.sh"
            )
        else:
            seen_pkgs: set[str] = set()
            for t in missing_conda:
                pkg = t["pkg"]
                if pkg in seen_pkgs:
                    continue
                seen_pkgs.add(pkg)
                channel = t.get("channel", "bioconda")
                cmd = [conda, "install", "-n", env_name, "-y", "-c", channel, "-c", "conda-forge", pkg]

                log_callback(f"📦 Installing: {pkg}")
                log_callback(f"$ {' '.join(cmd)}")
                try:
                    proc = await asyncio.create_subprocess_exec(
                        *cmd,
                        stdout=asyncio.subprocess.PIPE,
                        stderr=asyncio.subprocess.STDOUT,
                    )
                    assert proc.stdout is not None
                    async for raw in proc.stdout:
                        line = raw.decode(errors="replace").rstrip()
                        if line:
                            log_callback(line)
                    await proc.wait()
                    if proc.returncode == 0:
                        log_callback(f"✅ conda install: {pkg} done")
                    else:
                        log_callback(f"❌ conda install {pkg} exited {proc.returncode}")
                except Exception as exc:
                    log_callback(f"❌ conda error installing {pkg}: {exc}")

    # ── Python packages ──────────────────────────────────────────────────────
    pip_pkgs = sorted(set(missing_py + missing_pip_tools))
    if pip_pkgs:
        cmd_py = [sys.executable, "-m", "pip", "install"] + pip_pkgs
        log_callback(f"📦 pip install: {' '.join(pip_pkgs)}")
        try:
            proc_py = await asyncio.create_subprocess_exec(
                *cmd_py,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT,
            )
            assert proc_py.stdout is not None
            async for raw in proc_py.stdout:
                line = raw.decode(errors="replace").rstrip()
                if line:
                    log_callback(line)
            await proc_py.wait()
            log_callback(
                "✅ pip install: done" if proc_py.returncode == 0
                else f"❌ pip install exited {proc_py.returncode}"
            )
        except Exception as exc:
            log_callback(f"❌ pip error: {exc}")

    # ── Re-check ─────────────────────────────────────────────────────────────
    report2 = check_environment()
    still_missing = [
        t["name"] for t in report2.get("required_tools", []) if not t["available"]
    ] + [
        p["pkg"] for p in report2.get("python_packages", []) if not p["ok"]
    ]
    if still_missing:
        log_callback(f"⚠️  Still missing: {', '.join(still_missing)}")
        return False

    if report2.get("all_full_run_ok"):
        log_callback("🎉 Full-run environment installed and ready!")
    else:
        missing_full = [t["name"] for t in report2.get("missing_full_run_tools", [])]
        missing_full += [p["pkg"] for p in report2.get("missing_full_run_python", [])]
        if missing_full:
            log_callback(f"⚠️ Full-run optional capabilities still missing: {', '.join(missing_full)}")
    log_callback("🎉 All required tools installed and ready!")
    return True
