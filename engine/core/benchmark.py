"""
core/benchmark.py - App helpers for launching and monitoring benchmarks.

The exhaustive benchmark can run for hours or days, so the app launches it as
an independent process and reads progress from the runner's manifest/log files.
"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
from pathlib import Path
from typing import Any


DEFAULT_BENCHMARK_ROOT = Path.home() / "Documents" / "HMM-Discovery-Benchmark"


def resolve_path(value: str | Path, app_dir: Path | None = None) -> Path:
    """Resolve a user-entered path, expanding ~ and app-relative paths."""
    raw = Path(str(value).strip()).expanduser()
    if raw.is_absolute() or app_dir is None:
        return raw.resolve()
    return (app_dir / raw).resolve()


def runner_script(app_dir: Path) -> Path:
    path = app_dir / "scripts" / "run_all_database_benchmark.py"
    if not path.exists():
        raise FileNotFoundError(f"Benchmark runner not found: {path}")
    return path


def manifest_path(out_dir: str | Path) -> Path:
    return Path(out_dir).expanduser().resolve() / "benchmark_manifest.json"


def pid_path(out_dir: str | Path) -> Path:
    return Path(out_dir).expanduser().resolve() / "benchmark.pid"


def log_path(out_dir: str | Path) -> Path:
    return Path(out_dir).expanduser().resolve() / "logs" / "app_benchmark.log"


def read_manifest(out_dir: str | Path) -> dict[str, Any]:
    path = manifest_path(out_dir)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception as exc:
        return {"status": "manifest_error", "error": str(exc), "manifest_path": str(path)}


def tail_log(out_dir: str | Path, max_lines: int = 80) -> str:
    candidates = [
        log_path(out_dir),
        Path(out_dir).expanduser().resolve() / "logs" / "benchmark.log",
        Path(out_dir).expanduser().resolve() / "logs" / "full_benchmark.screen.log",
    ]
    for path in candidates:
        if path.exists():
            try:
                lines = path.read_text(errors="replace").splitlines()
                return "\n".join(lines[-max_lines:])
            except Exception as exc:
                return f"Could not read log {path}: {exc}"
    return "No benchmark log found yet."


def process_is_running(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    try:
        stat = subprocess.run(
            ["ps", "-p", str(pid), "-o", "stat="],
            capture_output=True,
            text=True,
            timeout=2,
        ).stdout.strip()
        if stat.startswith("Z"):
            return False
    except Exception:
        pass
    return True


def current_pid(out_dir: str | Path) -> int | None:
    path = pid_path(out_dir)
    if not path.exists():
        return None
    try:
        return int(path.read_text().strip())
    except Exception:
        return None


def benchmark_is_running(out_dir: str | Path) -> bool:
    manifest = read_manifest(out_dir)
    if manifest.get("status") in {"complete", "dry_run_complete", "failed"}:
        return False
    pid = current_pid(out_dir)
    return bool(pid and process_is_running(pid))


def launch_benchmark(
    *,
    app_dir: Path,
    fasta: str | Path,
    out_dir: str | Path,
    preset: str = "all",
    min_free_gb: float = 20.0,
    cpu: int = 4,
    nt_orf_mode: str = "sixframe",
    min_orf_aa: int = 30,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Start or resume a benchmark as a detached process."""
    app_dir = Path(app_dir).resolve()
    fasta_path = resolve_path(fasta, app_dir)
    out_path = resolve_path(out_dir, app_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    (out_path / "logs").mkdir(parents=True, exist_ok=True)

    if not fasta_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")

    pid = current_pid(out_path)
    if pid and process_is_running(pid):
        return {
            "status": "already_running",
            "pid": pid,
            "out": str(out_path),
            "log": str(log_path(out_path)),
        }

    cmd = [
        sys.executable,
        str(runner_script(app_dir)),
        "--preset",
        preset,
        "--fasta",
        str(fasta_path),
        "--out",
        str(out_path),
        "--min-free-gb",
        str(float(min_free_gb)),
        "--cpu",
        str(int(cpu)),
        "--nt-orf-mode",
        str(nt_orf_mode),
        "--min-orf-aa",
        str(int(min_orf_aa)),
    ]
    if dry_run:
        cmd.append("--dry-run")

    log_file = log_path(out_path)
    with log_file.open("ab") as fh:
        fh.write(("\n[app] launching: " + " ".join(cmd) + "\n").encode())
        proc = subprocess.Popen(
            cmd,
            cwd=str(app_dir),
            stdout=fh,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            close_fds=True,
            start_new_session=True,
        )

    pid_path(out_path).write_text(str(proc.pid))
    return {
        "status": "started",
        "pid": proc.pid,
        "out": str(out_path),
        "log": str(log_file),
        "dry_run": dry_run,
    }


def stop_benchmark(out_dir: str | Path) -> dict[str, Any]:
    pid = current_pid(out_dir)
    if not pid:
        return {"status": "not_running", "message": "No benchmark PID file found."}
    if not process_is_running(pid):
        return {"status": "not_running", "pid": pid}
    try:
        os.killpg(pid, signal.SIGTERM)
    except Exception:
        os.kill(pid, signal.SIGTERM)
    return {"status": "stopped", "pid": pid}
