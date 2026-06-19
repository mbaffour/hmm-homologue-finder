"""
core/storage_cleanup.py - Safe removal of bulky regenerable project files.

This module deliberately preserves final result tables, reports, figures, logs,
HMMs, alignments, input data, and reproducibility files. It removes only cache
and intermediate files that can be regenerated from the run settings.
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any


def _bytes(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        try:
            return path.stat().st_size
        except OSError:
            return 0
    total = 0
    for item in path.rglob("*"):
        if item.is_file():
            try:
                total += item.stat().st_size
            except OSError:
                pass
    return total


def _running_benchmark(project_dir: Path) -> tuple[bool, int | None]:
    try:
        from core.benchmark import current_pid, process_is_running

        pid = current_pid(project_dir)
        return bool(pid and process_is_running(pid)), pid
    except Exception:
        return False, None


def _candidate_paths(project_dir: Path, include_download_cache: bool) -> list[Path]:
    candidates: list[Path] = []
    for rel in [
        "runtime",
        "search_results/stream_cache",
        "search_results/translated",
        "results/synteny_context_cache",
        "tmp",
        "temp",
    ]:
        candidates.append(project_dir / rel)

    for parent_rel in ["search_results", "results"]:
        parent = project_dir / parent_rel
        if parent.exists():
            candidates.extend(parent.glob("*.gff"))
            candidates.extend(parent.glob("*.domtblout"))
            candidates.extend(parent.glob("*_sixframe.faa"))
            candidates.extend(parent.glob("*_prodigal.faa"))
            candidates.extend(parent.glob("translated_all.faa"))

    if include_download_cache:
        candidates.append(project_dir / "cache")
        candidates.append(project_dir / "databases" / "cache")

    # Keep prepared annotation setup DBs out of this list: those are expensive reusable
    # indexed assets, not ordinary temporary chunks.
    return [p for p in candidates if p.exists()]


def cleanup_preview(project_dir: str | Path, include_download_cache: bool = True) -> dict[str, Any]:
    root = Path(project_dir).expanduser().resolve()
    running, pid = _running_benchmark(root)
    rows = []
    total = 0
    for path in _candidate_paths(root, include_download_cache):
        size = _bytes(path)
        if size <= 0:
            continue
        total += size
        rows.append(
            {
                "path": str(path.relative_to(root) if path.is_relative_to(root) else path),
                "bytes": size,
                "kind": "directory" if path.is_dir() else "file",
            }
        )
    rows.sort(key=lambda r: r["bytes"], reverse=True)
    return {
        "project_dir": str(root),
        "running": running,
        "pid": pid,
        "include_download_cache": include_download_cache,
        "total_bytes": total,
        "items": rows,
    }


def cleanup_project(
    project_dir: str | Path,
    *,
    include_download_cache: bool = True,
    dry_run: bool = True,
) -> dict[str, Any]:
    preview = cleanup_preview(project_dir, include_download_cache=include_download_cache)
    if preview["running"]:
        preview["status"] = "blocked_running"
        preview["message"] = f"Benchmark process {preview['pid']} is still running; cleanup is disabled."
        return preview
    if dry_run:
        preview["status"] = "preview"
        return preview

    root = Path(project_dir).expanduser().resolve()
    removed = []
    freed = 0
    for item in preview["items"]:
        path = root / item["path"]
        size = _bytes(path)
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)
        else:
            path.unlink(missing_ok=True)
        if not path.exists():
            removed.append(item)
            freed += size
    return {
        **preview,
        "status": "cleaned",
        "removed": removed,
        "freed_bytes": freed,
    }


def format_bytes(num: int | float) -> str:
    value = float(num or 0)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if value < 1024 or unit == "TB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    return f"{value:.1f} TB"
