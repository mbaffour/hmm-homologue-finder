"""
core/runner.py — AsyncJobRunner: runs subprocesses with live log streaming.

Each step UI panel owns one AsyncJobRunner instance. Calling .start(cmd)
kicks off the process; the reactive log_lines value is updated line-by-line
so the UI updates in real time.
"""
from __future__ import annotations

import asyncio
import os
import time
from pathlib import Path
from typing import Optional

from shiny import reactive

from .logger import AuditLogger, file_hash


class AsyncJobRunner:
    """Non-blocking subprocess runner with reactive log streaming."""

    def __init__(self, step_name: str = "", audit: Optional[AuditLogger] = None):
        self.step_name = step_name
        self.audit = audit
        self.is_running: reactive.Value[bool] = reactive.value(False)
        self.log_lines: reactive.Value[list[str]] = reactive.value([])
        self.returncode: reactive.Value[Optional[int]] = reactive.value(None)
        self._task: Optional[asyncio.Task] = None
        self._start_time: float = 0.0

    def start(
        self,
        cmd: list,
        env: Optional[dict] = None,
        cwd: Optional[Path] = None,
        input_files: Optional[dict] = None,
        output_files: Optional[dict] = None,
    ) -> None:
        """Launch the command asynchronously."""
        import shutil as _shutil
        cmd = [str(c) for c in cmd]
        # Resolve the executable to its full path using the augmented PATH
        # so the subprocess never fails with "command not found"
        aug = self._augmented_path()
        if cmd and not os.path.isabs(cmd[0]):
            full = _shutil.which(cmd[0], path=aug)
            if full:
                cmd[0] = full
        self.log_lines.set([f"$ {' '.join(cmd)}", ""])
        self.returncode.set(None)
        self.is_running.set(True)
        self._start_time = time.monotonic()
        self._task = asyncio.ensure_future(
            self._run(cmd, env or {}, cwd, input_files or {}, output_files or {})
        )

    @staticmethod
    def _augmented_path() -> str:
        """Return PATH augmented with conda env bin dirs so tools like mafft/hmmbuild are found."""
        import shutil
        extras: list[str] = []
        home = Path.home()

        # Current conda env (set when running via conda run or conda activate)
        for var in ("CONDA_PREFIX", "VIRTUAL_ENV"):
            val = os.environ.get(var)
            if val:
                extras.append(str(Path(val) / "bin"))
                break

        # Fallback: common conda env names and locations
        for base in [home / "miniforge3", home / "miniconda3", home / "anaconda3",
                     Path("/opt/anaconda3"), Path("/opt/miniconda3"), Path("/opt/homebrew")]:
            for envname in ["hmm_env", "base", ""]:
                subdir = base / "envs" / envname / "bin" if envname else base / "bin"
                if subdir.is_dir() and shutil.which("mafft", path=str(subdir)):
                    extras.append(str(subdir))
                    break

        current = os.environ.get("PATH", "")
        all_paths = extras + [p for p in current.split(os.pathsep) if p not in extras]
        return os.pathsep.join(all_paths)

    async def _run(
        self,
        cmd: list[str],
        extra_env: dict,
        cwd: Optional[Path],
        input_files: dict,
        output_files: dict,
    ) -> int:
        augmented = self._augmented_path()
        merged_env = {**os.environ, "PATH": augmented, **extra_env}
        try:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                env=merged_env,
                cwd=str(cwd) if cwd else None,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT,
            )
            async for raw in proc.stdout:  # type: ignore[union-attr]
                line = raw.decode(errors="replace").rstrip()
                current = self.log_lines.get()
                self.log_lines.set(current + [line])
            rc = await proc.wait()
        except Exception as exc:
            self.log_lines.set(self.log_lines.get() + [f"[ERROR] {exc}"])
            rc = 1

        duration = time.monotonic() - self._start_time
        self.returncode.set(rc)
        self.is_running.set(False)

        if self.audit:
            in_hashes = {k: file_hash(Path(v)) for k, v in input_files.items() if v}
            out_hashes = {k: file_hash(Path(v)) for k, v in output_files.items() if v}
            self.audit.record(
                step=self.step_name,
                command=cmd,
                exit_code=rc,
                duration_sec=duration,
                input_files=in_hashes,
                output_files=out_hashes,
            )
        return rc

    def cancel(self) -> None:
        if self._task and not self._task.done():
            self._task.cancel()
        self.is_running.set(False)
        self.log_lines.set(self.log_lines.get() + ["[CANCELLED]"])

    def succeeded(self) -> bool:
        rc = self.returncode.get()
        return rc is not None and rc == 0

    def get_log(self) -> str:
        return "\n".join(self.log_lines.get())
