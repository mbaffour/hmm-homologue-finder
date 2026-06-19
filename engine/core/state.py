"""
core/state.py — PipelineState: persists step completion to .pipeline_state.json.

Usage:
    state = PipelineState(proj_dir)
    state.mark_complete("hmm_build", {"hmm_path": str(hmm)})
    if state.is_complete("hmm_build"):
        ...
    params = state.get_params("hmm_build")
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Optional


class PipelineState:
    STEPS = [
        "input",
        "msa",
        "hmm_build",
        "search",
        "validate",
        "iterate",
        "classify",
        "synteny",
        "taxonomy",
        "phylo",
        "matrix",
        "clusters",
        "motifs",
        "annotation",
        "export",
    ]

    def __init__(self, proj_dir: Path):
        self.proj_dir = Path(proj_dir)
        self._state_file = self.proj_dir / ".pipeline_state.json"
        self._data: dict[str, Any] = self._load()

    def _load(self) -> dict[str, Any]:
        if self._state_file.exists():
            try:
                return json.loads(self._state_file.read_text())
            except Exception:
                pass
        return {"steps": {}, "project": {}, "metadata": {}}

    def _save(self) -> None:
        self._data["metadata"]["last_updated"] = datetime.now().isoformat()
        self._state_file.write_text(json.dumps(self._data, indent=2))

    # ------------------------------------------------------------------
    # Step tracking
    # ------------------------------------------------------------------

    def mark_complete(self, step: str, params: Optional[dict] = None) -> None:
        self._data["steps"][step] = {
            "status": "complete",
            "completed_at": datetime.now().isoformat(),
            "params": params or {},
        }
        self._save()

    def mark_failed(self, step: str, error: str = "") -> None:
        self._data["steps"][step] = {
            "status": "failed",
            "failed_at": datetime.now().isoformat(),
            "error": error,
        }
        self._save()

    def mark_running(self, step: str) -> None:
        self._data["steps"][step] = {
            "status": "running",
            "started_at": datetime.now().isoformat(),
        }
        self._save()

    def is_complete(self, step: str) -> bool:
        return self._data["steps"].get(step, {}).get("status") == "complete"

    def is_running(self, step: str) -> bool:
        return self._data["steps"].get(step, {}).get("status") == "running"

    def get_status(self, step: str) -> str:
        return self._data["steps"].get(step, {}).get("status", "pending")

    def get_params(self, step: str) -> dict:
        return self._data["steps"].get(step, {}).get("params", {})

    def get_all_steps(self) -> dict[str, dict]:
        return {
            step: {
                "status": self.get_status(step),
                **self._data["steps"].get(step, {}),
            }
            for step in self.STEPS
        }

    def reset_step(self, step: str) -> None:
        self._data["steps"].pop(step, None)
        self._save()

    # ------------------------------------------------------------------
    # Project metadata
    # ------------------------------------------------------------------

    def set_project(self, key: str, value: Any) -> None:
        self._data["project"][key] = value
        self._save()

    def get_project(self, key: str, default: Any = None) -> Any:
        return self._data["project"].get(key, default)

    def set_input(
        self,
        input_path: str,
        seq_type: str,
        seq_count: int,
        mode: str = "generic",
    ) -> None:
        self._data["project"].update(
            {
                "input_path": input_path,
                "seq_type": seq_type,
                "seq_count": seq_count,
                "biology_mode": mode,
            }
        )
        self._save()
