"""
core/sessions.py — Session notes and recent projects helpers.

Session persistence in this app works at the project-directory level:
the pipeline state JSON already tracks step completion. Session notes
add a lightweight human-readable layer so users know where they left off.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any

# File that stores the list of recently-opened project paths
_RECENTS_FILE = Path.home() / ".hmm_discovery_recents.json"
_MAX_RECENTS = 10


# ---------------------------------------------------------------------------
# Recent projects
# ---------------------------------------------------------------------------

def load_recents() -> list[str]:
    """Return list of recently-used project paths (newest first)."""
    try:
        data = json.loads(_RECENTS_FILE.read_text())
        return [p for p in data if Path(p).exists()][:_MAX_RECENTS]
    except Exception:
        return []


def add_recent(proj_path: str) -> None:
    """Prepend *proj_path* to the recents list and save."""
    recents = load_recents()
    # Deduplicate while preserving order
    recents = [proj_path] + [p for p in recents if p != proj_path]
    try:
        _RECENTS_FILE.write_text(json.dumps(recents[:_MAX_RECENTS], indent=2))
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Session notes (stored inside the project folder)
# ---------------------------------------------------------------------------

def _notes_file(proj_dir: Path) -> Path:
    logs = proj_dir / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    return logs / "session_notes.jsonl"


def save_note(proj_dir: Path, note: str, step_statuses: dict) -> None:
    """Append a timestamped session note to the project's notes file."""
    entry = {
        "timestamp": datetime.now().isoformat(),
        "note": note.strip(),
        "steps": step_statuses,
    }
    nf = _notes_file(proj_dir)
    with nf.open("a") as fh:
        fh.write(json.dumps(entry) + "\n")


def load_notes(proj_dir: Path, last_n: int = 5) -> list[dict]:
    """Return the most recent *last_n* session notes, newest first."""
    nf = _notes_file(proj_dir)
    if not nf.exists():
        return []
    notes: list[dict] = []
    try:
        for line in nf.read_text().splitlines():
            line = line.strip()
            if line:
                try:
                    notes.append(json.loads(line))
                except Exception:
                    pass
    except Exception:
        pass
    return list(reversed(notes))[:last_n]
