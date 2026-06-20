#!/usr/bin/env python3
"""
db_catalog.py — read the searchable-database catalog from the bundled engine and
provide terminal helpers to list databases and pick them interactively.

The catalog (names + metadata) lives in engine/databases/builtin.py. This module
imports it on demand and degrades gracefully (name-only fallback) if the engine
cannot be imported.
"""
from __future__ import annotations

import sys
from pathlib import Path


def load_catalog(engine_dir: Path) -> list[dict]:
    """Return the list of database dicts from engine/databases/builtin.py.

    Returns [] if the engine package cannot be imported.
    """
    engine_dir = Path(engine_dir)
    try:
        if str(engine_dir) not in sys.path:
            sys.path.insert(0, str(engine_dir))
        from databases.builtin import BUILTIN_DATABASES  # type: ignore
        return list(BUILTIN_DATABASES)
    except Exception as e:  # pragma: no cover - defensive
        print(f"(could not load database catalog from {engine_dir}: {e})", file=sys.stderr)
        return []


def _fmt(db: dict) -> str:
    typ = db.get("type", "?")
    bits = [f"[{typ}]"]
    if db.get("size_hint"):
        bits.append(str(db["size_hint"]))
    if db.get("est_time"):
        bits.append(str(db["est_time"]))
    line = f"{db.get('name', '?')}  ({'  '.join(bits)})"
    rel = db.get("relevance") or db.get("notes")
    if rel:
        line += f"\n        {rel}"
    return line


def list_databases(engine_dir: Path) -> None:
    """Print the full catalog (for --list-databases)."""
    cat = load_catalog(engine_dir)
    if not cat:
        print("No database catalog available.")
        return
    print("\nAvailable databases (pass exact names to --databases):\n")
    for i, db in enumerate(cat, 1):
        print(f"  {i:2d}. {_fmt(db)}")
    print()


def pick_databases(engine_dir: Path, default_names: list[str]) -> str:
    """Interactive numbered multi-select; returns a comma-separated name string.

    Enter = defaults, 'all' = every database, otherwise space/comma-separated
    numbers. Falls back to the defaults if the catalog can't be loaded.
    """
    cat = load_catalog(engine_dir)
    if not cat:
        return ",".join(default_names)
    print("\nSelect databases to search:")
    print("  Enter = defaults (*)   |   'all' = everything   |   else: numbers e.g. 2 4 7\n")
    for i, db in enumerate(cat, 1):
        mark = "*" if db.get("name") in default_names else " "
        print(f" {mark}{i:2d}. {_fmt(db)}")
    raw = input("\n  databases > ").strip()
    if not raw:
        return ",".join(default_names)
    if raw.lower() == "all":
        return ",".join(db["name"] for db in cat)
    chosen: list[str] = []
    for tok in raw.replace(",", " ").split():
        if tok.isdigit() and 1 <= int(tok) <= len(cat):
            chosen.append(cat[int(tok) - 1]["name"])
        else:
            print(f"  (ignoring '{tok}')")
    if not chosen:
        print("  No valid selection; using defaults.")
        return ",".join(default_names)
    return ",".join(dict.fromkeys(chosen))  # de-dup, preserve order
