"""
databases/registry.py — DatabaseRegistry: manages databases.json per project.

Stores database metadata: name, type (protein/nucleotide), path/url,
streaming (bool), enabled (bool), download_url, size_bytes, notes.
"""

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from .builtin import BUILTIN_DATABASES

logger = logging.getLogger(__name__)

_DB_FILENAME = "databases.json"


class DatabaseRegistry:
    """
    Manages the set of databases available for HMM discovery searches.

    Persists state to ``<proj_dir>/databases.json``. On first use,
    populates the file from the built-in database catalogue.
    """

    def __init__(self, proj_dir: Path) -> None:
        """
        Load the registry from ``<proj_dir>/databases.json``.

        If the file does not exist the registry is seeded from
        :data:`~databases.builtin.BUILTIN_DATABASES` and saved immediately.

        Parameters
        ----------
        proj_dir:
            Root directory of the active project.  The file
            ``databases.json`` will be read from / written to this
            directory.
        """
        self._proj_dir = Path(proj_dir)
        self._db_file = self._proj_dir / _DB_FILENAME
        self._databases: dict[str, dict] = {}  # keyed by name

        if self._db_file.exists():
            self._load()
        else:
            self._seed_from_builtins()
            self.save()

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def _load(self) -> None:
        """Read databases.json into memory, then merge any new builtins."""
        try:
            raw = json.loads(self._db_file.read_text(encoding="utf-8"))
            if not isinstance(raw, list):
                raise ValueError("databases.json must be a JSON array")
            self._databases = {entry["name"]: entry for entry in raw}
            logger.debug("Loaded %d databases from %s", len(self._databases), self._db_file)
        except (json.JSONDecodeError, KeyError, ValueError) as exc:
            logger.warning("Could not parse %s: %s — reseeding from builtins", self._db_file, exc)
            self._seed_from_builtins()
            self.save()
            return

        # Merge new/updated builtins into existing registry.
        # - Add databases that are in builtins but missing from the file.
        # - Update URL/streaming/notes for existing ones (preserve user's path/enabled).
        changed = False
        for builtin in BUILTIN_DATABASES:
            name = builtin["name"]
            if name not in self._databases:
                # Brand new database — add it
                self._databases[name] = dict(builtin)
                changed = True
                logger.debug("Added new builtin database: %s", name)
            else:
                # Existing database — update URL/streaming/notes but keep user overrides
                existing = self._databases[name]
                for key in ("download_url", "url", "streaming", "notes", "relevance", "size_hint",
                            "est_time", "search_mode", "setup_handler", "annotation_url",
                            "release", "checksum_md5", "annotation_checksum_md5"):
                    if key in builtin and builtin[key] != existing.get(key):
                        existing[key] = builtin[key]
                        changed = True
                # Re-enable if builtin says enabled and user hasn't explicitly disabled
                if builtin.get("enabled", True) and not existing.get("enabled", True):
                    if not existing.get("_user_disabled"):
                        existing["enabled"] = True
                        changed = True

        # Remove databases that are no longer in builtins (name changed, etc.)
        # but only if they have no local path set (user hasn't customized them)
        builtin_names = {db["name"] for db in BUILTIN_DATABASES}
        stale = [
            name for name in self._databases
            if name not in builtin_names
            and not self._databases[name].get("path")
            and self._databases[name].get("download_url", "").startswith("http")
        ]
        for name in stale:
            del self._databases[name]
            changed = True
            logger.debug("Removed stale database: %s", name)

        # Seeds are inputs/positive controls, not searchable databases. Older
        # app builds could persist seed FASTA files as custom DB entries; remove
        # those automatically so projects stay aligned with the built-in DB
        # registry unless the user explicitly registers a non-seed custom DB.
        seed_like = [
            name for name in self._databases
            if name.lower() in {"seed", "seeds", "seed_sequences"}
            or name.lower().endswith("_seed_sequences")
        ]
        for name in seed_like:
            del self._databases[name]
            changed = True
            logger.debug("Removed seed input from database registry: %s", name)

        if changed:
            self.save()

    def _seed_from_builtins(self) -> None:
        """Populate registry from the built-in catalogue."""
        self._databases = {db["name"]: dict(db) for db in BUILTIN_DATABASES}
        logger.debug("Seeded registry with %d built-in databases", len(self._databases))

    def save(self) -> None:
        """Persist the current registry state to ``databases.json``."""
        self._proj_dir.mkdir(parents=True, exist_ok=True)
        data = list(self._databases.values())
        self._db_file.write_text(
            json.dumps(data, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        logger.debug("Saved %d databases to %s", len(data), self._db_file)

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def get_all(self) -> list[dict]:
        """Return a list of all registered database records (copies)."""
        return [dict(db) for db in self._databases.values()]

    def get_enabled(self) -> list[dict]:
        """Return only the databases that have ``enabled=True``."""
        return [dict(db) for db in self._databases.values() if db.get("enabled", False)]

    def get(self, name: str) -> Optional[dict]:
        """
        Return the database record for *name*, or ``None`` if not found.

        Parameters
        ----------
        name:
            Exact database name as stored in the registry.
        """
        entry = self._databases.get(name)
        return dict(entry) if entry is not None else None

    # ------------------------------------------------------------------
    # Mutations
    # ------------------------------------------------------------------

    def register(
        self,
        name: str,
        db_type: str,
        path: Optional[str] = None,
        url: Optional[str] = None,
        streaming: bool = False,
        notes: str = "",
        download_url: Optional[str] = None,
        optional: bool = True,
        size_hint: str = "",
        enabled: bool = True,
    ) -> None:
        """
        Add or overwrite a database entry.

        If a database with *name* already exists it is replaced in full.

        Parameters
        ----------
        name:
            Unique display name for the database.
        db_type:
            ``"protein"`` or ``"nucleotide"``.
        path:
            Absolute path to a local FASTA / HMM file, or ``None``.
        url:
            Streaming base URL, or ``None``.
        streaming:
            ``True`` if the database is accessed remotely without a full
            local download.
        notes:
            Free-text description.
        download_url:
            URL from which the database can be downloaded on demand.
        optional:
            ``True`` if a missing local file should only show a badge
            rather than blocking searches.
        size_hint:
            Human-readable size estimate, e.g. ``"~700 MB"``.
        enabled:
            Whether the database is active by default.
        """
        if db_type not in ("protein", "nucleotide"):
            raise ValueError(f"db_type must be 'protein' or 'nucleotide', got {db_type!r}")

        # Preserve per-run stats if the entry already exists.
        existing = self._databases.get(name, {})
        self._databases[name] = {
            "name": name,
            "type": db_type,
            "path": path,
            "url": url,
            "download_url": download_url,
            "streaming": streaming,
            "enabled": enabled,
            "optional": optional,
            "size_hint": size_hint,
            "notes": notes,
            "relevance": existing.get("relevance", ""),
            "est_time": existing.get("est_time", ""),
            "search_mode": existing.get("search_mode", "hmmsearch"),
            "setup_handler": existing.get("setup_handler"),
            "annotation_url": existing.get("annotation_url"),
            "release": existing.get("release", ""),
            "checksum_md5": existing.get("checksum_md5", ""),
            "annotation_checksum_md5": existing.get("annotation_checksum_md5", ""),
            "last_searched": existing.get("last_searched"),
            "last_hit_count": existing.get("last_hit_count"),
        }
        self.save()

    def set_enabled(self, name: str, enabled: bool) -> None:
        """
        Enable or disable a database by name.

        Parameters
        ----------
        name:
            Database name.
        enabled:
            New enabled state.

        Raises
        ------
        KeyError
            If *name* is not registered.
        """
        self._require(name)
        self._databases[name]["enabled"] = enabled
        self.save()

    def update_path(self, name: str, path: str) -> None:
        """Alias for :meth:`set_path` used by inline download handler."""
        if name in self._databases:
            self._databases[name]["path"] = path
            self.save()

    def set_path(self, name: str, path: str) -> None:
        """
        Update the local file path for a database (e.g. after download).

        Parameters
        ----------
        name:
            Database name.
        path:
            Absolute path to the downloaded file.

        Raises
        ------
        KeyError
            If *name* is not registered.
        """
        self._require(name)
        self._databases[name]["path"] = path
        self.save()

    def remove(self, name: str) -> None:
        """
        Remove a database from the registry.

        Does nothing if *name* is not registered (idempotent).
        """
        if name in self._databases:
            del self._databases[name]
            self.save()

    def record_search(self, name: str, hit_count: int) -> None:
        """
        Update ``last_searched`` timestamp and ``last_hit_count`` for *name*.

        Parameters
        ----------
        name:
            Database name.
        hit_count:
            Number of hits returned in the most recent search.
        """
        if name not in self._databases:
            return
        self._databases[name]["last_searched"] = datetime.now(timezone.utc).isoformat()
        self._databases[name]["last_hit_count"] = hit_count
        self.save()

    # ------------------------------------------------------------------
    # Availability & status
    # ------------------------------------------------------------------

    def is_available(self, name: str) -> bool:
        """
        Return ``True`` if the database can be used for a search.

        A database is considered available when:

        * It has a ``path`` pointing to an existing local file, **or**
        * ``streaming=True`` (remote access, no local file required).

        Parameters
        ----------
        name:
            Database name.  Returns ``False`` for unknown names.
        """
        entry = self._databases.get(name)
        if entry is None:
            return False
        if entry.get("streaming"):
            return True
        local = entry.get("path")
        if local and Path(local).exists():
            return True
        return False

    def status_badge(self, name: str) -> str:
        """
        Return a short status string for display in the UI.

        Possible return values:

        ``"local"``
            A local file exists and is ready to use.
        ``"streaming"``
            The database is accessed remotely; no download required.
        ``"download_available"``
            No local file, not streaming, but a ``download_url`` is set.
        ``"not_configured"``
            No local path, not streaming, and no download URL.

        Parameters
        ----------
        name:
            Database name.  Returns ``"not_configured"`` for unknown names.
        """
        entry = self._databases.get(name)
        if entry is None:
            return "not_configured"

        local = entry.get("path")
        if local and Path(local).exists():
            return "local"

        if entry.get("streaming"):
            return "streaming"

        if entry.get("download_url"):
            return "download_available"

        return "not_configured"

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _require(self, name: str) -> None:
        """Raise ``KeyError`` if *name* is not in the registry."""
        if name not in self._databases:
            raise KeyError(f"Database {name!r} is not registered")

    # ------------------------------------------------------------------
    # Compatibility aliases
    # ------------------------------------------------------------------

    def list_all(self) -> list[dict]:
        """Alias for :meth:`get_all` kept for UI compatibility."""
        return self.get_all()

    def __len__(self) -> int:
        return len(self._databases)

    def __repr__(self) -> str:
        return (
            f"DatabaseRegistry(proj_dir={self._proj_dir!r}, "
            f"n_databases={len(self._databases)})"
        )
