#!/usr/bin/env python3
"""
manage_cache.py — inspect / prune the shared database cache.

The pipeline keeps three kinds of data under the shared cache (default
~/.cache/hmm-homologue-finder): downloaded DB files, six-frame ORF translation
caches (<db>.sixframe.min<N>.faa, ~5 GB each, regenerable), and the VOGDB VFAM
annotation DB (~7 GB). On a server these accumulate, so this gives a quick size
report and safe pruning.

  python3 manage_cache.py                      # report sizes (default)
  python3 manage_cache.py --prune-translations # delete regenerable ORF caches
  python3 manage_cache.py --prune-vogdb        # delete the VOGDB annotation DB
  python3 manage_cache.py --prune-downloads    # delete downloaded DB files (re-fetched on use)
  python3 manage_cache.py --clear-all          # delete the ENTIRE cache (asks first on a terminal)
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path


def _size(p: Path) -> int:
    if p.is_file():
        return p.stat().st_size
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())


def _human(n: int) -> str:
    x = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if x < 1024 or unit == "TB":
            return f"{x:.1f} {unit}"
        x /= 1024


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", type=Path, default=Path.home() / ".cache" / "hmm-homologue-finder")
    ap.add_argument("--prune-translations", action="store_true",
                    help="delete six-frame ORF translation caches (regenerated on next run)")
    ap.add_argument("--prune-vogdb", action="store_true",
                    help="delete the VOGDB VFAM annotation DB (re-downloaded when needed)")
    ap.add_argument("--prune-downloads", action="store_true",
                    help="delete downloaded database files (re-downloaded on the next run)")
    ap.add_argument("--clear-all", action="store_true",
                    help="delete the ENTIRE shared cache (downloads + translations + VOGDB)")
    ap.add_argument("--yes", action="store_true",
                    help="skip the confirmation prompt for --clear-all")
    args = ap.parse_args()
    cache = args.cache.expanduser()
    if not cache.exists():
        print(f"No cache at {cache}")
        return

    trans = sorted(cache.rglob("*.sixframe.*.faa"))
    vogdb = cache / "annotation" / "vogdb"
    dbfiles = cache / "cache"

    print(f"Cache: {cache}  (total {_human(_size(cache))})")
    print(f"  downloaded DBs:        {_human(_size(dbfiles)) if dbfiles.exists() else '0 B'}")
    print(f"  translation caches:    {_human(sum(_size(t) for t in trans))}  ({len(trans)} file(s))")
    print(f"  VOGDB annotation DB:   {_human(_size(vogdb)) if vogdb.exists() else '0 B'}")

    if args.prune_translations:
        freed = sum(_size(t) for t in trans)
        for t in trans:
            t.unlink(missing_ok=True)
        print(f"Pruned {len(trans)} translation cache(s), freed {_human(freed)}.")
    if args.prune_vogdb and vogdb.exists():
        freed = _size(vogdb)
        shutil.rmtree(vogdb, ignore_errors=True)
        print(f"Removed VOGDB annotation DB, freed {_human(freed)}.")
    if args.prune_downloads and dbfiles.exists():
        freed = _size(dbfiles)
        shutil.rmtree(dbfiles, ignore_errors=True)
        print(f"Removed downloaded database files, freed {_human(freed)}.")
    if args.clear_all:
        total = _size(cache)
        if not args.yes and sys.stdin.isatty():
            resp = input(f"Delete the ENTIRE cache at {cache} ({_human(total)})? [y/N] ").strip().lower()
            if resp not in ("y", "yes"):
                print("Aborted — nothing deleted.")
                return
        shutil.rmtree(cache, ignore_errors=True)
        print(f"Cleared the whole cache, freed {_human(total)}.")


if __name__ == "__main__":
    main()
