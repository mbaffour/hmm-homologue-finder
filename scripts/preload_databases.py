#!/usr/bin/env python3
"""preload_databases.py - proactively DOWNLOAD the search databases and build the SIX-FRAME
TRANSLATION cache ahead of time, so later discovery runs spend their time searching rather than
downloading and translating.

Both the database downloads and the six-frame ORF translation of nucleotide databases are cached
PERSISTENTLY under the shared cache (<db-cache>/cache/...). The engine builds these during any run
with --keep-cache and REUSES them on every later iteration AND every future run - e.g. the INPHARED
six-frame translation drops from ~40 min the first time to ~20 s once cached. This command warms
those caches now, before you have your real seeds, by running a MINIMAL pass over the selected
databases driven by the bundled demo seeds: the engine downloads each database and six-frame-
translates the nucleotide ones into the shared cache (the demo seeds find ~nothing - that doesn't
matter, the download + translation caches are built regardless). It is IDEMPOTENT: anything already
cached is reused, not rebuilt, so you can run it whenever you have spare time/bandwidth.

It does NOT modify the engine - it warms exactly the caches a real discovery run consumes, by using
the same engine + the same --keep-cache mechanism, so the cached ORFs are guaranteed reusable.

USAGE (use the conda env's python so the tools are on PATH, no `conda activate` needed):
    <env>/bin/python scripts/preload_databases.py [--databases "A,B"] [--db-cache DIR]
        [--min-orf-aa N] [--cpu N] [--seed FASTA] [--dry-run]
Default databases = the full default discovery set; default cache = ~/.cache/hmm-homologue-finder.
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEPLOY = HERE.parent
BENCHMARK = DEPLOY / "engine" / "scripts" / "run_all_database_benchmark.py"
DEMO_SEED = DEPLOY / "examples" / "example_seeds.fasta"

# The full default discovery database set, kept in sync with hmm_finder.DATABASES.
try:
    sys.path.insert(0, str(HERE))
    from hmm_finder import DATABASES as DEFAULT_DATABASES
except Exception:
    DEFAULT_DATABASES = ("INPHARED genomes,INPHARED proteins,SwissProt,RefSeq viral proteins,"
                         "RefSeq viral genomes,Gut Phage Database (GPD),GVD-AVrC,"
                         "Pfam (sequences),Pfam (domain scan),VOGDB VFAM (annotation)")


def translation_cache(shared: Path, min_aa: int) -> dict:
    """{relpath: size_bytes} of every six-frame translation cache file under <shared>/cache.
    This is exactly the cache a discovery run reuses (`<dbfile>.sixframe.min<N>.faa`)."""
    cdir = shared / "cache"
    out = {}
    if cdir.exists():
        for f in cdir.rglob(f"*.sixframe.min{min_aa}.faa"):
            try:
                out[str(f.relative_to(cdir))] = f.stat().st_size
            except OSError:
                pass
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--databases", default=DEFAULT_DATABASES,
                    help="comma-separated database names to preload (default: the full discovery set; "
                         "run `python scripts/hmm_finder.py --list-databases` for the catalog)")
    ap.add_argument("--db-cache", type=Path, default=Path.home() / ".cache" / "hmm-homologue-finder",
                    help="persistent shared cache dir (default ~/.cache/hmm-homologue-finder)")
    ap.add_argument("--min-orf-aa", type=int, default=30,
                    help="min six-frame ORF length to translate/cache (default 30, matches the search floor)")
    ap.add_argument("--cpu", default="4")
    ap.add_argument("--seed", type=Path, default=DEMO_SEED,
                    help="seed FASTA that DRIVES the prep pass (any valid FASTA; default: bundled demo seeds). "
                         "Its hits are irrelevant - the download + six-frame caches are built regardless.")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the resolved command + current cache state, then exit without downloading")
    args = ap.parse_args()

    shared = args.db_cache.expanduser().resolve()
    (shared / "cache").mkdir(parents=True, exist_ok=True)
    (shared / "db_setup").mkdir(parents=True, exist_ok=True)
    if not args.seed.exists():
        sys.exit(f"[preload] seed FASTA not found: {args.seed}")

    before = translation_cache(shared, args.min_orf_aa)
    print(f"[preload] shared cache : {shared}")
    print(f"[preload] databases    : {args.databases}")
    print(f"[preload] six-frame translation caches already present (min{args.min_orf_aa}): {len(before)}")

    # A scratch run dir INSIDE the shared cache, with cache/db_setup symlinked to the persistent
    # dirs (exactly how hmm_finder wires the engine), so downloads + translations land in the
    # shared cache and survive when we delete the scratch run dir.
    import tempfile
    prep = Path(tempfile.mkdtemp(prefix="_preload_run_", dir=str(shared)))
    for name in ("cache", "db_setup"):
        link = prep / name
        if not link.exists():
            link.symlink_to(shared / name)
    cmd = [sys.executable, str(BENCHMARK), "--fasta", str(args.seed), "--out", str(prep),
           "--databases", args.databases, "--cpu", str(args.cpu), "--keep-cache", "--skip-tree",
           "--min-recovery", "0", "--min-orf-aa", str(args.min_orf_aa)]

    if args.dry_run:
        print("[preload] (dry-run) would run:\n  " + " ".join(map(str, cmd)))
        shutil.rmtree(prep, ignore_errors=True)
        return 0

    env = dict(os.environ)
    env["PATH"] = os.path.dirname(sys.executable) + os.pathsep + env.get("PATH", "")
    print("[preload] warming caches (downloads + six-frame translation - the slow part, done once):")
    print("  " + " ".join(map(str, cmd)), flush=True)
    rc = subprocess.call(cmd, env=env, stdin=subprocess.DEVNULL)

    after = translation_cache(shared, args.min_orf_aa)
    new = {k for k in after if k not in before}
    print(f"\n[preload] six-frame translation caches now: {len(after)} ({len(new)} new this run)")
    for k in sorted(after):
        print(f"  {'NEW ' if k in new else 'have'} {after[k] / 1e6:9.1f} MB  {k}")
    print("[preload] Done. Later discovery runs over these databases reuse the downloads and the "
          "six-frame translation - no re-download, no re-translation.")
    shutil.rmtree(prep, ignore_errors=True)   # the warmed caches live in <shared>/cache, not here
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
