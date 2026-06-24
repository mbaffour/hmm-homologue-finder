#!/usr/bin/env python3
"""run_pipeline.py - run the pipeline from Python, fully autonomously, with ZERO interactive input.

A thin, no-prompt launcher around hmm_finder.py (family discovery) and scan_genome.py (single-genome
scan). It:
  1. puts the running interpreter's own bin dir on PATH, so the conda env tools (hmmsearch, mafft,
     prodigal, iqtree, ...) are found WITHOUT `conda activate`;
  2. runs with stdin detached (DEVNULL) so no step can ever block on a prompt;
  3. injects no-prompt defaults (full database set; offline unless --email) and forwards every other
     flag straight through.

USAGE - invoke with the conda env's python (so the env's tools are on PATH):
    <env>/bin/python scripts/run_pipeline.py --fasta SEED.faa --out-dir DIR [--email you@inst.edu] [...]

  Family discovery (default):
    ~/miniforge3/envs/hmm-discovery/bin/python scripts/run_pipeline.py \
        --fasta seed.faa --find-interrupted --out-dir ~/hmm_runs/myrun --iterations 2 --cpu 4
  Single genome ("does THIS genome carry my gene?") - add --scan:
    ... run_pipeline.py --scan --hmm gene.hmm --accession NC_008720 --email you@inst.edu --out ~/scan1

LAUNCHER FLAGS (consumed here, not forwarded):
  --scan              route to scan_genome.py (single-genome mode) instead of hmm_finder.py
  --preset NAME       apply a named bundle of defaults (see --list-presets)
  --list-presets      print the presets and exit
  --dry-run           print the exact resolved command (with injected defaults) and exit, run nothing
  -h, --help          print this help and exit

Only required flag: --fasta (discovery) or --hmm/--seeds + --genome/--accession (scan).
Organism names from NCBI need --email; otherwise it runs fully offline. Any hmm_finder.py /
scan_genome.py flag works (run those with --help for the full list). Exits 0 on success, non-zero on
failure; on a real run it also writes the resolved command to <out-dir>/run_command.txt.
"""
import os
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Named default bundles. The user can always override any flag by passing it explicitly.
PRESETS = {
    "phage-discovery": ["--all-databases", "--find-interrupted"],   # the typical overprinting search
    "discovery":       ["--all-databases"],                         # clean family discovery
    "offline":         ["--all-databases", "--no-annotate"],        # no NCBI at all
    "smoke":           ["--smoke"],                                 # fast plumbing self-test
}


def _has(args, *names):
    return any(a == n or a.startswith(n + "=") for a in args for n in names)


def _out_dir(args):
    """The value of --out-dir/--out if present (so we can write run_command.txt there)."""
    for flag in ("--out-dir", "--out"):
        for i, a in enumerate(args):
            if a == flag and i + 1 < len(args):
                return args[i + 1]
            if a.startswith(flag + "="):
                return a.split("=", 1)[1]
    return None


def _validate_fasta(path, what):
    p = Path(path)
    if not p.exists():
        sys.exit(f"[run_pipeline] {what} not found: {path}")
    try:
        head = p.read_text(errors="ignore").lstrip()[:1]
    except Exception as e:
        sys.exit(f"[run_pipeline] could not read {what} {path}: {e}")
    if head != ">":
        sys.exit(f"[run_pipeline] {what} does not look like FASTA (no '>' header): {path}")


def main() -> int:
    args = list(sys.argv[1:])

    if not args or _has(args, "-h", "--help"):
        print(__doc__)
        return 0
    if "--list-presets" in args:
        print("presets (--preset NAME):")
        for k, v in PRESETS.items():
            print(f"  {k:16s} {' '.join(v)}")
        return 0

    # launcher-level flags (consumed, not forwarded)
    scan_mode = "--scan" in args
    if scan_mode:
        args.remove("--scan")
    dry_run = "--dry-run" in args
    if dry_run:
        args.remove("--dry-run")
    if "--preset" in args:
        i = args.index("--preset")
        if i + 1 >= len(args) or args[i + 1] not in PRESETS:
            sys.exit(f"[run_pipeline] --preset needs one of: {', '.join(PRESETS)}")
        preset = args[i + 1]
        del args[i:i + 2]
        if preset == "smoke":
            scan_mode = False
        args = PRESETS[preset] + args            # preset first; explicit user flags still win

    target = HERE / ("scan_genome.py" if scan_mode else "hmm_finder.py")

    # no-prompt defaults (only if the user did not already choose) + light validation
    if scan_mode:
        if not _has(args, "--genome", "--accession"):
            sys.exit("[run_pipeline] --scan needs --genome PATH or --accession ACC")
        if not _has(args, "--hmm", "--seeds"):
            sys.exit("[run_pipeline] --scan needs --hmm FILE or --seeds FASTA")
        # scan_genome already exits cleanly on a missing email for --accession; nothing to inject
    else:
        if not _has(args, "--fasta"):
            sys.exit("[run_pipeline] discovery needs --fasta SEED.faa (or use --scan for one genome)")
        for i, a in enumerate(args):
            if a == "--fasta" and i + 1 < len(args):
                _validate_fasta(args[i + 1], "seed FASTA")
        if not _has(args, "--databases", "--all-databases", "--pick-databases", "--list-databases", "--smoke"):
            args.append("--all-databases")          # full default set, no prompt
        if not _has(args, "--email", "--no-annotate"):
            args.append("--no-annotate")            # offline unless an email is given
        if not _has(args, "--skip-tool-check"):
            args.append("--skip-tool-check")

    cmd = [sys.executable, str(target)] + args
    pretty = shlex.join(cmd)
    if dry_run:
        print("[run_pipeline] (dry-run) would run:\n  " + pretty)
        return 0

    # 1) make the conda env's binaries discoverable (this interpreter lives in <env>/bin).
    env = dict(os.environ)
    env["PATH"] = os.path.dirname(sys.executable) + os.pathsep + env.get("PATH", "")
    # 2) record the exact command in the output dir (before running) for provenance.
    od = _out_dir(args)
    if od:
        try:
            Path(od).mkdir(parents=True, exist_ok=True)
            (Path(od) / "run_command.txt").write_text(pretty + "\n", encoding="utf-8")
        except Exception:
            pass
    print("[run_pipeline] " + pretty, flush=True)
    # 3) stdin detached so NOTHING can block on input(); the env carries the tools.
    return subprocess.call(cmd, env=env, stdin=subprocess.DEVNULL)


if __name__ == "__main__":
    raise SystemExit(main())
