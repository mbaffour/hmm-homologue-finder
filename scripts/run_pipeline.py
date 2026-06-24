#!/usr/bin/env python3
"""run_pipeline.py — run the WHOLE discovery pipeline from Python, fully autonomously,
with ZERO interactive input.

It is a thin, no-prompt launcher around hmm_finder.py that:
  1. puts the running interpreter's own bin directory on PATH, so the conda env's tools
     (hmmsearch, mafft, prodigal, iqtree, …) are found WITHOUT `conda activate`;
  2. runs with stdin detached (DEVNULL) so no step can ever block on a prompt;
  3. injects no-prompt defaults (full database set; offline unless you pass --email)
     and forwards every other flag straight to hmm_finder.py.

USAGE — invoke it with the conda env's python (so the env's tools are on PATH):
    ~/miniforge3/envs/hmm-discovery/bin/python scripts/run_pipeline.py \
        --fasta /path/to/seed.faa --find-interrupted --out-dir ~/hmm_runs/myrun --iterations 2 --cpu 4

  • Only required flag: --fasta (the seed protein/nucleotide FASTA).
  • Organism names from NCBI: add `--email you@inst.edu` (otherwise it runs fully offline).
  • Any hmm_finder.py flag works (--databases, --name, --no-controls, --color-by, …); run
    `python scripts/hmm_finder.py --help` for the full list. This launcher only adds defaults
    that you have NOT already specified.

The process exits 0 on success; the run writes its package, tables, figures and report under
--out-dir (default: a folder beside the seed). Nothing is ever asked of you.
"""
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> int:
    args = list(sys.argv[1:])

    # 1) make the conda env's binaries discoverable (this interpreter lives in <env>/bin,
    #    alongside hmmsearch/mafft/…). No `conda activate` needed.
    env = dict(os.environ)
    env["PATH"] = os.path.dirname(sys.executable) + os.pathsep + env.get("PATH", "")

    # 2) no-prompt defaults (only if the user did not already choose):
    def _has(*names):
        return any(a == n or a.startswith(n + "=") for a in args for n in names)

    # full default database set (avoids the "choose databases" prompt on a TTY)
    if not _has("--databases", "--all-databases", "--pick-databases", "--list-databases", "--smoke"):
        args.append("--all-databases")
    # offline organism annotation unless an email is given (avoids the email prompt; never sends
    # a placeholder address)
    if not _has("--email", "--no-annotate"):
        args.append("--no-annotate")
    if not _has("--skip-tool-check"):
        args.append("--skip-tool-check")

    cmd = [sys.executable, str(HERE / "hmm_finder.py")] + args
    print("[run_pipeline] " + " ".join(cmd), flush=True)
    # 3) stdin detached so NOTHING can block on input(); the env carries the tools.
    return subprocess.call(cmd, env=env, stdin=subprocess.DEVNULL)


if __name__ == "__main__":
    raise SystemExit(main())
