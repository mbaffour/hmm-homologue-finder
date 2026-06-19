#!/usr/bin/env python3
"""
build_tree_of_hits.py
=====================
Build a maximum-likelihood phylogenetic tree of the DISCOVERED homologues
(distinct from the seed-only tree the pipeline makes).

Steps: MAFFT align -> trimAl (-gt 0.5) -> IQ-TREE (ModelFinder + 1000 UFBoot).

INPUT  : a FASTA of unique, ORF-validated family domain proteins
         (e.g. runA/benchmark/validated/hits_unique_aa.faa)
OUTPUT : <out-dir>/hits.aln.faa, hits.aln.trim.faa, hits.treefile, hits.iqtree,
         and a PNG/SVG rendering if toytree is available.

USAGE
-----
  python3 build_tree_of_hits.py --faa <unique_aa.faa> --out-dir <dir> [--cpu 8]
"""
from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

# Ensure conda env tools (mafft, trimal, iqtree) are on PATH regardless of caller.
_cand = ([Path(os.environ["CONDA_PREFIX"]) / "bin"] if os.environ.get("CONDA_PREFIX") else []) \
    + [Path.home() / _n / "envs" / "hmm-discovery" / "bin"
       for _n in ("miniforge3", "mambaforge", "miniconda3", "anaconda3")]
for _b in _cand:
    if _b.is_dir():
        os.environ["PATH"] = f"{_b}{os.pathsep}{os.environ.get('PATH', '')}"
        break


def run(cmd: list[str], **kw) -> None:
    print("  $", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--faa", type=Path, required=True, help="unique family domain AA FASTA")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--cpu", default="8")
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    aln = out / "hits.aln.faa"
    trim = out / "hits.aln.trim.faa"
    prefix = out / "hits"

    n = sum(1 for ln in args.faa.read_text().splitlines() if ln.startswith(">"))
    print(f"Building tree from {n} unique family domains: {args.faa}")
    if n < 4:
        print("Fewer than 4 sequences; IQ-TREE needs >=4. Aborting tree build.")
        return

    # 1. Align
    with aln.open("w") as fh:
        run(["mafft", "--auto", "--thread", str(args.cpu), str(args.faa)], stdout=fh)
    # 2. Trim gappy columns (the -gt 0.5 that keeps the alignment compact)
    run(["trimal", "-in", str(aln), "-out", str(trim), "-gt", "0.5"])
    # 3. ML tree with model selection + ultrafast bootstrap
    run(["iqtree", "-s", str(trim), "-m", "MFP", "-B", "1000",
         "-T", str(args.cpu), "--prefix", str(prefix), "-redo"])

    # 4. Optional rendering
    try:
        import toytree  # noqa: F401
        import toyplot.png, toyplot.svg
        tre = toytree.tree(str(prefix.with_suffix(".treefile")))
        canvas, _, _ = tre.draw(width=900, height=max(400, 14 * n), tip_labels_align=True)
        toyplot.png.render(canvas, str(out / "hits_tree.png"))
        toyplot.svg.render(canvas, str(out / "hits_tree.svg"))
        print(f"  rendered hits_tree.png / .svg")
    except Exception as e:
        print(f"  (rendering skipped: {e}; tree file is at {prefix}.treefile)")

    print(f"Done. Tree: {prefix}.treefile")


if __name__ == "__main__":
    main()
