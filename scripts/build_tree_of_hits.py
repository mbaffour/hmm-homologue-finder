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
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

# Ensure conda env tools (mafft, trimal, iqtree) are on PATH regardless of caller.
from env_paths import ensure_env_on_path  # noqa: E402  (sibling helper in scripts/)
ensure_env_on_path()

# Reuse the engine's alignment quality + publication-figure helpers (and the
# accuracy-tiered MAFFT strategy) so the homolog alignment is a first-class,
# high-quality deliverable. Best-effort: degrade gracefully if the engine isn't
# importable (tree + raw alignment still produced).
_ENGINE = Path(__file__).resolve().parent.parent / "engine"
try:
    sys.path.insert(0, str(_ENGINE))
    from pipeline.alignment import accuracy_flags, alignment_quality, alignment_figure  # noqa: E402
except Exception:
    accuracy_flags = alignment_quality = alignment_figure = None


def _mafft_strategy(n: int, mode: str) -> "tuple[list, str]":
    """Return (mafft_flags, label) for the requested alignment mode."""
    presets = {
        "linsi": (["--localpair", "--maxiterate", "1000"], "L-INS-i"),
        "ginsi": (["--globalpair", "--maxiterate", "1000"], "G-INS-i"),
        "einsi": (["--genafpair", "--maxiterate", "1000"], "E-INS-i"),
        "auto":  (["--auto"], "auto"),
        "fftns": (["--retree", "2"], "FFT-NS-2"),
    }
    if mode in presets:
        return presets[mode]
    # "accurate" (default): L-INS-i where tractable, auto for very large sets.
    if accuracy_flags is not None:
        flags, lab = accuracy_flags(n)
        return (flags or ["--auto"], lab)
    return ((["--localpair", "--maxiterate", "1000"], "L-INS-i") if 2 <= n <= 500
            else (["--auto"], "auto"))


def run(cmd: list[str], **kw) -> None:
    print("  $", " ".join(str(c) for c in cmd), flush=True)
    # Capture stderr instead of inheriting it. MAFFT's wrapper writes to
    # /dev/stderr, which fails with "Permission denied" under WSL when fd 2 is an
    # inherited pipe; giving the child its own stderr handle avoids that. Captured
    # stderr is surfaced only if the command fails. (Harmless on macOS/Linux.)
    kw.setdefault("stderr", subprocess.PIPE)
    r = subprocess.run(cmd, **kw)
    if r.returncode != 0:
        if getattr(r, "stderr", None):
            try:
                print(r.stderr.decode(errors="replace"))
            except Exception:
                print(r.stderr)
        raise subprocess.CalledProcessError(r.returncode, cmd)


def _combine_with_seeds(hits_faa: Path, seeds_faa: Path, out_path: Path) -> int:
    """Write discovered hits + seeds to one FASTA so the seeds are placed *within*
    the homolog tree. Seed record IDs are prefixed ``SEED_`` so they are visibly
    marked in the alignment, the figure, and the Newick tip labels. Returns the
    combined (de-duplicated by ID) sequence count."""
    from Bio import SeqIO
    seen: set[str] = set()
    recs = []
    for rec in SeqIO.parse(str(hits_faa), "fasta"):
        if rec.id in seen:
            continue
        seen.add(rec.id)
        recs.append(rec)
    if seeds_faa and Path(seeds_faa).exists():
        for rec in SeqIO.parse(str(seeds_faa), "fasta"):
            rec.id = f"SEED_{rec.id}"
            rec.name = rec.id
            rec.description = ""
            if rec.id in seen:
                continue
            seen.add(rec.id)
            recs.append(rec)
    SeqIO.write(recs, str(out_path), "fasta")
    return len(recs)


def _organism_labels(hits_tsv: Path) -> dict:
    """Map hit_id -> 'Organism_accession' (newick-safe) from a hits.tsv."""
    mp = {}
    try:
        with hits_tsv.open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                hid = (row.get("hit_id") or "").strip()
                org = (row.get("organism") or "").strip()
                gid = (row.get("genome_id") or "").strip()
                lab = f"{org}_{gid}" if org else gid
                lab = "".join(c if c.isalnum() or c in "._-" else "_" for c in lab)[:60]
                if hid and lab:
                    mp[hid] = lab
    except Exception:
        pass
    return mp


def _relabel_newick(newick: str, mapping: dict) -> str:
    for hid in sorted(mapping, key=len, reverse=True):
        newick = re.sub(r"\b" + re.escape(hid) + r"\b", mapping[hid], newick)
    return newick


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--faa", type=Path, required=True, help="unique family domain AA FASTA")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--cpu", default="8")
    ap.add_argument("--hits-tsv", type=Path, default=None,
                    help="optional hits.tsv; relabels tree tips with organism names")
    ap.add_argument("--seeds", type=Path, default=None,
                    help="optional seed protein FASTA to include in the tree/alignment; "
                         "seed tips are marked 'SEED_*' so you can see where your starting "
                         "sequences fall among the discovered homologs")
    ap.add_argument("--mafft-mode",
                    choices=("accurate", "linsi", "ginsi", "einsi", "auto", "fftns"),
                    default="accurate",
                    help="MAFFT strategy for the homolog alignment. 'accurate' (default) "
                         "uses L-INS-i where tractable and falls back to --auto for very "
                         "large sets; or force a specific strategy.")
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    aln = out / "hits.aln.faa"
    trim = out / "hits.aln.trim.faa"
    prefix = out / "hits"

    # Optionally fold the seeds into the alignment/tree (marked 'SEED_*') so the
    # reader can see where the starting sequences fall among discovered homologs.
    if args.seeds and args.seeds.exists():
        input_faa = out / "tree_input.faa"
        n = _combine_with_seeds(args.faa, args.seeds, input_faa)
        print(f"Aligning {n} sequences (discovered homologs + seeds, marked SEED_*): {input_faa}")
    else:
        input_faa = args.faa
        n = sum(1 for ln in args.faa.read_text().splitlines() if ln.startswith(">"))
        print(f"Aligning {n} sequences: {args.faa}")
    if n < 2:
        print("Fewer than 2 sequences; nothing to align.")
        return

    # 1. High-quality alignment (accuracy-tiered MAFFT).
    flags, strat = _mafft_strategy(n, args.mafft_mode)
    print(f"  MAFFT strategy: {strat}")
    with aln.open("w") as fh:
        run(["mafft", *flags, "--thread", str(args.cpu), str(input_faa)], stdout=fh)

    # 1b. Alignment quality stats + publication-ready coloured MSA figure — the
    # alignment is a first-class deliverable, not just tree input.
    if alignment_quality is not None:
        try:
            q = alignment_quality(aln)
            (out / "hits.aln.stats.json").write_text(json.dumps(q, indent=2))
            print(f"  alignment: {q.get('n_sequences')} seqs x {q.get('aln_length')} cols; "
                  f"{q.get('conserved_columns')} conserved cols; "
                  f"mean pairwise id {q.get('avg_pairwise_id')}%")
        except Exception as e:
            print(f"  (alignment stats skipped: {e})")
    if alignment_figure is not None:
        for fmt in ("png", "svg", "pdf"):
            try:
                alignment_figure(aln, out, fmt=fmt)
            except Exception as e:
                print(f"  (alignment figure {fmt} skipped: {e})")

    if n < 4:
        print("Fewer than 4 sequences; alignment written, skipping IQ-TREE (needs >=4).")
        return

    # 2. Trim gappy columns (the -gt 0.5 that keeps the alignment compact)
    run(["trimal", "-in", str(aln), "-out", str(trim), "-gt", "0.5"])
    # 3. ML tree with model selection + ultrafast bootstrap.
    # -T AUTO lets IQ-TREE pick an optimal thread count <= physical cores; -ntmax
    # caps it at the requested cpu. (A fixed -T greater than the core count makes
    # IQ-TREE abort with "more threads than CPU cores available".)
    # -seed fixes the stochastic ML search + UFBoot resampling so reruns on the
    # same alignment yield an identical tree (needed for golden-output regression).
    run(["iqtree", "-s", str(trim), "-m", "MFP", "-B", "1000",
         "-T", "AUTO", "-ntmax", str(args.cpu), "-seed", "12345",
         "--prefix", str(prefix), "-redo"])

    # 3b. Organism-labelled tree (readable figure): rewrite cryptic ORF tip IDs
    # to 'Organism_accession'. Keep the original treefile untouched.
    treefile = Path(str(prefix) + ".treefile")
    render_tree = treefile
    if args.hits_tsv and args.hits_tsv.exists() and treefile.exists():
        mp = _organism_labels(args.hits_tsv)
        if mp:
            labeled = out / "hits.labeled.treefile"
            labeled.write_text(_relabel_newick(treefile.read_text(), mp))
            render_tree = labeled
            print(f"  wrote organism-labelled tree: {labeled}")

    # 4. Optional rendering (uses the labelled tree when available). Emit editable
    # vector formats (SVG for Inkscape, PDF for Illustrator) plus a 300-dpi PNG
    # preview. toyplot SVG keeps text as real text elements (editable, not paths).
    try:
        import toytree  # noqa: F401
        import toyplot.png, toyplot.svg
        tre = toytree.tree(str(render_tree))
        canvas, _, _ = tre.draw(width=1000, height=max(400, 16 * n), tip_labels_align=True)
        made = []
        toyplot.svg.render(canvas, str(out / "hits_tree.svg")); made.append("svg")
        try:
            toyplot.png.render(canvas, str(out / "hits_tree.png")); made.append("png")
        except Exception as e:
            print(f"  (tree PNG skipped: {e})")
        try:
            import toyplot.pdf
            toyplot.pdf.render(canvas, str(out / "hits_tree.pdf")); made.append("pdf")
        except Exception as e:
            print(f"  (tree PDF skipped: {e})")
        print(f"  rendered hits_tree.{{{','.join(made)}}}")
    except Exception as e:
        print(f"  (tree rendering skipped: {e}; Newick tree is at {prefix}.treefile)")

    print(f"Done. Tree: {prefix}.treefile")


if __name__ == "__main__":
    main()
