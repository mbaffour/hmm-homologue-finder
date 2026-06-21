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


def _safe(s: str, n: int = 55) -> str:
    """Newick/MAFFT-safe label: keep alnum . _ - ; everything else -> '_'."""
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(s)).strip("_")[:n] or "x"


def _short_acc(rec_id: str) -> str:
    """Extract a compact accession from a record id: 'sp|P49861|NAME' -> 'P49861',
    otherwise the id unchanged."""
    parts = str(rec_id).split("|")
    return parts[1] if len(parts) >= 2 and parts[0] in ("sp", "tr") else str(rec_id)


def _organism_from_desc(desc: str) -> str:
    """Best-effort organism name from a FASTA description.
    Handles UniProt ('... OS=Escherichia phage T5 OX=...'), NCBI bracket
    ('... [Escherichia phage phiKT]'), and NCBI genome titles
    ('NC_019520.1:.. Escherichia phage phiKT, complete genome'). '' if unknown."""
    d = str(desc or "")
    m = re.search(r"OS=(.+?)(?:\s+(?:OX|GN|PE|SV)=|$)", d)
    if m:
        return m.group(1).strip()
    m = re.search(r"\[([^\]]+)\]", d)
    if m:
        return m.group(1).strip()
    # strip leading accession / accession:coords tokens, then trailing genome boilerplate
    d = re.sub(r"^(?:[A-Za-z]{1,5}[_0-9][\w.]*(?::\d+-\d+)?\s+)+", "", d)
    d = re.sub(r",?\s*(complete genome|complete sequence|genomic sequence|complete cds|"
               r"genome assembly.*|DNA, complete.*)\s*$", "", d, flags=re.I)
    return d.strip().strip(",").strip()


def _canonical_organism(name: str, fallback: str = "") -> str:
    """Canonical key for COUNTING unique organisms. Collapses host-genus aliases of
    the same phage — 'Enterobacteria phage N4' and 'Escherichia phage N4' both ->
    'n4' — by keying on the name after 'phage'/'virus'. Unnamed/metagenomic entries
    ('uncultured virus', blank) are NOT collapsed: they fall back to the genome
    accession so each distinct genome still counts once. Returns '' if nothing."""
    s = re.sub(r"^(UNVERIFIED:?|MAG:?|TPA(?:_asm)?:?)\s*", "", str(name or "").strip(), flags=re.I).strip()
    if (not s) or re.search(r"uncultured|unclassified|metagenom|environmental", s, re.I):
        return (str(fallback).strip() or s).lower()
    m = re.search(r"\b(?:phage|virus)\b\s+(.+)$", s, flags=re.I)
    key = m.group(1).strip() if (m and m.group(1).strip()) else s
    return re.sub(r"\s+", " ", key).strip().lower()


def _uniquify(label: str, used: set) -> str:
    """Ensure a tip label is unique (Newick requires it) by appending _2, _3 …."""
    base, k, out = label, 2, label
    while out in used:
        out = f"{base}_{k}"
        k += 1
    used.add(out)
    return out


def _build_tree_input(hits_faa: Path, seeds_faa, hits_tsv, out_path: Path) -> int:
    """Write the alignment/tree input FASTA with ORGANISM-FIRST tip labels so every
    figure shows readable names (accession only as a fallback). Discovered hits are
    labelled 'Organism_accession_xN' where N is how many distinct ORGANISMS carry
    that EXACT domain sequence (so a single tip honestly conveys how widespread the
    gene is — the tree is deduplicated by sequence, but the occurrence count is
    not lost). Counting unique organisms (by name) rather than genome accessions
    avoids double-counting the same phage that appears in several databases under
    different accessions. Seeds are labelled 'Organism_accession_seed'."""
    from Bio import SeqIO
    hit_map = _organism_labels(hits_tsv) if (hits_tsv and Path(hits_tsv).exists()) else {}
    # How many distinct ORGANISMS carry each exact domain sequence — recovers the
    # discovery breadth that exact-sequence dedup would otherwise hide on the tree.
    occ: dict = {}
    if hits_tsv and Path(hits_tsv).exists():
        per_seq: dict = {}
        try:
            with Path(hits_tsv).open(newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    s = row.get("aa_sequence", "")
                    org = _canonical_organism(row.get("organism", ""),
                                              row.get("genome_id", ""))
                    if s and org:
                        per_seq.setdefault(s, set()).add(org)
            occ = {s: len(g) for s, g in per_seq.items()}
        except Exception:
            occ = {}
    used: set = set()
    recs = []
    for rec in SeqIO.parse(str(hits_faa), "fasta"):
        base = hit_map.get(rec.id) or _safe(rec.id)         # organism_genomeid, else id
        n = occ.get(str(rec.seq), 1)
        label = f"{base}_x{n}" if n > 1 else base           # _xN = carried by N organisms
        rec.id = _uniquify(label, used); rec.name = rec.id; rec.description = ""
        recs.append(rec)
    if seeds_faa and Path(seeds_faa).exists():
        for rec in SeqIO.parse(str(seeds_faa), "fasta"):
            acc = _short_acc(rec.id)
            org = _organism_from_desc(rec.description)
            label = f"{_safe(org)}_{_safe(acc)}_seed" if org else f"{_safe(acc)}_seed"
            rec.id = _uniquify(label, used); rec.name = rec.id; rec.description = ""
            recs.append(rec)
    SeqIO.write(recs, str(out_path), "fasta")
    return len(recs)


def _organism_labels(hits_tsv) -> dict:
    """Map hit_id -> 'Organism_accession' (newick-safe) from a hits.tsv.
    Accepts a str or Path."""
    mp = {}
    try:
        with Path(hits_tsv).open(newline="") as fh:
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


def _render_newick(newick_path: Path, out_dir: Path, stem: str, layout: str = "r") -> list:
    """Render a publication tree to {stem}.svg/.png/.pdf — editable vector (Inkscape
    / Illustrator) + 300-dpi raster, on an OPAQUE WHITE background (toytree's default
    canvas is transparent → black on dark viewers). Tips are coloured by host genus
    (first token of the organism label) with a legend, and well-supported nodes
    (UFBoot >= 80) carry a dot. layout 'r' = rectangular, 'c' = circular. Graceful
    if toytree/toyplot/matplotlib are absent."""
    try:
        import math
        import toytree  # noqa: F401
        import toyplot, toyplot.svg, toyplot.png, toyplot.marker
        import matplotlib
        import matplotlib.colors as mcolors
    except Exception as e:
        print(f"  ({stem} rendering skipped: {e}; Newick at {newick_path})")
        return []
    try:
        tre = toytree.tree(str(newick_path))
        tips = [str(t) for t in tre.get_tip_labels()]
        ntips = len(tips)
        # colour tips by host genus (first underscore token of the organism label)
        genera = [t.split("_")[0] for t in tips]
        uniq = sorted(set(genera))
        cmap = matplotlib.colormaps["tab20"].resampled(max(len(uniq), 1))
        gcol = {g: mcolors.to_hex(cmap(i)) for i, g in enumerate(uniq)}
        tip_colors = [gcol[g] for g in genera]
        # support dots (UFBoot >= 80) on internal nodes
        sizes, ncolors = [], []
        for v in list(tre.get_node_data("support")):
            ok = isinstance(v, (int, float)) and not math.isnan(v) and v >= 80
            sizes.append(7 if ok else 0)
            ncolors.append("#222222" if ok else "transparent")
        if layout == "c":
            dim = max(700, 9 * ntips)
            canvas = toyplot.Canvas(width=dim + 280, height=dim, style={"background-color": "white"})
        else:
            canvas = toyplot.Canvas(width=1300, height=max(400, 18 * ntips),
                                    style={"background-color": "white"})
        ax = canvas.cartesian(padding=25)
        ax.show = False
        tre.draw(axes=ax, layout=layout, tip_labels_align=(layout != "c"),
                 tip_labels_colors=tip_colors, node_sizes=sizes, node_colors=ncolors)
        # host-genus colour legend + a support-dot key
        try:
            entries = [(g, toyplot.marker.create(shape="o", size=11, mstyle={"fill": gcol[g]}))
                       for g in uniq]
            entries.append(("UFBoot >= 80", toyplot.marker.create(shape="o", size=9,
                                                                  mstyle={"fill": "#222222"})))
            canvas.legend(entries, corner=("right", 6, 150, 16 * len(entries)))
        except Exception as e:
            print(f"  ({stem} legend skipped: {e})")
        made = []
        toyplot.svg.render(canvas, str(out_dir / f"{stem}.svg")); made.append("svg")
        try:
            toyplot.png.render(canvas, str(out_dir / f"{stem}.png")); made.append("png")
        except Exception as e:
            print(f"  ({stem} PNG skipped: {e})")
        try:
            import toyplot.pdf
            toyplot.pdf.render(canvas, str(out_dir / f"{stem}.pdf")); made.append("pdf")
        except Exception as e:
            print(f"  ({stem} PDF skipped: {e})")
        print(f"  rendered {stem}.{{{','.join(made)}}} ({ntips} tips, {len(uniq)} host genera)")
        return made
    except Exception as e:
        print(f"  ({stem} render skipped: {e})")
        return []


def _homologs_only_newick(treefile: Path, out_path: Path):
    """Write a copy of the tree with the seed tips (labels ending '_seed') pruned,
    so the discovered homologs can be shown in a legible 'result' figure separate
    from the seed-context tree. Returns out_path, or None if not applicable."""
    try:
        from Bio import Phylo
        t = Phylo.read(str(treefile), "newick")
        seeds = [tp for tp in t.get_terminals() if (tp.name or "").endswith("_seed")]
        if not seeds or (len(t.get_terminals()) - len(seeds)) < 3:
            return None
        for tp in seeds:
            t.prune(tp)
        Phylo.write(t, str(out_path), "newick")
        return out_path
    except Exception as e:
        print(f"  (homologs-only prune skipped: {e})")
        return None


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

    # Build the alignment/tree input with ORGANISM-FIRST tip labels (hits from
    # hits.tsv; optional seeds parsed from their headers, marked '*_seed'). Doing
    # this up front means BOTH the alignment figure and the tree show readable
    # organism names, with accession only as a fallback.
    input_faa = out / "tree_input.faa"
    n = _build_tree_input(args.faa, args.seeds, args.hits_tsv, input_faa)
    extra = " (+ seeds, marked *_seed)" if (args.seeds and args.seeds.exists()) else ""
    print(f"Aligning {n} organism-labelled sequences{extra}: {input_faa}")
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

    # Tip labels are already organism-first (set when the input was built), so the
    # Newick tree, the alignment, and the figure all carry readable names directly.
    treefile = Path(str(prefix) + ".treefile")
    render_tree = treefile

    # 4. Render the full tree (homologs + any seeds, organism-labelled), then ALSO
    # a legible homologs-only tree (seed tips pruned) — large trees are dense once
    # the seeds are added, so the result figure is easier to read on its own.
    _render_newick(render_tree, out, "hits_tree")
    ho = _homologs_only_newick(treefile, out / "hits_tree_homologs_only.treefile")
    if ho:
        _render_newick(ho, out, "hits_tree_homologs_only")
        # (circular layout is left to iTOL/FigTree on the exported Newick — toytree
        # 3.0.10 doesn't implement a circular layout.)

    print(f"Done. Tree: {prefix}.treefile")


if __name__ == "__main__":
    main()
