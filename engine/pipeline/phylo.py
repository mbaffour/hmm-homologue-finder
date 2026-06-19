"""
pipeline/phylo.py — IQ-TREE phylogenetic tree + toytree figure rendering.

Runs IQ-TREE for maximum-likelihood tree inference and renders
tip-labelled, colour-coded trees with toytree.  Falls back gracefully
when optional dependencies (toytree / toyplot) are absent.
"""
from __future__ import annotations

import io
import json
import shutil
import subprocess
import sys
from pathlib import Path

from .utils import find_tool, run_cmd
from typing import Optional

import pandas as pd


# ---------------------------------------------------------------------------
# Confidence-tier colour map (mirrors confidence.py tiers)
# ---------------------------------------------------------------------------

_TIER_COLORS: dict[str, str] = {
    "high_confidence": "#4CAF50",   # green
    "putative":        "#2196F3",   # blue
    "divergent":       "#FF9800",   # orange
    "likely_fp":       "#F44336",   # red
    "seed":            "#9C27B0",   # purple (seed sequences)
    "unknown":         "#9E9E9E",   # gray
}


def _iqtree_binary() -> Optional[str]:
    """Return the full path to the first IQ-TREE binary found, or None."""
    for name in ("iqtree2", "iqtree"):
        path = find_tool(name)
        if path:
            return path
    return None


# ---------------------------------------------------------------------------
# IQ-TREE runner
# ---------------------------------------------------------------------------

def run_iqtree(
    aln_path: Path,
    out_dir: Path,
    model: str = "TEST",
    bootstrap: int = 1000,
    cpu: int = 4,
) -> dict:
    """Run IQ-TREE on an alignment file.

    Parameters
    ----------
    aln_path : Path
        Input alignment in FASTA or PHYLIP format.
    out_dir : Path
        Output directory; IQ-TREE files are prefixed ``iqtree``.
    model : str
        Substitution model or ``"TEST"`` for automatic model selection.
    bootstrap : int
        Number of ultra-fast bootstrap replicates (UFBoot, -B flag).
    cpu : int
        Thread count.

    Returns
    -------
    dict
        {treefile: Path | None, logfile: Path | None,
         success: bool, model_used: str, error: str}
    """
    aln_path = Path(aln_path)
    out_dir  = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    result: dict = {
        "treefile":   None,
        "logfile":    None,
        "success":    False,
        "model_used": model,
        "error":      "",
    }

    binary = _iqtree_binary()
    if binary is None:
        result["error"] = "iqtree2/iqtree not found on PATH"
        print(f"WARNING: {result['error']}", file=sys.stderr)
        return result

    if not aln_path.exists():
        result["error"] = f"Alignment file not found: {aln_path}"
        print(f"ERROR: {result['error']}", file=sys.stderr)
        return result

    # ── Sanitise sequence IDs ─────────────────────────────────────────
    # IQ-TREE fails on duplicated or truncated IDs. Pre-process the
    # alignment to guarantee unique IDs ≤ 50 chars.
    try:
        from Bio import SeqIO, AlignIO
        import re as _re

        aln = list(SeqIO.parse(str(aln_path), "fasta"))
        seen: dict[str, int] = {}
        for rec in aln:
            # Sanitise: keep alphanum + _ + - only, max 40 chars
            safe = _re.sub(r"[^A-Za-z0-9._\-]", "_", rec.id)[:40]
            count = seen.get(safe, 0)
            seen[safe] = count + 1
            rec.id = safe if count == 0 else f"{safe}_{count}"
            rec.description = ""

        sanitised_aln = out_dir / "iqtree_input.faa"
        SeqIO.write(aln, str(sanitised_aln), "fasta")
        aln_path = sanitised_aln
    except Exception as _exc:
        print(f"WARNING: ID sanitisation failed ({_exc}), using original file", file=sys.stderr)

    prefix = str(out_dir / "iqtree")
    cmd = [
        binary,
        "-s", str(aln_path),
        "-m", model,
        "-T", str(cpu),
        "--prefix", prefix,
        "-redo",
    ]
    if bootstrap and bootstrap > 0:
        cmd.extend(["-B", str(max(int(bootstrap), 1000))])

    print(f"INFO: Running: {' '.join(cmd)}", file=sys.stderr)
    try:
        proc = run_cmd(cmd, timeout=7200)
    except subprocess.TimeoutExpired:
        result["error"] = "IQ-TREE timed out after 2 hours"
        print(f"ERROR: {result['error']}", file=sys.stderr)
        return result
    except Exception as exc:
        result["error"] = f"IQ-TREE subprocess error: {exc}"
        print(f"ERROR: {result['error']}", file=sys.stderr)
        return result

    logfile   = out_dir / "iqtree.log"
    treefile  = out_dir / "iqtree.treefile"

    # IQ-TREE writes log to stderr
    logfile.write_text(proc.stdout + "\n" + proc.stderr)

    result["logfile"] = logfile

    if proc.returncode != 0:
        result["error"] = f"IQ-TREE exited with code {proc.returncode}"
        print(f"ERROR: {result['error']}\n{proc.stderr[-2000:]}", file=sys.stderr)
        return result

    if not treefile.exists():
        result["error"] = "IQ-TREE finished but treefile not found"
        print(f"ERROR: {result['error']}", file=sys.stderr)
        return result

    result["treefile"] = treefile
    result["success"]  = True

    # Parse best-fit model from IQ-TREE log
    for line in proc.stdout.splitlines():
        if "Best-fit model:" in line:
            try:
                result["model_used"] = line.split("Best-fit model:")[1].strip().split()[0]
            except IndexError:
                pass
            break

    return result


# ---------------------------------------------------------------------------
# Tree rendering
# ---------------------------------------------------------------------------

def _render_tree_body(treefile, hits_df, out_dir, result):
    """Inner body of render_tree — split out for exception safety.

    Renders the IQ-TREE Newick tree with matplotlib via Bio.Phylo. This is far
    more robust and controllable than coordinate-fragile interactive libraries:
    the canvas scales to the tip count, tip labels are coloured by confidence
    tier, and a proper legend lists only the tiers actually present. Both PNG
    (raster, for the report) and SVG (vector, for manuscripts) are written.
    """
    # Build protein_id → tier lookup so tip labels can be coloured by tier.
    tier_map: dict[str, str] = {}
    if hits_df is not None and not hits_df.empty:
        for _, row in hits_df.iterrows():
            pid  = str(row.get("protein_id", "") or "")
            tier = str(row.get("confidence_tier", "unknown") or "unknown")
            if pid:
                tier_map[pid] = tier

    try:
        import matplotlib
        matplotlib.use("Agg")          # headless raster backend
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from Bio import Phylo
    except ImportError as exc:
        msg = f"matplotlib/Bio.Phylo not installed ({exc}); writing placeholder files"
        print(f"WARNING: {msg}", file=sys.stderr)
        result["error"] = msg
        png_path = out_dir / "tree.png"
        svg_path = out_dir / "tree.svg"
        png_path.write_bytes(b"")
        svg_path.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg">'
            "<text>matplotlib/Bio.Phylo not installed</text></svg>"
        )
        result["png_path"] = png_path
        result["svg_path"] = svg_path
        return result

    try:
        tree = Phylo.read(str(treefile), "newick")
    except Exception as exc:
        result["error"] = f"Failed to parse treefile: {exc}"
        print(f"ERROR: {result['error']}", file=sys.stderr)
        return result

    # Ladderise for a tidy, conventional ordering.
    try:
        tree.ladderize()
    except Exception:
        pass

    # Count tips and size the canvas so labels never overlap.
    terminals = tree.get_terminals()
    n_tips = max(1, len(terminals))
    fig_w = 10.0
    fig_h = max(3.5, n_tips * 0.32 + 1.2)   # ~0.32 inch per tip

    def _tier_for(name: str) -> str:
        return tier_map.get(str(name or ""), "unknown")

    # Bio.Phylo colours labels via a {clade_name: color} mapping.
    label_colors = {
        t.name: _TIER_COLORS.get(_tier_for(t.name), _TIER_COLORS["unknown"])
        for t in terminals if t.name
    }

    # Show tip names only (suppress internal-node support labels as text — they
    # would clutter a small tree); colour each tip by its tier.
    def _label(clade):
        return clade.name if clade.is_terminal() and clade.name else ""

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 8,
        "svg.fonttype": "none",   # editable text in the SVG
    })
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    try:
        Phylo.draw(
            tree,
            axes=ax,
            do_show=False,
            label_func=_label,
            label_colors=label_colors,
            show_confidence=False,
        )
    except Exception as exc:
        # Last-resort minimal draw (no per-tip colours).
        try:
            ax.clear()
            Phylo.draw(tree, axes=ax, do_show=False, label_func=_label,
                       show_confidence=False)
        except Exception as exc2:
            plt.close(fig)
            result["error"] = f"Bio.Phylo draw failed: {exc2}"
            print(f"ERROR: {result['error']}", file=sys.stderr)
            return result

    ax.set_title(f"Maximum-likelihood tree (n={n_tips} sequences)",
                 fontsize=10, fontweight="bold")
    ax.set_xlabel("Substitutions per site", fontsize=8)
    ax.set_ylabel("")
    # Bio.Phylo prints a y tick per tip; hide the numeric ticks for a clean look.
    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    # Legend: only the tiers that actually appear among the tips.
    present = []
    seen = set()
    for t in terminals:
        tier = _tier_for(t.name)
        if tier not in seen:
            seen.add(tier)
            present.append(tier)
    order = ["high_confidence", "putative", "divergent", "likely_fp", "seed", "unknown"]
    present.sort(key=lambda x: order.index(x) if x in order else 99)
    handles = [
        mpatches.Patch(color=_TIER_COLORS.get(t, _TIER_COLORS["unknown"]),
                       label=t.replace("_", " ").title())
        for t in present
    ]
    if handles:
        ax.legend(handles=handles, loc="lower right", fontsize=7,
                  framealpha=0.9, edgecolor="#cccccc", title="Confidence tier",
                  title_fontsize=7)

    fig.tight_layout()

    png_path = out_dir / "tree.png"
    svg_path = out_dir / "tree.svg"
    try:
        fig.savefig(str(png_path), dpi=200, bbox_inches="tight")
        result["png_path"] = png_path
    except Exception as exc:
        print(f"WARNING: PNG export failed: {exc}", file=sys.stderr)
    try:
        fig.savefig(str(svg_path), bbox_inches="tight")
        result["svg_path"] = svg_path
    except Exception as exc:
        print(f"WARNING: SVG export failed: {exc}", file=sys.stderr)
    plt.close(fig)

    result["success"] = bool(result.get("png_path") or result.get("svg_path"))
    return result

# Internal alias

def render_tree(
    treefile: Path,
    hits_df: pd.DataFrame,
    out_dir: Path,
) -> dict:
    """Render an IQ-TREE treefile with toytree, coloured by confidence tier.

    Parameters
    ----------
    treefile : Path
        Newick treefile produced by IQ-TREE.
    hits_df : pd.DataFrame
        Hit table with columns ``protein_id`` and ``confidence_tier``.
        Seed sequences should have tier ``"seed"`` or will fall back to
        ``"unknown"``.
    out_dir : Path
        Directory where PNG and SVG are written.

    Returns
    -------
    dict
        {png_path: Path | None, svg_path: Path | None, success: bool, error: str}
    """
    treefile = Path(treefile)
    out_dir  = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    result: dict = {
        "png_path": None,
        "svg_path": None,
        "success":  False,
        "error":    "",
    }

    if not treefile.exists():
        result["error"] = f"Treefile not found: {treefile}"
        print(f"ERROR: {result['error']}", file=sys.stderr)
        return result

    # Wrap the whole rendering in a safety net so ghostscript /
    # toytree dependency issues never crash the caller.
    try:
        return _render_tree_body(treefile, hits_df, out_dir, result)
    except Exception as exc:
        result["error"] = f"render_tree failed: {exc}"
        print(f"WARNING: {result['error']}", file=sys.stderr)
        return result



# ---------------------------------------------------------------------------
# Newick → JSON (minimal, for optional frontend display)
# ---------------------------------------------------------------------------

def newick_to_dict(treefile: Path) -> dict:
    """Convert a Newick treefile to a minimal JSON-serialisable dict.

    Returns an empty dict if toytree is unavailable or parsing fails.
    The structure is: {name, branch_length, children: [...]} (recursive).

    Parameters
    ----------
    treefile : Path
        Path to a Newick treefile.

    Returns
    -------
    dict
        Tree as nested dict, or {} on failure.
    """
    treefile = Path(treefile)
    if not treefile.exists():
        return {}

    try:
        import toytree
        tt = toytree.tree(str(treefile))
        return _toytree_to_dict(tt.treenode)
    except Exception:
        pass

    # Fallback: try Bio.Phylo
    try:
        from Bio import Phylo
        from io import StringIO

        newick_str = treefile.read_text()
        tree = Phylo.read(StringIO(newick_str), "newick")

        def _clade_to_dict(clade) -> dict:
            return {
                "name":          clade.name or "",
                "branch_length": clade.branch_length or 0.0,
                "confidence":    clade.confidence,
                "children":      [_clade_to_dict(c) for c in clade.clades],
            }

        return _clade_to_dict(tree.root)
    except Exception:
        return {}


def _toytree_to_dict(node) -> dict:
    """Recursively convert a toytree TreeNode to a plain dict."""
    try:
        return {
            "name":          node.name or "",
            "branch_length": float(node.dist) if node.dist else 0.0,
            "support":       float(node.support) if node.support else None,
            "children":      [_toytree_to_dict(c) for c in node.children],
        }
    except Exception:
        return {}
