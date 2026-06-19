"""
pipeline/matrix.py — Genome presence/absence matrix.

Builds binary presence/absence matrices from hit tables and produces
clustered Plotly heatmaps for visualising gene distribution across genomes.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import pandas as pd

try:
    import plotly.graph_objects as go
    _PLOTLY_AVAILABLE = True
except ImportError:
    _PLOTLY_AVAILABLE = False
    go = None  # type: ignore


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_matrix(
    hits_df: pd.DataFrame,
    confidence_tiers: Optional[list] = None,
) -> pd.DataFrame:
    """Build a binary genome × gene presence/absence matrix.

    Parameters
    ----------
    hits_df : pd.DataFrame
        Must contain ``genome_id``, ``hit_name`` (gene column),
        and optionally ``confidence_tier``.
    confidence_tiers : list, optional
        If provided, only rows with matching confidence_tier are included.
        E.g. ``["high_confidence", "putative"]``.

    Returns
    -------
    pd.DataFrame
        Binary matrix with genomes as rows, gene/hit names as columns.
        Values are 0 or 1.  Empty DataFrame if inputs are insufficient.
    """
    if hits_df.empty:
        return pd.DataFrame()

    # Accept flexible column names: genome_id / target_name, hit_name / protein_id
    df_in = hits_df.copy()
    if "genome_id" not in df_in.columns:
        for alt in ("target_name", "accession", "protein_id"):
            if alt in df_in.columns:
                df_in = df_in.rename(columns={alt: "genome_id"})
                break
    if "hit_name" not in df_in.columns:
        # For a single-gene family, every hit represents the same profile.
        if "target_name" in df_in.columns and df_in.get("genome_id", pd.Series()).equals(df_in.get("target_name", pd.Series())):
            df_in["hit_name"] = "protein_family"
        else:
            for alt in ("protein_id", "accession", "target_name"):
                if alt in df_in.columns:
                    df_in["hit_name"] = df_in[alt]
                    break
            else:
                df_in["hit_name"] = "gene"

    required = {"genome_id", "hit_name"}
    if not required.issubset(df_in.columns):
        missing = required - set(df_in.columns)
        print(f"WARNING: matrix.build_matrix — missing columns: {missing}", file=sys.stderr)
        return pd.DataFrame()

    hits_df = df_in

    df = hits_df.copy()

    # Filter by confidence tier if requested
    if confidence_tiers and "confidence_tier" in df.columns:
        df = df[df["confidence_tier"].isin(confidence_tiers)]

    if df.empty:
        return pd.DataFrame()

    # Add a temporary presence indicator
    df = df[["genome_id", "hit_name"]].drop_duplicates()
    df["present"] = 1

    try:
        matrix = df.pivot_table(
            index="genome_id",
            columns="hit_name",
            values="present",
            aggfunc="max",
            fill_value=0,
        ).astype(int)
        matrix.columns.name = None
        matrix.index.name   = "genome_id"
    except Exception as exc:
        print(f"ERROR: Could not build presence/absence matrix: {exc}", file=sys.stderr)
        return pd.DataFrame()

    return matrix


def matrix_stats(matrix_df: pd.DataFrame) -> dict:
    """Compute summary statistics for a presence/absence matrix.

    Parameters
    ----------
    matrix_df : pd.DataFrame
        Output of :func:`build_matrix`.

    Returns
    -------
    dict
        {n_genomes, n_genes, avg_genes_per_genome, avg_genomes_per_gene,
         core_genes, accessory_genes}
        core_genes: gene names present in >90% of genomes.
        accessory_genes: gene names present in <10% of genomes.
    """
    empty = {
        "n_genomes":            0,
        "n_genes":              0,
        "avg_genes_per_genome": 0.0,
        "avg_genomes_per_gene": 0.0,
        "core_genes":           [],
        "accessory_genes":      [],
    }

    if matrix_df.empty:
        return empty

    n_genomes = len(matrix_df)
    n_genes   = len(matrix_df.columns)

    if n_genomes == 0 or n_genes == 0:
        return empty

    genes_per_genome = matrix_df.sum(axis=1)
    genomes_per_gene = matrix_df.sum(axis=0)

    prevalence = genomes_per_gene / n_genomes

    core      = list(prevalence[prevalence >= 0.90].index)
    accessory = list(prevalence[prevalence <  0.10].index)

    return {
        "n_genomes":            n_genomes,
        "n_genes":              n_genes,
        "avg_genes_per_genome": round(float(genes_per_genome.mean()), 2),
        "avg_genomes_per_gene": round(float(genomes_per_gene.mean()), 2),
        "core_genes":           sorted(core),
        "accessory_genes":      sorted(accessory),
    }


def heatmap_figure(
    matrix_df: pd.DataFrame,
    max_genomes: int = 100,
) -> "go.Figure":  # type: ignore[name-defined]
    """Produce a clustered Plotly heatmap of the presence/absence matrix.

    Rows (genomes) are ordered by Jaccard-distance hierarchical clustering
    when scipy is available; otherwise sorted alphabetically.

    Parameters
    ----------
    matrix_df : pd.DataFrame
        Binary matrix from :func:`build_matrix`.
    max_genomes : int
        Cap the number of displayed genomes for legibility.

    Returns
    -------
    go.Figure
        Plotly Figure object. Returns an empty Figure on failure.
    """
    if not _PLOTLY_AVAILABLE:
        print("ERROR: plotly not installed.", file=sys.stderr)
        return go.Figure() if go else None  # type: ignore[union-attr]

    if matrix_df.empty:
        return go.Figure()

    # Limit rows
    df = matrix_df.head(max_genomes).copy()

    # Attempt clustering
    df = _cluster_matrix(df)

    z      = df.values.tolist()
    y_labs = list(df.index)
    x_labs = list(df.columns)

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=x_labs,
            y=y_labs,
            colorscale=[
                [0.0, "#f0f0f0"],
                [1.0, "#2563eb"],
            ],
            showscale=False,
            hoverongaps=False,
            hovertemplate="Genome: %{y}<br>Gene: %{x}<br>Present: %{z}<extra></extra>",
        )
    )

    n_shown = len(y_labs)
    fig.update_layout(
        title=dict(
            text=f"Presence/Absence Matrix ({n_shown} genomes × {len(x_labs)} genes)",
            font=dict(size=14),
        ),
        xaxis=dict(
            tickangle=-45,
            tickfont=dict(size=max(6, min(11, 300 // max(len(x_labs), 1)))),
            title="Gene",
        ),
        yaxis=dict(
            tickfont=dict(size=max(6, min(10, 300 // max(n_shown, 1)))),
            title="Genome",
            autorange="reversed",
        ),
        margin=dict(l=160, r=20, t=60, b=120),
        height=max(400, min(1200, 12 * n_shown + 150)),
    )

    return fig


def heatmap_png(
    matrix_df: pd.DataFrame,
    out_path: Optional[Path] = None,
    dpi: int = 300,
) -> bytes:
    """Render the presence/absence matrix as a publication PNG (bytes).

    Falls back to a matplotlib static figure so no plotly/orca dependency
    is needed for image export.
    """
    import io
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    matplotlib.rcParams.update({
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "font.family":  "sans-serif",
    })

    if matrix_df.empty:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return buf.getvalue()

    df = matrix_df.head(100).copy()
    df = _cluster_matrix(df)

    n_g, n_genes = df.shape
    fig_h = max(4, min(20, 0.2 * n_g + 2))
    fig_w = max(6, min(20, 0.4 * n_genes + 2))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(df.values, aspect="auto", cmap="Blues", interpolation="nearest",
              vmin=0, vmax=1)

    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(df.columns, rotation=45, ha="right",
                       fontsize=max(5, min(8, 80 // max(n_genes, 1))))
    ax.set_yticks(range(n_g))
    ax.set_yticklabels(df.index, fontsize=max(5, min(8, 80 // max(n_g, 1))))
    ax.set_xlabel("Gene", fontsize=9)
    ax.set_ylabel("Genome", fontsize=9)
    ax.set_title(f"Presence/Absence ({n_g} genomes × {n_genes} genes)", fontsize=11)

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    data = buf.getvalue()

    if out_path:
        _op = Path(out_path)
        _op.parent.mkdir(parents=True, exist_ok=True)
        _op.write_bytes(data)

    return data


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _cluster_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Re-order rows by hierarchical clustering (Jaccard distance) if scipy available."""
    if len(df) < 3 or len(df.columns) < 2:
        return df.sort_index()

    try:
        from scipy.spatial.distance import pdist
        from scipy.cluster.hierarchy import linkage, leaves_list

        dist = pdist(df.values, metric="jaccard")
        # Guard against all-zero or all-same rows
        if any(d != d for d in dist):  # NaN check
            return df.sort_index()
        Z = linkage(dist, method="average")
        order = leaves_list(Z)
        return df.iloc[order]
    except ImportError:
        # scipy not available — fall back to alphabetical sort
        return df.sort_index()
    except Exception as exc:
        print(f"WARNING: Clustering failed, using alphabetical order: {exc}", file=sys.stderr)
        return df.sort_index()
