"""
pipeline/taxonomy.py — Taxonomy parsing + distribution data.

Extracts taxonomy information from sequence IDs and description strings,
and produces data structures suitable for Plotly Sankey and treemap figures.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# NCBI accession prefix → broad host-type heuristic
# ---------------------------------------------------------------------------
_ACCESSION_HOST_HINTS: dict[str, str] = {
    "NP_":  "bacteria",
    "WP_":  "bacteria",
    "YP_":  "bacteria",
    "AP_":  "bacteria",
    "NC_":  "bacteria",
    "AE":   "bacteria",
    "CP":   "bacteria",
    "NZ_":  "bacteria",
}

_PHAGE_KEYWORDS = re.compile(
    r"phage|virus|viru|bacteriophage|siphovir|podo|myovir|inovir|microvir",
    re.IGNORECASE,
)

_EUKARYOTE_HOSTS = re.compile(
    r"human|homo sapiens|mus musculus|murine|bovine|chicken|plant|arabidopsis",
    re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def extract_taxonomy_from_id(seq_id: str) -> dict:
    """Parse taxonomy hints from an NCBI-style or INPHARED accession.

    For INPHARED IDs, pipe-delimited taxonomy strings are parsed directly.
    For standard NCBI accessions, host type is inferred from prefix patterns.

    Parameters
    ----------
    seq_id : str
        Sequence identifier (may include description after first whitespace).

    Returns
    -------
    dict
        {accession, organism_name, taxonomy_string, host_type}
    """
    if not seq_id or not isinstance(seq_id, str):
        return _empty_tax()

    # Split off description
    parts = seq_id.strip().split(None, 1)
    accession   = parts[0]
    description = parts[1] if len(parts) > 1 else ""

    # ---- INPHARED format: ACC|taxonomy_string|... ----
    if "|" in accession:
        pipe_parts = accession.split("|")
        acc = pipe_parts[0]
        # Look for taxonomy in any pipe-delimited field
        tax_str = ""
        organism = ""
        for p in pipe_parts[1:]:
            if ";" in p:
                tax_str = p
            elif p and not re.match(r"^\d+$", p):
                organism = organism or p
        host_type = _infer_host_type(tax_str + " " + organism + " " + description)
        return {
            "accession":       acc,
            "organism_name":   organism or _extract_organism(description),
            "taxonomy_string": tax_str,
            "host_type":       host_type,
        }

    # ---- 6-frame translated ORF: CONTIG|frame_F1|100_200 ----
    if "|frame_" in accession:
        contig = accession.split("|frame_")[0]
        return {
            "accession":       contig,
            "organism_name":   "",
            "taxonomy_string": "",
            "host_type":       "unknown",
        }

    # ---- Standard NCBI accession ----
    host_type = "unknown"
    for prefix, hint in _ACCESSION_HOST_HINTS.items():
        if accession.startswith(prefix):
            host_type = hint
            break

    organism = _extract_organism(description)
    if host_type == "unknown":
        host_type = _infer_host_type(description + " " + organism)

    return {
        "accession":       accession,
        "organism_name":   organism,
        "taxonomy_string": "",
        "host_type":       host_type,
    }


def taxonomy_table(hits_df: pd.DataFrame) -> pd.DataFrame:
    """Add taxonomy columns to the hits DataFrame.

    Parameters
    ----------
    hits_df : pd.DataFrame
        Must contain ``protein_id`` column; optionally ``description``.

    Returns
    -------
    pd.DataFrame
        hits_df with added columns:
        accession, organism_name, taxonomy_string, host_type.
    """
    # Accept target_name as alias for protein_id (output of parse_tblout)
    df = hits_df.copy()
    if "protein_id" not in df.columns and "target_name" in df.columns:
        df["protein_id"] = df["target_name"]

    if df.empty or "protein_id" not in df.columns:
        return df

    # Combine protein_id with description for richer parsing
    def _combined(row: pd.Series) -> str:
        pid  = str(row.get("protein_id", ""))
        desc = str(row.get("description", ""))
        return f"{pid} {desc}".strip()

    tax_rows = df.apply(lambda row: extract_taxonomy_from_id(_combined(row)), axis=1)
    tax_df = pd.DataFrame(list(tax_rows))

    for col in ["accession", "organism_name", "taxonomy_string", "host_type"]:
        df[col] = tax_df[col].values

    return df


def sankey_data(taxonomy_df: pd.DataFrame) -> dict:
    """Produce a Plotly Sankey-compatible data structure.

    The flow is: host_type → phage_family → confidence_tier.

    Parameters
    ----------
    taxonomy_df : pd.DataFrame
        Output of :func:`taxonomy_table`; needs host_type, confidence_tier,
        and optionally phage_family.

    Returns
    -------
    dict
        Plotly Sankey ``link`` + ``node`` dicts:
        {node_labels, node_colors, link_source, link_target, link_value}
    """
    if taxonomy_df.empty:
        return {
            "node_labels": [],
            "node_colors": [],
            "link_source": [],
            "link_target": [],
            "link_value":  [],
        }

    df = taxonomy_df.copy()

    # Ensure required columns exist with defaults
    if "host_type" not in df.columns:
        df["host_type"] = "unknown"
    if "phage_family" not in df.columns:
        df["phage_family"] = df.apply(_infer_phage_family, axis=1)
    if "confidence_tier" not in df.columns:
        df["confidence_tier"] = "unknown"

    # Build node list (unique labels in flow order)
    host_types   = sorted(df["host_type"].fillna("unknown").unique())
    families     = sorted(df["phage_family"].fillna("unknown").unique())
    conf_tiers   = sorted(df["confidence_tier"].fillna("unknown").unique())

    all_nodes = host_types + families + conf_tiers
    node_idx  = {label: i for i, label in enumerate(all_nodes)}

    # Assign colours by tier
    _TIER_COLOURS = {
        "high_confidence": "#2ca02c",
        "putative":        "#1f77b4",
        "divergent":       "#ff7f0e",
        "likely_fp":       "#d62728",
        "unknown":         "#aec7e8",
    }
    _HOST_COLOURS = {
        "bacteria": "#8c564b",
        "unknown":  "#c7c7c7",
    }

    node_colors = []
    for label in all_nodes:
        if label in _TIER_COLOURS:
            node_colors.append(_TIER_COLOURS[label])
        elif label in _HOST_COLOURS:
            node_colors.append(_HOST_COLOURS[label])
        else:
            node_colors.append("#aec7e8")

    # Build links
    link_source: list[int] = []
    link_target: list[int] = []
    link_value:  list[int] = []

    # host_type → phage_family
    for (h, f), grp in df.groupby(
        [df["host_type"].fillna("unknown"), df["phage_family"].fillna("unknown")]
    ):
        if h in node_idx and f in node_idx:
            link_source.append(node_idx[h])
            link_target.append(node_idx[f])
            link_value.append(len(grp))

    # phage_family → confidence_tier
    for (f, c), grp in df.groupby(
        [df["phage_family"].fillna("unknown"), df["confidence_tier"].fillna("unknown")]
    ):
        if f in node_idx and c in node_idx:
            link_source.append(node_idx[f])
            link_target.append(node_idx[c])
            link_value.append(len(grp))

    return {
        "node_labels": all_nodes,
        "node_colors": node_colors,
        "link_source": link_source,
        "link_target": link_target,
        "link_value":  link_value,
    }


def treemap_data(
    taxonomy_df: pd.DataFrame,
    groupby: str = "host_type",
) -> pd.DataFrame:
    """Produce a Plotly treemap-ready DataFrame.

    Parameters
    ----------
    taxonomy_df : pd.DataFrame
        Output of :func:`taxonomy_table`.
    groupby : str
        Primary grouping column (e.g. ``"host_type"``).

    Returns
    -------
    pd.DataFrame
        Columns: labels, parents, values — suitable for ``go.Treemap``.
    """
    if taxonomy_df.empty:
        return pd.DataFrame(columns=["labels", "parents", "values"])

    df = taxonomy_df.copy()

    if "phage_family" not in df.columns:
        df["phage_family"] = df.apply(_infer_phage_family, axis=1)
    if "host_type" not in df.columns:
        df["host_type"] = "unknown"

    group_col  = groupby if groupby in df.columns else ("host_type" if "host_type" in df.columns else df.columns[0])
    conf_col   = "confidence_tier" if "confidence_tier" in df.columns else None

    df[group_col]  = df[group_col].fillna("unknown")
    df["phage_family"] = df["phage_family"].fillna("unknown")

    rows: list[dict] = [{"labels": "All Hits", "parents": "", "values": len(df)}]

    for grp_val, grp_df in df.groupby(group_col):
        rows.append({
            "labels":  str(grp_val),
            "parents": "All Hits",
            "values":  len(grp_df),
        })
        for fam, fam_df in grp_df.groupby("phage_family"):
            label = f"{grp_val} / {fam}"
            rows.append({
                "labels":  label,
                "parents": str(grp_val),
                "values":  len(fam_df),
            })
            if conf_col:
                for conf, conf_df in fam_df.groupby(conf_col):
                    rows.append({
                        "labels":  f"{label} / {conf}",
                        "parents": label,
                        "values":  len(conf_df),
                    })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _empty_tax() -> dict:
    return {
        "accession":       "",
        "organism_name":   "",
        "taxonomy_string": "",
        "host_type":       "unknown",
    }


def _extract_organism(description: str) -> str:
    """Try to extract '[Organism name]' from NCBI-style descriptions."""
    match = re.search(r"\[([^\[\]]+)\]", description)
    return match.group(1) if match else ""


def _infer_host_type(text: str) -> str:
    """Infer whether a hit comes from a phage, bacterium, or eukaryote."""
    if _PHAGE_KEYWORDS.search(text):
        return "phage"
    if _EUKARYOTE_HOSTS.search(text):
        return "eukaryote"
    if re.search(r"bacteri|archae|prokaryot", text, re.IGNORECASE):
        return "bacteria"
    return "unknown"


def _infer_phage_family(row: pd.Series) -> str:
    """Infer phage family from taxonomy_string or organism_name."""
    combined = " ".join([
        str(row.get("taxonomy_string", "")),
        str(row.get("organism_name",   "")),
        str(row.get("description",     "")),
    ])
    family_patterns = {
        "Siphoviridae":  r"siphovi|siphovir",
        "Myoviridae":    r"myovir",
        "Podoviridae":   r"podovir",
        "Inoviridae":    r"inovir",
        "Microviridae":  r"microvir",
        "Herelleviridae":r"herellevi",
        "Autographiviridae": r"autographivi",
        "Drexlerviridae":r"drexleri",
        "Demerecviridae":r"demerecvi",
    }
    for family, pattern in family_patterns.items():
        if re.search(pattern, combined, re.IGNORECASE):
            return family
    return "unclassified"
