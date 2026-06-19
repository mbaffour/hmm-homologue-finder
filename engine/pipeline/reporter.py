"""
pipeline/reporter.py — HTML report (Jinja2) + reproducibility JSON + METHODS_TEXT.

Produces:
  - reports/reproducibility.json  (full audit record)
  - reports/METHODS_TEXT.txt      (manuscript-ready methods paragraph)
  - reports/summary_report.html   (rendered Jinja2 HTML summary)
  - results/export_YYYYMMDD.zip   (all output directories archived)
"""
from __future__ import annotations

import json
import sys
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd


def _compact_database_provenance(benchmark: dict) -> list[dict]:
    rows = []
    for name, rec in (benchmark.get("databases", {}) or {}).items():
        if not isinstance(rec, dict):
            continue
        files = rec.get("files") or {}
        file_records = files.values() if isinstance(files, dict) else files
        urls = []
        checksums = []
        accessed = []
        size_bytes = 0
        for info in file_records:
            if not isinstance(info, dict):
                continue
            url = info.get("url")
            checksum = info.get("sha256")
            access_time = info.get("accessed_at") or info.get("downloaded_at")
            if url and url not in urls:
                urls.append(url)
            if checksum and checksum not in checksums:
                checksums.append(checksum)
            if access_time:
                accessed.append(str(access_time))
            try:
                size_bytes += int(info.get("size_bytes") or 0)
            except Exception:
                pass
        rows.append(
            {
                "database": name,
                "status": rec.get("status", ""),
                "type": rec.get("type", ""),
                "search_mode": rec.get("search_mode", ""),
                "nt_orf_mode": rec.get("nt_orf_mode", ""),
                "file_count": rec.get("file_count", len(urls)),
                "source_url_count": len(urls),
                "source_urls": urls[:10],
                "source_size_bytes": size_bytes,
                "source_accessed_first": min(accessed) if accessed else "",
                "source_accessed_last": max(accessed) if accessed else "",
                "source_sha256_count": len(checksums),
                "source_sha256s": checksums[:10],
                "hit_count": rec.get("hit_count", ""),
                "strict_count": rec.get("strict_count", rec.get("strict_hit_count", "")),
                "runtime_seconds": rec.get("runtime_seconds", ""),
                "error": rec.get("error", ""),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Reproducibility JSON
# ---------------------------------------------------------------------------

def build_reproducibility_json(
    proj_dir: Path,
    hits_df: Optional[pd.DataFrame],
    state: dict,
    tools: dict,
) -> dict:
    """
    Assemble complete reproducibility record and save to reports/reproducibility.json.
    """
    proj_dir = Path(proj_dir)
    reports_dir = proj_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Audit trail
    trail_file = proj_dir / "logs" / "audit_trail.jsonl"
    commands = []
    if trail_file.exists():
        for line in trail_file.read_text().splitlines():
            line = line.strip()
            if line:
                try:
                    commands.append(json.loads(line))
                except Exception:
                    pass

    # Hit summary
    hit_summary: dict = {}
    if hits_df is not None and not hits_df.empty:
        if "confidence_tier" in hits_df.columns:
            hit_summary = hits_df["confidence_tier"].value_counts().to_dict()
        hit_summary["total"] = len(hits_df)
        if "genome_id" in hits_df.columns:
            hit_summary["n_genomes"] = hits_df["genome_id"].nunique()

    # Input summary from state
    input_params = state.get("input", {}).get("params", {}) if isinstance(state, dict) else {}
    benchmark = state.get("benchmark", {}) if isinstance(state, dict) else {}
    if not benchmark:
        manifest_file = proj_dir / "benchmark_manifest.json"
        if manifest_file.exists():
            try:
                benchmark = json.loads(manifest_file.read_text())
            except Exception:
                benchmark = {}

    citation_guidance = {
        "summary": (
            "Cite HMM Discovery plus every external tool and database selected "
            "or used in this run. Installed-but-unused tools do not need to be "
            "cited for a specific analysis."
        ),
        "core_tools": [
            "HMMER (hmmbuild, hmmsearch, hmmscan, hmmpress)",
            "MAFFT",
            "trimAl",
            "Clustal Omega when used",
            "Prodigal when conventional gene prediction or synteny rescue is used",
            "seqkit when large FASTA files are split",
            "IQ-TREE / ModelFinder / ultrafast bootstrap when phylogeny is run",
            "MEME Suite (MEME/FIMO) when motif analysis is run",
            "CD-HIT and MMseqs2 when clustering is run",
            "DIAMOND when reciprocal/background similarity checks are run",
            "clinker, pyGenomeViz, or Easyfig-compatible exports when synteny visualization is run",
            "toytree, toyplot, Ghostscript, matplotlib, Plotly, and Kaleido when figures are exported",
        ],
        "database_catalog": [
            "INPHARED genomes",
            "INPHARED vConTACT2 proteins",
            "UniProt Swiss-Prot",
            "NCBI RefSeq viral proteins",
            "NCBI RefSeq viral genomes",
            "NCBI RefSeq bacterial proteins",
            "Gut Phage Database (GPD)",
            "GVD-AVrC / Aggregated Gut Viral Catalogue",
            "Pfam-A sequences",
            "Pfam-A HMM/domain library",
            "VOGDB VFAM HMMs and annotations",
            "NCBI Entrez/GenBank/RefSeq records when remote synteny context is fetched",
        ],
        "release_files": ["ACKNOWLEDGEMENTS.md", "CITATION.cff"],
    }

    repro = {
        "generated_at":  datetime.now().isoformat(),
        "project_dir":   str(proj_dir),
        "tool_versions": {
            k: {"version": v.get("version"), "description": v.get("description")}
            for k, v in (tools or {}).items()
            if v.get("available")
        },
        "pipeline_steps": state if isinstance(state, dict) else {},
        "input_summary":  input_params,
        "hit_summary":    hit_summary,
        "database_provenance": _compact_database_provenance(benchmark)
        if isinstance(benchmark, dict)
        else [],
        "citation_guidance": citation_guidance,
        "all_commands":   commands,
        "n_commands":     len(commands),
    }

    out = reports_dir / "reproducibility.json"
    try:
        out.write_text(json.dumps(repro, indent=2, default=str))
    except Exception as exc:
        print(f"WARNING: Could not write reproducibility.json: {exc}", file=sys.stderr)
    return repro


# ---------------------------------------------------------------------------
# Methods text
# ---------------------------------------------------------------------------

def generate_methods_text(proj_dir: Path, repro: dict) -> str:
    """
    Generate a manuscript-ready Methods paragraph from repro dict.
    Saves to reports/METHODS_TEXT.txt and returns the string.
    """
    proj_dir = Path(proj_dir)
    tv = repro.get("tool_versions", {})

    def ver(tool, fallback="unknown version"):
        """Return a clean, short version string (number only when possible)."""
        import re as _re
        raw = tv.get(tool, {}).get("version") or fallback
        # Try to extract just the version number from strings like
        # "trimAl v1.5.rev1 build[...]" or "IQ-TREE version 3.1.1 for MacOS..."
        m = _re.search(r"v?(\d+\.\d+[\w.\-]*)", str(raw))
        if m:
            return m.group(1)
        return raw

    # Determine aligner used
    aligner         = "MAFFT"
    aligner_version = ver("mafft")
    if "clustalo" in tv:
        aligner         = "Clustal Omega"
        aligner_version = ver("clustalo")
    elif "muscle" in tv:
        aligner         = "MUSCLE"
        aligner_version = ver("muscle")

    # Determine IQ-TREE version label
    iqtree_ver = ver("iqtree2") if "iqtree2" in tv else ver("iqtree")

    steps = repro.get("pipeline_steps", {})
    inp   = repro.get("input_summary", {})
    hs    = repro.get("hit_summary", {})

    # Extract search parameters from state
    search_params  = steps.get("search", {}).get("params", {})
    if not search_params:
        # Also try top-level pipeline_steps for flat state representations
        search_params = repro.get("search_params", {})

    n_seeds        = inp.get("seq_count", inp.get("n_seeds", "N"))
    evalue         = search_params.get("evalue", "1e-5")
    strict         = search_params.get("strict", search_params.get("strict_threshold", 45))
    moderate       = search_params.get("moderate", search_params.get("moderate_threshold", 30))
    hmm_cov        = float(search_params.get("hmm_cov_floor", search_params.get("hmm_cov", 0.30)))
    profile_length = (
        inp.get("hmm_length", inp.get("hmm_leng", inp.get("profile_length", "N")))
    )
    trimal_ver     = ver("trimal")
    hmmer_ver      = ver("hmmbuild", ver("hmmsearch", "unknown"))

    # Databases searched
    db_list = search_params.get("databases", [])
    if isinstance(db_list, str):
        db_list = [db_list]
    dbs_str = ", ".join(db_list) if db_list else "the registered databases"

    n_hc      = hs.get("high_confidence", 0)
    n_put     = hs.get("putative", 0)
    n_genomes = hs.get("n_genomes", 0)

    # Format evalue for readability
    try:
        evalue_f = float(str(evalue).replace("×10", "e").replace("⁻", "-"))
        evalue_str = f"{evalue_f:.0e}".replace("e-0", "e-").replace("e+0", "e+")
    except (ValueError, TypeError):
        evalue_str = str(evalue)

    text = (
        f"Profile HMM searches were conducted using HMMER {hmmer_ver} (Eddy 2011). "
        f"Seed sequences (n = {n_seeds}) were aligned with {aligner} {aligner_version} "
        f"and poorly conserved alignment columns were removed with trimAl {trimal_ver} "
        f"(Capella-Gutierrez et al. 2009) using the 'automated1' heuristic. "
        f"The resulting profile HMM ({profile_length} match positions) was searched "
        f"against {dbs_str} using hmmsearch with an E-value threshold of {evalue_str}. "
        f"Hits were classified using a multi-evidence confidence engine combining bit "
        f"score (strict ≥ {strict} bits, moderate ≥ {moderate} bits), "
        f"HMM profile coverage (≥ {int(hmm_cov * 100)} %), reciprocal "
        f"hmmsearch validation against the seed set, and taxonomic consistency scoring. "
        f"{n_hc} high-confidence and {n_put} putative homologs were identified across "
        f"{n_genomes} genomes. Phylogenetic analysis was performed using IQ-TREE "
        f"{iqtree_ver} (Nguyen et al. 2015) with the best-fit substitution model "
        f"selected by ModelFinder and ultrafast bootstrap support (1,000 replicates). "
        f"Run-specific tool versions, selected databases, source URLs, access dates, "
        f"and checksums are reported in reproducibility.json; database and software "
        f"citations should include every tool and database selected for the run."
    )

    out = proj_dir / "reports" / "METHODS_TEXT.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        out.write_text(text)
    except Exception as exc:
        print(f"WARNING: Could not write METHODS_TEXT.txt: {exc}", file=sys.stderr)
    return text


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

def render_html_report(proj_dir: Path, context: dict) -> Path:
    """
    Render the Jinja2 HTML report template.
    Saves to reports/summary_report.html and returns the path.
    """
    proj_dir = Path(proj_dir)
    reports_dir = proj_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Find template relative to this file
    app_dir = Path(__file__).parent.parent
    template_path = app_dir / "www" / "templates" / "report.html.j2"

    if not template_path.exists():
        print(
            f"WARNING: Report template not found at {template_path}; "
            f"writing minimal fallback HTML.",
            file=sys.stderr,
        )
        html = _minimal_html_report(context)
    else:
        try:
            from jinja2 import Environment, FileSystemLoader, select_autoescape
            env = Environment(
                loader=FileSystemLoader(str(template_path.parent)),
                autoescape=select_autoescape(["html"]),
            )
            tmpl = env.get_template("report.html.j2")
            html = tmpl.render(**context)
        except ImportError:
            print("WARNING: jinja2 not installed; writing minimal HTML report.", file=sys.stderr)
            html = _minimal_html_report(context)
        except Exception as exc:
            print(f"WARNING: Jinja2 render failed: {exc}; writing minimal fallback.", file=sys.stderr)
            html = _minimal_html_report(context)

    out = reports_dir / "summary_report.html"
    try:
        out.write_text(html)
    except Exception as exc:
        print(f"WARNING: Could not write summary_report.html: {exc}", file=sys.stderr)
    return out


def _minimal_html_report(context: dict) -> str:
    """Return a bare-minimum HTML report when Jinja2 is unavailable."""
    hs      = context.get("hit_stats", {})
    name    = context.get("project_name", "HMM Project")
    date    = context.get("generated_at", datetime.now().strftime("%Y-%m-%d"))
    methods = context.get("methods_text", "")

    rows = ""
    for tier in ("high_confidence", "putative", "divergent", "likely_fp"):
        rows += f"<tr><td>{tier.replace('_', ' ').title()}</td><td>{hs.get(tier, 0)}</td></tr>\n"

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<title>{name} — HMM Discovery Report</title>
<style>body{{font-family:sans-serif;margin:2em;}}
table{{border-collapse:collapse;}}td,th{{border:1px solid #ccc;padding:.4em .8em;}}
pre{{background:#f5f5f5;padding:1em;white-space:pre-wrap;}}</style>
</head><body>
<h1>{name}</h1><p>Generated: {date}</p>
<h2>Hit Summary</h2>
<table><thead><tr><th>Tier</th><th>Count</th></tr></thead><tbody>
{rows}
<tr><td><b>Total</b></td><td>{hs.get('total', 0)}</td></tr>
<tr><td><b>Genomes</b></td><td>{hs.get('n_genomes', 0)}</td></tr>
</tbody></table>
<h2>Methods</h2><pre>{methods}</pre>
</body></html>"""


def build_report_context(
    proj_dir: Path,
    hits_df: Optional[pd.DataFrame],
    repro: dict,
    methods_text: str,
    tools: dict,
) -> dict:
    """Build the Jinja2 template context dict."""
    proj_dir = Path(proj_dir)
    hs = repro.get("hit_summary", {})

    # Databases searched (from audit trail)
    db_stats = []
    commands = repro.get("all_commands", [])
    for cmd_rec in commands:
        if cmd_rec.get("step") == "search":
            db_name = cmd_rec.get("extra", {}).get("db_name", "unknown")
            db_stats.append({
                "name": db_name,
                "type": cmd_rec.get("extra", {}).get("db_type", "protein"),
                "hits": cmd_rec.get("extra", {}).get("hit_count", 0),
                "high_confidence": cmd_rec.get("extra", {}).get("strict_count", 0),
            })

    # Key figures.
    # The app writes figures to figures/; the command-line pipeline writes them
    # to results/ (sometimes under slightly different names). Try every known
    # location/alias so the report embeds figures regardless of which pipeline
    # produced them. The first existing candidate wins; the stored path is
    # relative to the reports/ directory where the HTML is written.
    key_figures = []
    figure_candidates: list[tuple[str, list[str]]] = [
        ("HMM Profile Logo",        ["figures/hmm_logo.png", "results/hmm_logo.png"]),
        ("Score Distribution",      ["figures/score_calibration.png",
                                     "results/score_calibration.png"]),
        ("Phylogenetic Tree",       ["figures/tree.png", "results/tree.png"]),
        ("Presence/Absence Heatmap",["figures/heatmap.png",
                                     "results/heatmap.png",
                                     "results/presence_absence_heatmap.png"]),
        ("Synteny Map",             ["figures/synteny_map.png",
                                     "results/synteny_map.png"]),
        ("Taxonomic Distribution",  ["figures/taxonomy_sankey.png",
                                     "results/taxonomy_sankey.png"]),
    ]
    for title, candidates in figure_candidates:
        chosen = next((rel for rel in candidates if (proj_dir / rel).exists()), None)
        key_figures.append({
            "title": title,
            "path": f"../{chosen}" if chosen else None,
        })

    return {
        "project_name": proj_dir.name,
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M UTC"),
        "hit_stats": {
            "high_confidence": hs.get("high_confidence", 0),
            "putative": hs.get("putative", 0),
            "divergent": hs.get("divergent", 0),
            "likely_fp": hs.get("likely_fp", 0),
            "total": hs.get("total", 0),
            "n_genomes": hs.get("n_genomes", 0),
        },
        "databases_searched": db_stats,
        "tool_versions": tools,
        "methods_text": methods_text,
        "key_figures": key_figures,
    }


# ---------------------------------------------------------------------------
# Export ZIP
# ---------------------------------------------------------------------------

def create_export_zip(proj_dir: Path) -> Path:
    """
    Package results/, figures/, reports/, trees/ into a ZIP file.
    Returns path to the ZIP.
    """
    proj_dir = Path(proj_dir)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    zip_path = proj_dir / "results" / f"export_{ts}.zip"
    zip_path.parent.mkdir(parents=True, exist_ok=True)

    dirs_to_include = ["results", "figures", "reports", "trees", "logs"]

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for dirname in dirs_to_include:
            d = proj_dir / dirname
            if d.exists():
                for fpath in sorted(d.rglob("*")):
                    if fpath.is_file() and fpath != zip_path:
                        arcname = fpath.relative_to(proj_dir)
                        try:
                            zf.write(fpath, arcname)
                        except Exception:
                            pass

    return zip_path
