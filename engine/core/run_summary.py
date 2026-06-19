"""
core/run_summary.py - Project-level run summary generation.

Builds a compact, reproducible summary from files already written by the app or
benchmark runner. The summary is intentionally factual: it reports inputs,
settings, database statuses, hit counts, confidence tiers, and output files
without adding biological interpretation.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any


def _read_json(path: Path) -> dict[str, Any]:
    try:
        if path.exists():
            return json.loads(path.read_text())
    except Exception as exc:
        return {"_error": str(exc), "_path": str(path)}
    return {}


def _read_tsv(path: Path):
    try:
        if path.exists():
            import pandas as pd

            return pd.read_csv(path, sep="\t")
    except Exception:
        pass
    return None


def _first_existing(root: Path, rels: list[str]) -> Path | None:
    for rel in rels:
        path = root / rel
        if path.exists():
            return path
    return None


def _value_counts(df, column: str, limit: int | None = None) -> dict[str, int]:
    if df is None or df.empty or column not in df.columns:
        return {}
    counts = df[column].fillna("unknown").astype(str).value_counts()
    if limit:
        counts = counts.head(limit)
    return {str(k): int(v) for k, v in counts.to_dict().items()}


def _numeric_summary(df, column: str) -> dict[str, float]:
    if df is None or df.empty or column not in df.columns:
        return {}
    try:
        series = df[column].dropna().astype(float)
    except Exception:
        return {}
    if series.empty:
        return {}
    return {
        "min": float(series.min()),
        "median": float(series.median()),
        "max": float(series.max()),
    }


def _format_bytes(value: Any) -> str:
    try:
        n = float(value or 0)
    except Exception:
        return "0 B"
    units = ["B", "KB", "MB", "GB", "TB"]
    idx = 0
    while n >= 1024 and idx < len(units) - 1:
        n /= 1024
        idx += 1
    if idx == 0:
        return f"{int(n)} {units[idx]}"
    return f"{n:.1f} {units[idx]}"


def _database_urls_from_manifest_record(rec: dict[str, Any]) -> list[str]:
    urls: list[str] = []
    for info in _database_file_records(rec):
        url = info.get("url") if isinstance(info, dict) else ""
        if url and url not in urls:
            urls.append(str(url))
    return urls


def _database_file_records(rec: dict[str, Any]) -> list[dict[str, Any]]:
    files = rec.get("files") or {}
    if isinstance(files, dict):
        return [info for info in files.values() if isinstance(info, dict)]
    if isinstance(files, list):
        return [info for info in files if isinstance(info, dict)]
    return []


def _source_access_window(files: list[dict[str, Any]]) -> dict[str, str]:
    dates = sorted(
        str(info.get("accessed_at") or info.get("downloaded_at") or "")
        for info in files
        if info.get("accessed_at") or info.get("downloaded_at")
    )
    if not dates:
        return {}
    return {"first": dates[0], "last": dates[-1]}


def _source_size_bytes(files: list[dict[str, Any]]) -> int:
    total = 0
    for info in files:
        try:
            total += int(info.get("size_bytes") or 0)
        except Exception:
            pass
    return total


def _source_sha256s(files: list[dict[str, Any]], limit: int = 5) -> list[str]:
    hashes = []
    for info in files:
        value = str(info.get("sha256") or "")
        if value and value not in hashes:
            hashes.append(value)
    return hashes[:limit]


def _top_hits(df, limit: int = 10) -> list[dict[str, Any]]:
    if df is None or df.empty:
        return []
    sort_col = "bit_score" if "bit_score" in df.columns else None
    work = df.copy()
    if sort_col:
        try:
            work[sort_col] = work[sort_col].astype(float)
            work = work.sort_values(sort_col, ascending=False)
        except Exception:
            pass
    cols = [
        "target_name",
        "query_name",
        "protein_id",
        "genome_id",
        "database_source",
        "db_type",
        "evalue",
        "bit_score",
        "confidence_tier",
        "description",
    ]
    keep = [c for c in cols if c in work.columns]
    rows = []
    for rec in work.head(limit)[keep].fillna("").to_dict(orient="records"):
        rows.append({str(k): _json_safe(v) for k, v in rec.items()})
    return rows


def _json_safe(value: Any) -> Any:
    try:
        import math

        if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
            return None
    except Exception:
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return value


def _existing_outputs(root: Path) -> list[dict[str, Any]]:
    candidates = [
        "results/hits_main.tsv",
        "results/hits_best_per_genome.tsv",
        "results/per_database_metrics.tsv",
        "results/all_database_summary.tsv",
        "results/presence_absence_matrix.tsv",
        "results/taxonomy_table.tsv",
        "results/synteny_table.tsv",
        "results/synteny_placement_report.tsv",
        "results/synteny_neighborhoods.gff3",
        "results/hits_proteins.faa",
        "reports/summary_report.html",
        "reports/METHODS_TEXT.txt",
        "reports/reproducibility.json",
        "reports/all_database_benchmark_report.html",
        "figures/synteny_map.png",
        "figures/synteny_map.svg",
        "figures/tree.png",
        "figures/heatmap.png",
        "benchmark_manifest.json",
    ]
    rows = []
    for rel in candidates:
        path = root / rel
        if path.exists():
            rows.append(
                {
                    "path": rel,
                    "bytes": path.stat().st_size,
                    "modified": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
                }
            )
    return rows


def summarize_project(project_dir: str | Path) -> dict[str, Any]:
    root = Path(project_dir).expanduser().resolve()
    reports = root / "reports"
    results = root / "results"
    repro = _read_json(reports / "reproducibility.json")
    manifest = _read_json(root / "benchmark_manifest.json")
    state = _read_json(root / ".pipeline_state.json")

    hits_path = _first_existing(root, ["results/hits_main.tsv", "hits_main.tsv", "output/hits_main.tsv"])
    hits = _read_tsv(hits_path) if hits_path else None
    best_path = _first_existing(root, ["results/hits_best_per_genome.tsv", "hits_best_per_genome.tsv"])
    best = _read_tsv(best_path) if best_path else None
    db_metrics_path = _first_existing(root, ["results/per_database_metrics.tsv", "results/all_database_summary.tsv"])
    db_metrics = _read_tsv(db_metrics_path) if db_metrics_path else None
    synteny_path = _first_existing(root, ["results/synteny_table.tsv", "synteny_table.tsv"])
    synteny = _read_tsv(synteny_path) if synteny_path else None
    placement_path = _first_existing(root, ["results/synteny_placement_report.tsv"])
    placement = _read_tsv(placement_path) if placement_path else None

    active_command = manifest.get("active_command", {}) if isinstance(manifest, dict) else {}
    input_summary = repro.get("input_summary", {}) if isinstance(repro, dict) else {}
    if not input_summary:
        input_summary = manifest.get("core", {}).get("input_summary", {}) if isinstance(manifest, dict) else {}
    if not input_summary:
        input_summary = state.get("input_summary", {}) if isinstance(state, dict) else {}

    hit_summary = {
        "total_hits": int(len(hits)) if hits is not None else 0,
        "best_per_genome_hits": int(len(best)) if best is not None else 0,
        "confidence_tiers": _value_counts(hits, "confidence_tier"),
        "databases": _value_counts(hits, "database_source", limit=30),
        "db_types": _value_counts(hits, "db_type"),
        "bit_score": _numeric_summary(hits, "bit_score"),
        "evalue": _numeric_summary(hits, "evalue"),
        "top_hits": _top_hits(hits),
    }

    database_rows = []
    if isinstance(manifest, dict) and manifest.get("databases"):
        for name, rec in (manifest.get("databases", {}) or {}).items():
            files = _database_file_records(rec)
            urls = _database_urls_from_manifest_record(rec)
            access_window = _source_access_window(files)
            source_size = _source_size_bytes(files)
            sha256s = _source_sha256s(files)
            database_rows.append(
                {
                    "database": name,
                    "status": rec.get("status", ""),
                    "type": rec.get("type", ""),
                    "search_mode": rec.get("search_mode", ""),
                    "file_count": rec.get("file_count", len(urls) if urls else ""),
                    "downloaded_bytes": rec.get("downloaded_bytes", ""),
                    "downloaded_human": _format_bytes(rec.get("downloaded_bytes", 0)),
                    "source_url_count": len(urls),
                    "source_urls": urls[:5],
                    "source_size_bytes": source_size,
                    "source_size_human": _format_bytes(source_size),
                    "source_accessed_first": access_window.get("first", ""),
                    "source_accessed_last": access_window.get("last", ""),
                    "source_sha256_count": len(
                        {
                            str(info.get("sha256"))
                            for info in files
                            if info.get("sha256")
                        }
                    ),
                    "source_sha256s": sha256s,
                    "nt_orf_mode": rec.get("nt_orf_mode", ""),
                    "hit_count": rec.get("hit_count", ""),
                    "strict_count": rec.get("strict_count", rec.get("strict_hit_count", "")),
                    "runtime_seconds": rec.get("runtime_seconds", ""),
                    "error": str(rec.get("error", ""))[:240],
                }
            )
    elif db_metrics is not None and not db_metrics.empty:
        for rec in db_metrics.fillna("").to_dict(orient="records"):
            database_rows.append(
                {
                    "database": str(rec.get("database", "")),
                    "status": str(rec.get("status", "")),
                    "type": str(rec.get("type", rec.get("db_type", ""))),
                    "search_mode": str(rec.get("search_mode", "")),
                    "file_count": _json_safe(rec.get("file_count", "")),
                    "downloaded_bytes": _json_safe(rec.get("downloaded_bytes", "")),
                    "downloaded_human": _format_bytes(rec.get("downloaded_bytes", 0)),
                    "source_url_count": "",
                    "source_urls": [],
                    "source_size_bytes": "",
                    "source_size_human": "",
                    "source_accessed_first": "",
                    "source_accessed_last": "",
                    "source_sha256_count": "",
                    "source_sha256s": [],
                    "nt_orf_mode": str(rec.get("nt_orf_mode", "")),
                    "hit_count": _json_safe(rec.get("hit_count", "")),
                    "strict_count": _json_safe(rec.get("strict_count", rec.get("strict_hit_count", ""))),
                    "runtime_seconds": _json_safe(rec.get("runtime_seconds", "")),
                    "error": str(rec.get("error", ""))[:240],
                }
            )

    synteny_summary = {
        "neighborhood_rows": int(len(synteny)) if synteny is not None else 0,
        "placement_rows": int(len(placement)) if placement is not None else 0,
        "placement_statuses": _value_counts(placement, "status") or _value_counts(placement, "placement_status"),
    }

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "project_dir": str(root),
        "project_name": root.name,
        "run_status": manifest.get("status", "unknown") if isinstance(manifest, dict) and manifest else "project",
        "active_command": active_command,
        "input_summary": input_summary,
        "hit_summary": hit_summary,
        "database_summary": database_rows,
        "synteny_summary": synteny_summary,
        "outputs": _existing_outputs(root),
        "notes": [
            "Nucleotide six-frame mode is the discovery mode for short, overlapping, noncanonical, or annotation-missed genes.",
            "Prodigal mode is a faster conventional annotation baseline and should not be treated as exhaustive ORF discovery.",
        ],
    }
    return summary


def render_summary_markdown(summary: dict[str, Any]) -> str:
    hs = summary.get("hit_summary", {})
    lines = [
        f"# Run Summary: {summary.get('project_name', 'HMM Discovery Project')}",
        "",
        f"- Generated: {summary.get('generated_at', '')}",
        f"- Project folder: `{summary.get('project_dir', '')}`",
        f"- Run status: `{summary.get('run_status', 'unknown')}`",
        "",
        "## Input",
    ]
    inp = summary.get("input_summary", {}) or {}
    if inp:
        for key, value in inp.items():
            lines.append(f"- {key}: {value}")
    else:
        lines.append("- No input summary file was found yet.")

    active = summary.get("active_command", {}) or {}
    lines.extend(["", "## Run Settings"])
    if active:
        for key in ["preset", "databases", "nt_orf_mode", "min_orf_aa", "cpu", "evalue", "min_free_gb", "fasta", "out"]:
            if key in active:
                lines.append(f"- {key}: `{active.get(key)}`")
    else:
        lines.append("- No benchmark command metadata was found.")

    lines.extend(
        [
            "",
            "## Hits",
            f"- Total hits: {hs.get('total_hits', 0)}",
            f"- Best-per-genome hits: {hs.get('best_per_genome_hits', 0)}",
        ]
    )
    tiers = hs.get("confidence_tiers", {}) or {}
    if tiers:
        lines.append("- Confidence tiers: " + ", ".join(f"{k}={v}" for k, v in tiers.items()))
    dbs = hs.get("databases", {}) or {}
    if dbs:
        lines.append("- Hits by database: " + ", ".join(f"{k}={v}" for k, v in dbs.items()))
    bits = hs.get("bit_score", {}) or {}
    if bits:
        lines.append(
            f"- Bit score range: min={bits.get('min'):.2f}, median={bits.get('median'):.2f}, max={bits.get('max'):.2f}"
        )

    top_hits = hs.get("top_hits", []) or []
    lines.extend(["", "## Top Hits"])
    if top_hits:
        for idx, rec in enumerate(top_hits, start=1):
            label = rec.get("target_name") or rec.get("protein_id") or rec.get("query_name") or "hit"
            db = rec.get("database_source", "")
            bit = rec.get("bit_score", "")
            ev = rec.get("evalue", "")
            tier = rec.get("confidence_tier", "")
            lines.append(f"{idx}. {label} | db={db} | bit={bit} | evalue={ev} | tier={tier}")
    else:
        lines.append("- No hits table was found yet.")

    lines.extend(["", "## Databases"])
    db_rows = summary.get("database_summary", []) or []
    if db_rows:
        for rec in db_rows:
            details = []
            if rec.get("type"):
                details.append(f"type={rec.get('type')}")
            if rec.get("search_mode"):
                details.append(f"search_mode={rec.get('search_mode')}")
            if rec.get("file_count") not in ("", None):
                details.append(f"files={rec.get('file_count')}")
            if rec.get("downloaded_bytes") not in ("", None):
                details.append(f"downloaded={rec.get('downloaded_human')}")
            if rec.get("source_size_bytes"):
                details.append(f"source_size={rec.get('source_size_human')}")
            if rec.get("source_accessed_first"):
                if rec.get("source_accessed_first") == rec.get("source_accessed_last"):
                    details.append(f"accessed={rec.get('source_accessed_first')}")
                else:
                    details.append(
                        f"accessed={rec.get('source_accessed_first')}..{rec.get('source_accessed_last')}"
                    )
            mode = f", nt_orf_mode={rec.get('nt_orf_mode')}" if rec.get("nt_orf_mode") else ""
            err = f", error={rec.get('error')}" if rec.get("error") else ""
            lines.append(
                f"- {rec.get('database')}: {rec.get('status')} "
                f"(hits={rec.get('hit_count')}, strict={rec.get('strict_count')}"
                f"{', ' + ', '.join(details) if details else ''}{mode}{err})"
            )
            urls = rec.get("source_urls") or []
            if urls:
                shown = "; ".join(urls[:3])
                suffix = f"; ... +{len(urls) - 3} more" if len(urls) > 3 else ""
                lines.append(f"  - Source URL(s): {shown}{suffix}")
            sha256s = rec.get("source_sha256s") or []
            if sha256s:
                shown = ", ".join(value[:16] for value in sha256s[:3])
                more = rec.get("source_sha256_count", len(sha256s))
                suffix = f", ... {more} total" if more and int(more) > len(sha256s[:3]) else ""
                lines.append(f"  - SHA256 prefix(es): {shown}{suffix}")
    else:
        lines.append("- No database summary was found yet.")

    syn = summary.get("synteny_summary", {}) or {}
    lines.extend(
        [
            "",
            "## Synteny",
            f"- Neighborhood rows: {syn.get('neighborhood_rows', 0)}",
            f"- Placement rows: {syn.get('placement_rows', 0)}",
        ]
    )
    statuses = syn.get("placement_statuses", {}) or {}
    if statuses:
        lines.append("- Placement statuses: " + ", ".join(f"{k}={v}" for k, v in statuses.items()))

    lines.extend(["", "## Output Files"])
    outputs = summary.get("outputs", []) or []
    if outputs:
        for rec in outputs:
            lines.append(f"- `{rec.get('path')}` ({rec.get('bytes')} bytes)")
    else:
        lines.append("- No standard output files were found yet.")

    lines.extend(["", "## Caveats"])
    for note in summary.get("notes", []):
        lines.append(f"- {note}")
    lines.append("")
    return "\n".join(lines)


def write_run_summary(project_dir: str | Path) -> dict[str, Any]:
    root = Path(project_dir).expanduser().resolve()
    reports = root / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    summary = summarize_project(root)
    md = render_summary_markdown(summary)
    md_path = reports / "RUN_SUMMARY.md"
    json_path = reports / "run_summary.json"
    md_path.write_text(md)
    json_path.write_text(json.dumps(summary, indent=2, default=str))
    return {"summary": summary, "markdown": md, "markdown_path": md_path, "json_path": json_path}
