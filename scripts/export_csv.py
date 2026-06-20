#!/usr/bin/env python3
"""
export_csv.py — Excel-friendly CSV copies + merged tables for a discovery run.

For a <name>_discovery directory this writes:
  - run{i}/benchmark/validated/hits.csv   CSV copy of each per-run hit table
  - <root>/all_runs_hits.csv              every hit from every run (run_label column)
  - <root>/hit_summary.csv                per-run counts (hits, passed, six-frame vs
                                          protein-DB, unique sequences/organisms, DBs)
  - <root>/database_summary.csv           per-run database provenance, merged
Mirrors the merged CSVs into PACKAGE/00_tables/ and per-run hits.csv into PACKAGE
when a PACKAGE/ exists. Never raises on a single bad file.
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import pandas as pd


def _host_from_organism(org: str) -> str:
    """'Escherichia phage X' / 'Klebsiella virus Y' -> host genus 'Escherichia'."""
    m = re.match(r"^([A-Z][a-z]+)\s+(phage|virus)\b", str(org or ""))
    return m.group(1) if m else ""


def _read_tsv(p: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(p, sep="\t", dtype=str).fillna("")
    except Exception:
        return pd.DataFrame()


def _paper_table(df: pd.DataFrame) -> pd.DataFrame:
    """Compact main-paper table: one row per unique homolog (collapsed by exact
    sequence), with the key columns a reader needs. The full 37-column table is
    the supplementary file."""
    d = df.copy()
    d["_ev"] = pd.to_numeric(d["evalue"], errors="coerce")
    d["_bs"] = pd.to_numeric(d["bit_score"], errors="coerce")
    d["_dl"] = pd.to_numeric(d["domain_aa_len"], errors="coerce")
    rows = []
    for _, g in d.groupby("aa_sequence", sort=False):
        rep = g.sort_values("_bs", ascending=False).iloc[0]
        rows.append({
            "representative_organism": rep.get("organism", ""),
            "accession": rep.get("genome_id", ""),
            "database": rep.get("db_name", ""),
            "copies": len(g),
            "n_organisms": int(g[g["organism"] != ""]["organism"].nunique()),
            "domain_aa_len": "" if pd.isna(rep["_dl"]) else int(rep["_dl"]),
            "domain_coverage": rep.get("domain_coverage", ""),
            "best_evalue": "" if pd.isna(g["_ev"].min()) else f"{g['_ev'].min():.2g}",
            "best_bit_score": "" if pd.isna(g["_bs"].max()) else round(float(g["_bs"].max()), 1),
            "confidence_tier": rep.get("confidence_tier", ""),
            "example_hit_id": rep.get("hit_id", ""),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values("best_bit_score", ascending=False, na_position="last").reset_index(drop=True)
        out.insert(0, "rank", range(1, len(out) + 1))
    return out


def _dedup_hits(allh: pd.DataFrame) -> pd.DataFrame:
    """Collapse hits that are the SAME protein found across multiple databases or
    iterations into one row per unique sequence, recording how many — and which —
    databases/runs/genomes recovered it. A hit found in several databases is
    stronger evidence, and a single 'unique homologs' table is what most analyses
    want. Identity = exact amino-acid sequence (unambiguous for both six-frame and
    protein-database hits)."""
    if allh.empty or "aa_sequence" not in allh.columns:
        return pd.DataFrame()
    d = allh.copy()
    d["_ev"] = pd.to_numeric(d.get("evalue"), errors="coerce")
    d["_bs"] = pd.to_numeric(d.get("bit_score"), errors="coerce")

    def _uniq(series) -> list:
        return sorted({str(x) for x in series if str(x).strip() and str(x) != "nan"})

    rows = []
    for seq, g in d.groupby("aa_sequence", sort=False):
        rep = g.sort_values("_bs", ascending=False, na_position="last").iloc[0]
        dbs = _uniq(g.get("db_name", pd.Series(dtype=str)))
        runs = _uniq(g.get("run_label", pd.Series(dtype=str)))
        genomes = _uniq(g.get("genome_id", pd.Series(dtype=str)))
        orgs = _uniq(g.get("organism", pd.Series(dtype=str)))
        rows.append({
            "representative_organism": rep.get("organism", ""),
            "representative_genome": rep.get("genome_id", ""),
            "representative_db": rep.get("db_name", ""),
            "source_type": rep.get("source_type", ""),
            "n_databases": len(dbs), "databases": ";".join(dbs),
            "n_runs": len(runs), "runs": ";".join(runs),
            "n_copies": len(g), "n_genomes": len(genomes), "n_organisms": len(orgs),
            "domain_aa_len": rep.get("domain_aa_len", ""),
            "best_evalue": "" if pd.isna(g["_ev"].min()) else f"{g['_ev'].min():.2g}",
            "best_bit_score": "" if pd.isna(g["_bs"].max()) else round(float(g["_bs"].max()), 1),
            "confidence_tier": rep.get("confidence_tier", ""),
            "aa_sequence": seq,
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values("best_bit_score", ascending=False, na_position="last").reset_index(drop=True)
        out.insert(0, "homolog_id", [f"H{i:04d}" for i in range(1, len(out) + 1)])
    return out


def export(discovery: Path) -> list[str]:
    discovery = Path(discovery)
    written: list[str] = []
    pkg = discovery / "PACKAGE"

    run_frames = []
    for tsv in sorted(discovery.glob("run*/benchmark/validated/hits.tsv")):
        df = _read_tsv(tsv)
        if df.empty:
            continue
        csv = tsv.with_suffix(".csv")
        df.to_csv(csv, index=False)
        written.append(str(csv))
        run_frames.append(df)
        # mirror per-run hits.csv into the package if present
        run_id = tsv.parts[-4]  # run1 / run2 / run3
        if pkg.exists():
            dst = pkg / "02_sequences_per_run" / run_id
            dst.mkdir(parents=True, exist_ok=True)
            shutil.copy2(csv, dst / "hits.csv")

    if run_frames:
        allh = pd.concat(run_frames, ignore_index=True)
        allh.to_csv(discovery / "all_runs_hits.csv", index=False)
        written.append(str(discovery / "all_runs_hits.csv"))

        # Deduplicated view: one row per unique homolog with cross-database/run
        # provenance ("found in N databases"). The full all_runs_hits.csv is kept
        # as the complete, un-collapsed record.
        dedup = _dedup_hits(allh)
        if not dedup.empty:
            dedup.to_csv(discovery / "hits_deduplicated.csv", index=False)
            written.append(str(discovery / "hits_deduplicated.csv"))

        rows = []
        for rl, g in allh.groupby("run_label"):
            rows.append({
                "run": rl,
                "total_hits": len(g),
                "passed_filter": int((g["passes_orf_filter"] == "True").sum()),
                "six_frame_hits": int((g["source_type"] == "six_frame_orf").sum()),
                "protein_db_hits": int((g["source_type"] == "annotated_protein").sum()),
                "unique_sequences": int(g["aa_sequence"].nunique()),
                "unique_organisms": int(g[g["organism"] != ""]["organism"].nunique()),
                "databases": ";".join(sorted(x for x in g["db_name"].unique() if x)),
            })
        pd.DataFrame(rows).to_csv(discovery / "hit_summary.csv", index=False)
        written.append(str(discovery / "hit_summary.csv"))

        # Compact main-paper table from the most complete single run (collapsed
        # to unique homologs). The full all_runs_hits.csv stays as supplementary.
        best_run = max(run_frames, key=len)
        paper = _paper_table(best_run)
        if not paper.empty:
            paper.to_csv(discovery / "paper_main_table.csv", index=False)
            written.append(str(discovery / "paper_main_table.csv"))

        # Supplementary Table S1 — genome metadata (one row per genome/source)
        meta = []
        for gid, g in allh.groupby("genome_id"):
            org = next((o for o in g["organism"] if o), "")
            meta.append({
                "genome_id": gid, "organism": org, "host": _host_from_organism(org),
                "databases": ";".join(sorted(x for x in g["db_name"].unique() if x)),
                "source_type": ";".join(sorted(g["source_type"].unique())),
                "n_hits": len(g),
            })
        pd.DataFrame(meta).to_csv(discovery / "genome_metadata.csv", index=False)
        written.append(str(discovery / "genome_metadata.csv"))

        # Supplementary Table S3 — per-hit homology statistics
        s3_cols = ["hit_id", "organism", "genome_id", "db_name", "source_type",
                   "run_label", "evalue", "bit_score", "domain_aa_len",
                   "domain_coverage", "confidence_tier"]
        allh[[c for c in s3_cols if c in allh.columns]].to_csv(
            discovery / "homolog_stats.csv", index=False)
        written.append(str(discovery / "homolog_stats.csv"))

    db_frames = []
    for s in sorted(discovery.glob("run*/benchmark/results/all_database_summary.tsv")):
        df = _read_tsv(s)
        if df.empty:
            continue
        df.insert(0, "run", s.parts[-4])
        db_frames.append(df)
    if db_frames:
        pd.concat(db_frames, ignore_index=True).to_csv(discovery / "database_summary.csv", index=False)
        written.append(str(discovery / "database_summary.csv"))

    if pkg.exists():
        tables = pkg / "00_tables"
        tables.mkdir(exist_ok=True)
        for name in ("paper_main_table.csv", "hits_deduplicated.csv", "hit_summary.csv",
                     "database_summary.csv", "genome_metadata.csv", "homolog_stats.csv",
                     "all_runs_hits.csv"):
            src = discovery / name
            if src.exists():
                shutil.copy2(src, tables / name)
    return written


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discovery-dir", type=Path, required=True)
    args = ap.parse_args()
    for f in export(args.discovery_dir):
        print(f"  wrote {f}")


if __name__ == "__main__":
    main()
