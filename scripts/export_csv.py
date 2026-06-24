#!/usr/bin/env python3
"""
export_csv.py — Excel-friendly CSV copies + merged tables for a discovery run.

For a <name>_discovery directory this writes:
  - run{i}/benchmark/validated/hits.csv   CSV copy of each per-run hit table
  - <root>/all_runs_hits.csv              every hit from every run (run_label column)
  - <root>/hit_summary.csv                per-run counts (hits, passed, six-frame vs
                                          protein-DB, unique sequences/organisms, DBs)
  - <root>/database_summary.csv           per-run database provenance, merged
Mirrors the merged CSVs into PACKAGE/01_summary_tables/ and per-run hits.csv into PACKAGE
when a PACKAGE/ exists. Never raises on a single bad file.
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import pandas as pd

from canonical import canonical_organism as _canonical_organism  # shared source of truth
from package_layout import DIRS, PER_RUN  # single source of truth for PACKAGE/ layout


def _canon_org_set(g) -> set:
    """Set of canonical organisms in a hit group (collapses same-phage accessions;
    metagenomic genomes kept distinct). Used for every 'unique organisms' count."""
    s = {_canonical_organism(r.get("organism", ""), r.get("genome_id", "")) for _, r in g.iterrows()}
    s.discard("")
    return s


def _host_from_organism(org: str) -> str:
    """'Escherichia phage X' / 'Klebsiella virus Y' -> host genus 'Escherichia'."""
    m = re.match(r"^([A-Z][a-z]+)\s+(phage|virus)\b", str(org or ""))
    return m.group(1) if m else ""


def _base_acc(gid: str) -> str:
    """Strip the version suffix from an accession so the SAME physical genome catalogued
    under versioned + unversioned ids (e.g. NC_023589.1 from RefSeq AND NC_023589 from
    INPHARED) collapses to one genome. Used so genome counts aren't inflated by cross-database
    accession aliases."""
    return re.sub(r"\.\d+$", "", str(gid or ""))


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
            # number of DATABASE RECORDS carrying this exact sequence (the same gene catalogued
            # under several accessions/DBs), NOT biological gene copies/paralogs; n_genomes and
            # n_organisms are the physical counts.
            "database_records": len(g),
            "n_genomes": g["genome_id"].map(_base_acc).nunique(),
            "n_organisms": len(_canon_org_set(g)),
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
        # count PHYSICAL genomes (base accession) so versioned+unversioned aliases of the same
        # genome across databases don't inflate n_genomes
        genomes = _uniq(_base_acc(x) for x in g.get("genome_id", pd.Series(dtype=str)))
        orgs = _uniq(g.get("organism", pd.Series(dtype=str)))   # raw names (display)
        # unique organisms by canonical identity (host-genus aliases collapsed;
        # metagenomic/unnamed fall back to genome accession)
        canon = {_canonical_organism(r.get("organism", ""), r.get("genome_id", ""))
                 for _, r in g.iterrows()}
        canon.discard("")
        rows.append({
            "representative_organism": rep.get("organism", ""),
            "representative_genome": rep.get("genome_id", ""),
            "representative_db": rep.get("db_name", ""),
            "source_type": rep.get("source_type", ""),
            # breadth = how many UNIQUE ORGANISMS carry this exact sequence (the
            # headline discovery metric; immune to the same phage appearing in
            # several databases under different accessions)
            "n_organisms": len(canon), "organisms": ";".join(orgs),
            "n_databases": len(dbs), "databases": ";".join(dbs),
            "n_genomes": len(genomes), "n_runs": len(runs), "runs": ";".join(runs),
            "n_copies": len(g),
            "domain_aa_len": rep.get("domain_aa_len", ""),
            "best_evalue": "" if pd.isna(g["_ev"].min()) else f"{g['_ev'].min():.2g}",
            "best_bit_score": "" if pd.isna(g["_bs"].max()) else round(float(g["_bs"].max()), 1),
            "confidence_tier": rep.get("confidence_tier", ""),
            "aa_sequence": seq,
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        # Order by breadth (most widespread variant first) — the discovery story.
        out = out.sort_values("n_organisms", ascending=False, kind="stable").reset_index(drop=True)
        out.insert(0, "homolog_id", [f"H{i:04d}" for i in range(1, len(out) + 1)])
    return out


def _db_hit_summary(hits: pd.DataFrame, dbsum: pd.DataFrame) -> pd.DataFrame:
    """Complete per-database summary over EVERY database searched (including those
    with zero hits — their absence is informative, e.g. a gene missing from Pfam /
    Swiss-Prot). `dbsum` is the engine's per-database record for the representative
    run (database, status, hit_count, strict_count, nt_orf_mode, runtime_seconds);
    unique sequence/organism counts come from the validated `hits`. A trailing ALL
    row is deduplicated across databases."""
    if dbsum is None or dbsum.empty:
        return pd.DataFrame()
    hit_stats = {}
    if hits is not None and not hits.empty and "db_name" in hits.columns:
        for db, g in hits.groupby("db_name", sort=False):
            useq = int(g["aa_sequence"].nunique()) if "aa_sequence" in g.columns else 0
            hit_stats[db] = (useq, len(_canon_org_set(g)))
    rows = []
    for _, d in dbsum.iterrows():
        db = str(d.get("database", ""))
        useq, uorg = hit_stats.get(db, (0, 0))
        try:
            hc = int(float(d.get("hit_count", 0) or 0))
            sc = int(float(d.get("strict_count", 0) or 0))
        except (TypeError, ValueError):
            hc = sc = 0
        rows.append({
            "database": db,
            "type": "nucleotide (six-frame)" if str(d.get("nt_orf_mode", "")) == "sixframe" else "protein",
            "status": d.get("status", ""),
            "hits": hc,
            "strict_hits": sc,
            "unique_sequences": useq,
            "unique_organisms": uorg,
            "runtime_seconds": d.get("runtime_seconds", ""),
        })
    out = pd.DataFrame(rows).sort_values("hits", ascending=False, kind="stable")
    if hits is not None and not hits.empty:
        out = pd.concat([out, pd.DataFrame([{
            "database": "ALL (deduplicated across databases)",
            "type": "", "status": "", "hits": len(hits), "strict_hits": "",
            "unique_sequences": int(hits["aa_sequence"].nunique()) if "aa_sequence" in hits.columns else 0,
            "unique_organisms": len(_canon_org_set(hits)), "runtime_seconds": "",
        }])], ignore_index=True)
    return out


def _db_barplot(dbsum: pd.DataFrame, out_dir: Path) -> list:
    """Horizontal bar chart of hits per database (every database searched; 0-hit
    DBs included). Editable SVG + PDF + 300-dpi PNG. Empty list if matplotlib absent."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        matplotlib.rcParams["svg.fonttype"] = "none"
        matplotlib.rcParams["pdf.fonttype"] = 42
        import matplotlib.pyplot as plt
    except Exception:
        return []
    rows = dbsum[dbsum["database"] != "ALL (deduplicated across databases)"] if not dbsum.empty else dbsum
    if rows is None or rows.empty:
        return []
    labels = list(rows["database"])
    vals = [int(v) for v in rows["hits"]]
    colors = ["#4C72B0" if v > 0 else "#cccccc" for v in vals]
    fig, ax = plt.subplots(figsize=(8, max(2.2, 0.42 * len(labels) + 1)))
    y = range(len(labels))
    ax.barh(list(y), vals, color=colors, edgecolor="#33373d", linewidth=0.4)
    ax.set_yticks(list(y)); ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("hits", fontsize=9)
    ax.set_title("Hits per database searched (grey = searched, 0 hits)", fontsize=10, fontweight="bold")
    for i, v in zip(y, vals):
        ax.text(v + max(vals + [1]) * 0.01, i, str(v), va="center", fontsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    made = []
    for ext in ("png", "svg", "pdf"):
        p = out_dir / f"database_hits.{ext}"
        fig.savefig(p, dpi=300 if ext == "png" else None, bbox_inches="tight")
        made.append(str(p))
    plt.close(fig)
    return made


def _write_multifastas(allh: pd.DataFrame, dedup: pd.DataFrame, discovery: Path, pkg: Path) -> list:
    """Write combined multi-FASTAs so the whole hit set is usable in one file by
    other tools (Jalview, MEGA, CD-HIT, BLAST, etc.):
      - all_hits_aa.faa / all_hits_nt.fna : every validated hit (all runs)
      - unique_homologs_aa.faa            : one per unique sequence, rich header
    Mirrored into PACKAGE/02_sequences/per_run/ when a package exists."""
    written: list = []

    def _w(path: Path, frame: pd.DataFrame, seqcol: str, header) -> int:
        if seqcol not in frame.columns:
            return -1
        n = 0
        with open(path, "w") as fh:
            for _, r in frame.iterrows():
                seq = str(r.get(seqcol, "") or "").strip().replace("*", "")
                if not seq:
                    continue
                fh.write(f">{header(r)}\n{seq}\n")
                n += 1
        return n

    h_all = lambda r: f"{r.get('hit_id','')} {r.get('organism','')} [{r.get('db_name','')}]".strip()
    if _w(discovery / "all_hits_aa.faa", allh, "aa_sequence", h_all) >= 0:
        written.append(str(discovery / "all_hits_aa.faa"))
    if _w(discovery / "all_hits_nt.fna", allh, "nt_sequence", h_all) >= 0:
        written.append(str(discovery / "all_hits_nt.fna"))
    if dedup is not None and not dedup.empty:
        h_uniq = lambda r: (f"{r.get('homolog_id','')} {r.get('representative_organism','')} "
                            f"n_organisms={r.get('n_organisms','')} n_databases={r.get('n_databases','')}").strip()
        if _w(discovery / "unique_homologs_aa.faa", dedup, "aa_sequence", h_uniq) >= 0:
            written.append(str(discovery / "unique_homologs_aa.faa"))

    if pkg and pkg.exists():
        dst = pkg / DIRS["sequences"]
        dst.mkdir(parents=True, exist_ok=True)
        for name in ("all_hits_aa.faa", "all_hits_nt.fna", "unique_homologs_aa.faa"):
            src = discovery / name
            if src.exists():
                shutil.copy2(src, dst / name)
    return written


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
            dst = pkg / DIRS["sequences"] / PER_RUN / run_id
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

        # combined multi-FASTAs (all hits + unique homologs) for downstream tools
        written += _write_multifastas(allh, dedup, discovery, pkg)

        rows = []
        for rl, g in allh.groupby("run_label"):
            rows.append({
                "run": rl,
                "total_hits": len(g),
                "passed_filter": int((g["passes_orf_filter"] == "True").sum()),
                "six_frame_hits": int((g["source_type"] == "six_frame_orf").sum()),
                "protein_db_hits": int((g["source_type"] == "annotated_protein").sum()),
                "unique_sequences": int(g["aa_sequence"].nunique()),
                "unique_organisms": len(_canon_org_set(g)),
                "databases": ";".join(sorted(x for x in g["db_name"].unique() if x)),
            })
        pd.DataFrame(rows).to_csv(discovery / "hit_summary.csv", index=False)
        written.append(str(discovery / "hit_summary.csv"))

        # Compact main-paper table from the canonical run = the one recovering the most UNIQUE
        # homologs; ties break toward the LATER (converged) round so this matches the manifest's
        # converged headline (run_frames are in run1..runN order). Selecting by raw row count with
        # a first-tie-wins max() previously picked an earlier round with fewer unique homologs.
        best_run = max(enumerate(run_frames),
                       key=lambda t: (int(t[1]["aa_sequence"].nunique()), t[0]))[1]
        paper = _paper_table(best_run)
        if not paper.empty:
            paper.to_csv(discovery / "paper_main_table.csv", index=False)
            written.append(str(discovery / "paper_main_table.csv"))

        # Complete per-database summary over EVERY database searched (incl. 0-hit
        # ones), from the engine record for the most-complete run, joined with the
        # validated-hit unique counts. Plus a bar-chart graph.
        best_label = str(best_run["run_label"].iloc[0]) if "run_label" in best_run.columns else ""
        eng = _read_tsv(discovery / f"run{best_label}" / "benchmark" / "results" / "all_database_summary.tsv")
        dbsum = _db_hit_summary(best_run, eng)
        if not dbsum.empty:
            dbsum.to_csv(discovery / "database_hit_summary.csv", index=False)
            written.append(str(discovery / "database_hit_summary.csv"))
            written += _db_barplot(dbsum, discovery)

        # Supplementary Table S1 — genome metadata (one row per PHYSICAL genome). Collapse by
        # base accession so the same genome under versioned + unversioned ids (NC_023589.1 +
        # NC_023589, RefSeq + INPHARED) is ONE row, not two — otherwise the genome count is
        # inflated (~1.3x on gp75) by cross-database accession aliases.
        meta = []
        _allh = allh.copy()
        _allh["_base"] = _allh["genome_id"].map(_base_acc)
        for base, g in _allh.groupby("_base"):
            # `.get` so offline runs (no NCBI annotation -> no 'organism' column)
            # still export this table instead of aborting the whole CSV export.
            org = next((o for o in g.get("organism", pd.Series(dtype=str)) if o), "")
            accs = ";".join(sorted(x for x in g["genome_id"].unique() if x))
            meta.append({
                "genome_id": base, "accessions": accs, "n_accessions": g["genome_id"].nunique(),
                "organism": org, "host": _host_from_organism(org),
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
        tables = pkg / DIRS["tables"]
        tables.mkdir(parents=True, exist_ok=True)
        for name in ("paper_main_table.csv", "hits_deduplicated.csv", "database_hit_summary.csv",
                     "database_hits.png", "database_hits.svg", "database_hits.pdf",
                     "hit_summary.csv",
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
