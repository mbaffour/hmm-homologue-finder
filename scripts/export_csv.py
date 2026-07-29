#!/usr/bin/env python3
"""
export_csv.py — Excel-friendly CSV copies + merged tables for a discovery run.

For a <name>_discovery directory this writes:
  - run{i}/benchmark/validated/hits.csv   CSV copy of each per-run hit table
  - <root>/all_runs_hits.csv              every hit from every run (run_label column)
  - <root>/hit_summary.csv                per-run counts (hits, passed, six-frame vs
                                          protein-DB, unique sequences/organisms, DBs)
  - <root>/database_summary.csv           per-run database provenance, merged
  - <root>/stage*_summary.csv             one table per pipeline stage (stage_summary.py)
  - <root>/pipeline_stage_summary.csv     all stage tables concatenated
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
from run_selection import best_run_index, locus_ids  # one rule for locus id + canonical run
from package_layout import DIRS, PER_RUN  # single source of truth for PACKAGE/ layout
import stage_summary  # per-stage summary tables (same schema, concatenable)

# Files mirrored from the run root into PACKAGE/01_summary_tables/. Some are written by
# OTHER modules (family census, overprinting), so anything absent is simply skipped —
# listing a file here is a request, not a promise that it exists.
TABLE_EXPORTS = (
    "paper_main_table.csv", "hits_deduplicated.csv", "database_hit_summary.csv",
    "database_hits.png", "database_hits.svg", "database_hits.pdf",
    "hit_summary.csv", "database_summary.csv", "genome_metadata.csv",
    "homolog_stats.csv", "all_runs_hits.csv",
    "family_census.csv", "family_census_members.csv",
    "overprinted_loci.csv", "overprinting_summary.csv",
)


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


def _locus_ids(d: pd.DataFrame) -> list:
    """Physical-locus label per hit row — see `run_selection.locus_ids` for the rule.
    Kept as a thin adapter so the DataFrame and the pipeline paths share ONE
    implementation of what makes two hits the same homolog."""
    return locus_ids(d.to_dict("records"))


def _read_tsv(p: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(p, sep="\t", dtype=str).fillna("")
    except Exception:
        return pd.DataFrame()


def _paper_table(dedup: pd.DataFrame) -> pd.DataFrame:
    """Compact main-paper table: a projection of the deduplicated homolog table.

    It is derived from `_dedup_hits` output rather than recomputed from the raw hits, so
    the two exports cannot disagree about how many homologs were found — they previously
    used different groupings and different rounds, and shipped 55 next to 71."""
    if dedup is None or dedup.empty:
        return pd.DataFrame()
    keep = [
        ("representative_organism", "representative_organism"),
        ("representative_genome", "accession"),
        ("representative_db", "database"),
        ("n_loci", "n_loci"),                    # genomic copies of the gene
        ("n_copies", "database_records"),        # raw records, NOT independent support
        ("n_genomes", "n_genomes"),
        ("n_organisms", "n_organisms"),
        ("domain_aa_len", "domain_aa_len"),
        ("max_domain_aa_len_any_round", "max_domain_aa_len_any_round"),
        ("full_length_aa", "full_length_aa"),
        ("domain_coverage", "domain_coverage"),
        ("best_evalue", "best_evalue"),
        ("best_bit_score", "best_bit_score"),
        ("confidence_tier", "confidence_tier"),
        ("example_hit_id", "example_hit_id"),
    ]
    out = dedup[[src for src, _ in keep if src in dedup.columns]].copy()
    out.columns = [dst for src, dst in keep if src in dedup.columns]
    out = out.sort_values("best_bit_score", ascending=False, na_position="last").reset_index(drop=True)
    out.insert(0, "rank", range(1, len(out) + 1))
    return out


def _dedup_hits(allh: pd.DataFrame, all_rounds: pd.DataFrame = None) -> pd.DataFrame:
    """Collapse hits that are the SAME gene — found in several databases, or re-called
    in a later iteration — into one row per homolog LOCUS, recording how many, and
    which, databases/runs/genomes recovered it.

    Identity is the genomic locus (see `_locus_ids`), NOT the amino-acid string: the
    string is the HMM envelope slice and is re-trimmed whenever the model is refined,
    so grouping on it counted one gene up to three times.

    On independence: `n_databases` counts database *records*, and it is NOT a measure of
    corroboration — INPHARED redistributes RefSeq `NC_` records, so the same physical
    genome routinely appears in both. `n_organisms` (distinct phages) is the honest
    breadth metric and is what the table is sorted by.

    `allh` should be the CANONICAL run (see `run_selection.best_run_index`) so this table
    describes the same iteration as the deposited HMM, the controls and the figures.
    Pass every round as `all_rounds` to additionally report, per gene, the longest
    envelope any round called (`max_domain_aa_len_any_round`) — membership and the
    reported sequences still come from the canonical run, so nothing drifts, but a gene
    that an earlier round called longer is no longer silently published short."""
    if allh.empty or "aa_sequence" not in allh.columns:
        return pd.DataFrame()
    d = allh.copy()
    d["_ev"] = pd.to_numeric(d.get("evalue"), errors="coerce")
    d["_bs"] = pd.to_numeric(d.get("bit_score"), errors="coerce")
    d["_dl"] = pd.to_numeric(d.get("domain_aa_len"), errors="coerce")

    # Assign loci over EVERY round so a gene keeps one identity across iterations, then
    # keep only the loci the canonical run recovered.
    if all_rounds is not None and not all_rounds.empty:
        u = all_rounds.copy()
        u["_locus"] = _locus_ids(u)
        u["_dl_any"] = pd.to_numeric(u.get("domain_aa_len"), errors="coerce")
        longest = u.groupby("_locus")["_dl_any"].max()
        key = ["genome_id", "contig", "nt_start", "nt_end", "strand", "run_label"]
        key = [c for c in key if c in u.columns and c in d.columns]
        d = d.merge(u[key + ["_locus"]].drop_duplicates(subset=key), on=key, how="left")
        d["_locus"] = d["_locus"].fillna(pd.Series(_locus_ids(d), index=d.index))
        d["_longest_any"] = d["_locus"].map(longest)
    else:
        d["_locus"] = _locus_ids(d)
        d["_longest_any"] = d["_dl"]

    def _uniq(series) -> list:
        return sorted({str(x) for x in series if str(x).strip() and str(x) != "nan"})

    # ---- stage 1: one representative per physical gene copy (locus) ----------------
    # Collapses the two things that used to inflate the count: the same genome catalogued
    # in more than one database, and the same gene re-trimmed by a refined model in a
    # later round. Each locus keeps its BEST call (highest bit score, longest domain on
    # a tie) so a later, shorter envelope cannot shrink a reported protein.
    reps = []
    for _locus, g in d.groupby("_locus", sort=False):
        rep = g.sort_values(["_bs", "_dl"], ascending=False, na_position="last").iloc[0].copy()
        rep["_records"] = len(g)
        rep["_max_dl"] = g["_longest_any"].max()   # longest envelope in ANY round
        rep["_min_ev"] = g["_ev"].min()
        rep["_max_bs"] = g["_bs"].max()
        rep["_dbs"] = _uniq(g.get("db_name", pd.Series(dtype=str)))
        rep["_runs"] = _uniq(g.get("run_label", pd.Series(dtype=str)))
        rep["_genomes"] = _uniq(_base_acc(x) for x in g.get("genome_id", pd.Series(dtype=str)))
        reps.append(rep)
    R = pd.DataFrame(reps)
    if R.empty:
        return pd.DataFrame()

    # ---- stage 2: group gene copies that are the SAME PROTEIN ----------------------
    # One row per distinct homolog protein; `n_loci` says how many genomic copies carry
    # it. Because stage 1 already fixed the envelope drift, this count is stable across
    # iterations — grouping the raw rows on the string was what produced the inflated
    # "unique homolog" figure.
    rows = []
    for seq, g in R.groupby("aa_sequence", sort=False):
        rep = g.sort_values(["_max_bs", "_max_dl"], ascending=False, na_position="last").iloc[0]
        dbs = _uniq(x for lst in g["_dbs"] for x in lst)
        runs = _uniq(x for lst in g["_runs"] for x in lst)
        genomes = _uniq(x for lst in g["_genomes"] for x in lst)
        orgs = _uniq(g.get("organism", pd.Series(dtype=str)))    # raw names (display)
        canon = {_canonical_organism(r.get("organism", ""), r.get("genome_id", ""))
                 for _, r in g.iterrows()}
        canon.discard("")
        rows.append({
            "representative_organism": rep.get("organism", ""),
            "representative_genome": rep.get("genome_id", ""),
            "representative_db": rep.get("db_name", ""),
            "source_type": rep.get("source_type", ""),
            # breadth = how many distinct PHAGES carry this protein. This is the honest
            # corroboration metric: it is immune to one phage appearing in several
            # databases under different accessions.
            "n_organisms": len(canon), "organisms": ";".join(orgs),
            # genomic copies of the gene (one per locus)
            "n_loci": len(g),
            "n_genomes": len(genomes),
            # NOTE: database *records*, NOT independent corroboration — INPHARED
            # redistributes RefSeq NC_ records, so one physical genome routinely appears
            # in both. Use n_organisms for breadth claims.
            "n_databases": len(dbs), "databases": ";".join(dbs),
            "n_runs": len(runs), "runs": ";".join(runs),
            "n_copies": int(g["_records"].sum()),
            "domain_aa_len": rep.get("domain_aa_len", ""),
            # longest envelope any round called for this gene. Publishing only the final
            # round's call reported two Erwinia proteins ~43 % short of this.
            "max_domain_aa_len_any_round": "" if pd.isna(g["_max_dl"].max()) else int(g["_max_dl"].max()),
            "full_length_aa": rep.get("orf_aa_len", ""),
            "domain_coverage": rep.get("domain_coverage", ""),
            "best_evalue": "" if pd.isna(g["_min_ev"].min()) else f"{g['_min_ev'].min():.2g}",
            "best_bit_score": "" if pd.isna(g["_max_bs"].max()) else round(float(g["_max_bs"].max()), 1),
            "confidence_tier": rep.get("confidence_tier", ""),
            "example_hit_id": rep.get("hit_id", ""),
            "aa_sequence": seq,
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        # Order by breadth (most widespread protein first) — the discovery story.
        out = out.sort_values(["n_organisms", "best_bit_score"], ascending=False,
                              kind="stable").reset_index(drop=True)
        out.insert(0, "homolog_id", [f"H{i:04d}" for i in range(1, len(out) + 1)])
    return out


def _db_hit_summary(hits: pd.DataFrame, dbsum: pd.DataFrame) -> pd.DataFrame:
    """Complete per-database summary over EVERY database searched (including those
    with zero hits — their absence is informative, e.g. a gene missing from Pfam /
    Swiss-Prot). `dbsum` is the engine's per-database record for the representative
    run (database, status, hit_count, strict_count, nt_orf_mode, runtime_seconds);
    unique sequence/organism counts come from the validated `hits`. A trailing ALL
    row is deduplicated across databases.

    Every "unique_sequences" figure here goes through `_dedup_hits`, so this file, the
    homolog table and the paper table all count homologs the same way. Counting distinct
    `aa_sequence` strings directly (as this used to) counts one gene once per envelope
    the model happened to cut, which is how the package ended up quoting 55 in one file
    and 71 in another."""
    if dbsum is None or dbsum.empty:
        return pd.DataFrame()
    hit_stats = {}
    if hits is not None and not hits.empty and "db_name" in hits.columns:
        for db, g in hits.groupby("db_name", sort=False):
            hit_stats[db] = (len(_dedup_hits(g)), len(_canon_org_set(g)))
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
            "unique_sequences": len(_dedup_hits(hits)),
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

        # Canonical run = the iteration the WHOLE package describes, chosen by the shared
        # rule in run_selection so the deposited HMM, the controls, the figures and the
        # tables can no longer disagree about which round they came from.
        # run_frames are in run1..runN order, so index i -> run i+1.
        _by_idx = {i + 1: f.to_dict("records") for i, f in enumerate(run_frames)}
        best_i = best_run_index(_by_idx)
        best_run = run_frames[best_i - 1]

        # Deduplicated view: one row per distinct homolog protein, with cross-database
        # provenance. Built from the CANONICAL run so its count matches the deposited HMM,
        # the tree and the controls; `all_rounds` only adds the longest-envelope audit
        # column. The full all_runs_hits.csv stays as the complete, un-collapsed record.
        dedup = _dedup_hits(best_run, all_rounds=allh)
        if not dedup.empty:
            dedup.to_csv(discovery / "hits_deduplicated.csv", index=False)
            written.append(str(discovery / "hits_deduplicated.csv"))

        # combined multi-FASTAs (all hits + unique homologs) for downstream tools
        written += _write_multifastas(allh, dedup, discovery, pkg)

        rows = []
        for rl, g in allh.groupby("run_label"):
            # Optional columns are OMITTED when absent, never defaulted. `g["passes_orf_filter"]`
            # used to be indexed directly, so a hits.tsv written without it raised KeyError and
            # took the whole export down — no stage tables, no package, from one missing column.
            # Filling 0 instead would be worse: it would report "no hit passed validation" out of
            # a check that never ran (the defect stage_summary._stage03 guards against).
            row = {"run": rl, "total_hits": len(g)}
            if "passes_orf_filter" in g.columns:
                row["passed_filter"] = int((g["passes_orf_filter"] == "True").sum())
            if "source_type" in g.columns:
                row["six_frame_hits"] = int((g["source_type"] == "six_frame_orf").sum())
                row["protein_db_hits"] = int((g["source_type"] == "annotated_protein").sum())
            row["unique_sequences"] = len(_dedup_hits(g))   # same definition as every other table
            row["unique_organisms"] = len(_canon_org_set(g))
            if "db_name" in g.columns:
                row["databases"] = ";".join(sorted(x for x in g["db_name"].unique() if x))
            rows.append(row)
        pd.DataFrame(rows).to_csv(discovery / "hit_summary.csv", index=False)
        written.append(str(discovery / "hit_summary.csv"))

        # Projection of the SAME deduplicated table, so the two exports cannot disagree.
        paper = _paper_table(dedup)
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
                   "run_label", "evalue", "bit_score", "orf_aa_len", "domain_aa_len",
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

    # Per-stage summary tables (one shared schema) + the concatenated pipeline table.
    # Written BEFORE the mirror block so the freshly-built stage*_summary.csv files are
    # picked up by the glob below on this same pass, not one export late.
    written += stage_summary.build(discovery)

    if pkg.exists():
        tables = pkg / DIRS["tables"]
        tables.mkdir(parents=True, exist_ok=True)

        # The mirror is AUTHORITATIVE for the stage tables: whatever is in the run root is
        # what PACKAGE ships, including the absence of a table.
        #
        # This deletion pass is load-bearing, not tidiness. stage_summary.build() unlinks the
        # table of a stage that produced nothing (controls removed, --find-interrupted
        # dropped, a stage that failed this time), but nothing else in the pipeline ever
        # removes a file from PACKAGE: the loop below only copies, assemble_package does
        # pkg.mkdir(exist_ok=True) rather than wiping, and exporting into an existing output
        # directory is allowed. Without this, a re-export left the PREVIOUS run's
        # stage05_controls_summary.csv sitting in PACKAGE/01_summary_tables/ — a confident
        # table of numbers from a step that did not run.
        stale = list(tables.glob("stage*_summary.csv")) + [tables / stage_summary.SUMMARY_NAME]
        for p in stale:
            if p.exists() and not (discovery / p.name).exists():
                try:
                    p.unlink()
                except OSError:
                    pass

        # Named tables + every stage table, discovered by glob so adding a stage to
        # stage_summary.py never needs an edit here. sorted() keeps stage00..stage08 order.
        names = list(TABLE_EXPORTS) + sorted(
            p.name for p in discovery.glob("stage*_summary.csv")) + [stage_summary.SUMMARY_NAME]
        for name in dict.fromkeys(names):        # de-dup, preserve order
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
