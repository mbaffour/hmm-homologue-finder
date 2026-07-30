#!/usr/bin/env python3
"""
stage_summary.py — one summary table per pipeline stage, all in ONE shared schema.

    stage, stage_name, metric, value, units, note, source_file

Every stage table has exactly those seven columns. That is the whole point: it makes the
tables *dictionary-able* (one entry in the data dictionary covers all of them) and
*concatenable*, so `pipeline_stage_summary.csv` — the single table a reviewer reads — is
just the per-stage tables stacked, with nothing reformatted on the way.

Stages (a `_stageNN(discovery) -> list[dict]` each):

    00 input       parameters, seed checksum, how many sequences went in
    01 model       per-iteration counts, HMM length, why iteration stopped, canonical run
    02 search      rollup of database_hit_summary.csv (every database, incl. 0-hit ones)
    03 validation  ORF filter + confidence tiers on the canonical run's hits.tsv
    04 homologs    hits_deduplicated.csv — proteins, loci, organisms, genomes
    05 controls    control_report.json + sixframe_decoy_control.json (incl. the decoy FDR)
    06 seeds       seed_qc/seed_recovery.csv (+ family_census.csv when present)
    07 overprint   overprinted_loci.csv / overprinting_summary.csv (or interrupted_homologs.tsv)
    08 downstream  alignment stats, tree tips, synteny clusters

Each stage is wrapped in its own try/except returning `[]`, so a run without controls, or
without `--find-interrupted`, or stopped before the tree, still exports every stage it DOES
have instead of failing the whole export.

HARD CONSTRAINT — `build()` is pure filesystem: no network, no subprocesses, no re-search.
`export_csv.export()` calls it and export() runs TWICE per pipeline run (once mid-run, once
at the end), so anything expensive or side-effecting here would be paid for twice and could
change the run it is supposed to be describing.
"""
from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

import pandas as pd

from canonical import canonical_organism            # one definition of "same phage"
from run_selection import best_run_index, n_unique_loci, read_hits_rows

# The seven canonical columns, in order. Every row every stage emits has exactly these keys —
# `build()` reindexes to this tuple so a stage that forgets one still writes a valid table.
COLUMNS = ("stage", "stage_name", "metric", "value", "units", "note", "source_file")

SUMMARY_NAME = "pipeline_stage_summary.csv"   # the concatenation of every stage table

__all__ = ["COLUMNS", "SUMMARY_NAME", "build", "stage_rows"]


# --------------------------------------------------------------------------------------
# small, total helpers — none of these raise
# --------------------------------------------------------------------------------------
def _fmt(v):
    """Render a value for CSV, ALWAYS as text.

    Text on purpose: pandas types a column by its contents, so a stage whose values happen
    to be all-numeric with one float in the mix would re-render every integer count as
    '3.0'. The same metric must not change appearance depending on what else is in the
    table. Floats are rounded first (a bit score is not meaningful to 12 decimals) and
    NaN/inf become blank rather than the literal 'nan', which would read as a measurement."""
    if v is None:
        return ""
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if v != v or v in (float("inf"), float("-inf")):
            return ""
        return str(round(v, 6))
    return str(v)


def _row(stage, stage_name, metric, value, units="", note="", source_file=""):
    return {"stage": stage, "stage_name": stage_name, "metric": metric,
            "value": _fmt(value), "units": units, "note": note, "source_file": source_file}


def _rel(discovery: Path, p: Path) -> str:
    """Path relative to the run directory (POSIX), so the tables read the same on every OS
    and never leak an absolute user path into a shipped file."""
    try:
        return Path(p).resolve().relative_to(Path(discovery).resolve()).as_posix()
    except Exception:
        return Path(p).name


def _read_json(p: Path) -> dict:
    try:
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            d = json.load(fh)
        return d if isinstance(d, dict) else {}
    except Exception:
        return {}


def _read_table(p: Path, sep: str = ",") -> pd.DataFrame:
    """Read a CSV/TSV as strings. Everything here is counted or re-parsed explicitly, and
    string dtype keeps accessions like `NC_008720.10` from being mangled into floats."""
    try:
        return pd.read_csv(p, sep=sep, dtype=str).fillna("")
    except Exception:
        return pd.DataFrame()


def _num(x):
    try:
        v = float(x)
        return None if v != v else v
    except (TypeError, ValueError):
        return None


def _int(x, default=0) -> int:
    v = _num(x)
    return default if v is None else int(v)


def _slug(s: str) -> str:
    """Database/control names become metric-name-safe: 'RefSeq viral genomes' ->
    'refseq_viral_genomes'. The unmangled name is kept in the `note` column."""
    return re.sub(r"[^a-z0-9]+", "_", str(s or "").lower()).strip("_") or "unnamed"


def _counts(series) -> list:
    """(value, count) pairs sorted by count desc, blanks dropped."""
    vals = [str(x).strip() for x in series if str(x).strip()]
    out = {}
    for v in vals:
        out[v] = out.get(v, 0) + 1
    return sorted(out.items(), key=lambda kv: (-kv[1], kv[0]))


def _hmm_stats(hmm: Path) -> dict:
    """LENG (match states) and NSEQ (sequences the model was built from) out of an HMMER
    profile header. Read line-wise and bailed out of early — a profile is megabytes of
    emission table we have no reason to load."""
    out = {}
    try:
        with open(hmm, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("HMM "):
                    break                      # header is over; the matrix starts here
                m = re.match(r"^(LENG|NSEQ|NAME)\s+(\S+)", line)
                if m:
                    out[m.group(1)] = m.group(2)
    except Exception:
        return {}
    return out


def _canonical_run(discovery: Path):
    """(run_index, hits.tsv path, hits frame) for the iteration the package describes.

    Uses `run_selection.best_run_index` — the SAME rule as `export_csv` and `hmm_finder` —
    so a stage table can never quote a different round than the deposited HMM, the controls
    and the homolog tables. Returns (None, None, empty) when there are no runs yet."""
    runs = {}
    for tsv in sorted(discovery.glob("run*/benchmark/validated/hits.tsv")):
        m = re.search(r"run(\d+)", tsv.parts[-4])
        if m:
            runs[int(m.group(1))] = tsv
    if not runs:
        return (None, None, pd.DataFrame())
    rows = {i: read_hits_rows(p) for i, p in runs.items()}
    best = best_run_index(rows)
    return (best, runs.get(best), pd.DataFrame(rows.get(best) or []))


def _base_acc(gid: str) -> str:
    """Drop a version suffix so NC_023589.1 and NC_023589 are ONE genome (the same physical
    sequence is catalogued both ways by RefSeq and INPHARED)."""
    return re.sub(r"\.\d+$", "", str(gid or ""))


def _gbk_organism(p: Path) -> str:
    """ORGANISM name out of a GenBank record, '' if there is none.

    Only the header block is read (the loop bails at FEATURES, ~line 9): these records carry
    the whole neighbourhood nucleotide sequence and `build()` must stay cheap."""
    try:
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("FEATURES"):
                    break
                if line.startswith("  ORGANISM  "):
                    return line[len("  ORGANISM  "):].strip()
    except OSError:
        return ""
    return ""


# Engine status words for a database search. A database that FAILED or was SKIPPED returned
# zero hits because it was never searched — that zero is a tooling outcome, not evidence the
# gene is absent — so the two are counted separately from the databases that completed.
_DB_OK = {"complete", "completed", "ok", "success", "succeeded", "done", "finished"}
_DB_FAILED = {"failed", "fail", "error", "errored", "timeout", "timed_out", "crashed"}
_DB_SKIPPED = {"skipped", "skip", "missing", "absent", "unavailable", "not_found", "disabled"}


# --------------------------------------------------------------------------------------
# stage 00 — input
# --------------------------------------------------------------------------------------
def _stage00(discovery: Path) -> list:
    S, N = "00", "input"
    man = _read_json(discovery / "run_manifest.json")
    if not man:
        return []                       # manifest is written at the END of a run
    src = "run_manifest.json"
    p = man.get("parameters") or {}
    inp = man.get("input") or {}
    dbs = [d.strip() for d in str(p.get("databases", "")).split(",") if d.strip()]

    out = [
        _row(S, N, "run_label", p.get("label", ""), "", "user label for this run", src),
        _row(S, N, "tool", man.get("tool", ""), "", man.get("code_git_commit", "") and
             f"git commit {man.get('code_git_commit')}", src),
        _row(S, N, "n_input_seeds", _int(man.get("n_input_seeds")), "sequences",
             "FASTA records submitted — headers, NOT distinct proteins (exact duplicates "
             "are counted here and collapsed by the family census)", src),
        _row(S, N, "seed_fasta", inp.get("fasta", ""), "path",
             "user paths are redacted before the manifest ships", src),
        _row(S, N, "seed_fasta_sha256", inp.get("sha256", ""), "sha256",
             "checksum of the exact input file; a rerun that does not match is not a rerun", src),
        _row(S, N, "iterations_requested", p.get("iterations", ""), "iterations",
             "upper bound — the search may converge earlier (see stage 01)", src),
        _row(S, N, "databases_requested", len(dbs), "databases", "; ".join(dbs), src),
        _row(S, N, "translation_table", p.get("trans_table", ""), "NCBI table",
             "used for six-frame translation of nucleotide databases", src),
        _row(S, N, "prodigal_gate", p.get("prodigal_gate", ""), "flag",
             "False = hits are validated on six-frame ORF evidence, not on agreement with a "
             "gene caller (the target gene is antisense inside another CDS, so a caller "
             "does not annotate it)", src),
        _row(S, N, "min_recovery", p.get("min_recovery", ""), "fraction",
             "minimum fraction of seeds the refined model must still recover", src),
        _row(S, N, "cpu", p.get("cpu", ""), "threads", "", src),
        _row(S, N, "python", man.get("python", ""), "version", man.get("conda_env", ""), src),
        _row(S, N, "started_at", man.get("started_at", ""), "timestamp", "", src),
        _row(S, N, "finished_at", man.get("finished_at", ""), "timestamp", "", src),
    ]
    # wall-clock is derived, not stored — only report it when BOTH stamps parse
    try:
        t0 = datetime.strptime(str(man.get("started_at", "")), "%Y-%m-%d %H:%M:%S")
        t1 = datetime.strptime(str(man.get("finished_at", "")), "%Y-%m-%d %H:%M:%S")
        out.append(_row(S, N, "wall_clock", round((t1 - t0).total_seconds() / 60.0, 1),
                        "minutes", "finished_at - started_at", src))
    except (ValueError, TypeError):
        pass
    return [r for r in out if str(r["value"]) != ""]


# --------------------------------------------------------------------------------------
# stage 01 — model
# --------------------------------------------------------------------------------------
def _stage01(discovery: Path) -> list:
    S, N = "01", "model"
    man = _read_json(discovery / "run_manifest.json")
    out = []

    # Per-iteration counts are recomputed from the run directories rather than copied from
    # the manifest, so this table stays correct for a run that was interrupted before the
    # manifest was written (export() runs mid-pipeline too).
    per_run = {}
    for tsv in sorted(discovery.glob("run*/benchmark/validated/hits.tsv")):
        m = re.search(r"run(\d+)", tsv.parts[-4])
        if m:
            per_run[int(m.group(1))] = tsv
    for i in sorted(per_run):
        tsv = per_run[i]
        rows = read_hits_rows(tsv)
        src = _rel(discovery, tsv)
        out.append(_row(S, N, f"iteration{i}_hits", len(rows), "hits",
                        "validated hit ROWS (one genome copy can appear in several databases)", src))
        out.append(_row(S, N, f"iteration{i}_unique_loci", n_unique_loci(rows), "loci",
                        "distinct genomic loci — the homolog identity used everywhere "
                        "(same organism + strand + overlapping ORF interval)", src))
        hmm = discovery / f"run{i}" / "benchmark" / "hmm" / "benchmark_profile.hmm"
        st = _hmm_stats(hmm)
        if st.get("LENG"):
            out.append(_row(S, N, f"iteration{i}_hmm_length", _int(st["LENG"]), "match states",
                            "profile HMM length (LENG) for this iteration", _rel(discovery, hmm)))
        if st.get("NSEQ"):
            out.append(_row(S, N, f"iteration{i}_hmm_nseq", _int(st["NSEQ"]), "sequences",
                            "sequences the iteration's model was built from (NSEQ)",
                            _rel(discovery, hmm)))

    if per_run:
        out.append(_row(S, N, "iterations_completed", len(per_run), "iterations", "", ""))

    best, best_tsv, _ = _canonical_run(discovery)
    if best is not None:
        out.append(_row(S, N, "canonical_run", best, "iteration",
                        "the iteration the WHOLE package describes: most distinct loci, ties "
                        "to the later (converged) round — run_selection.best_run_index",
                        _rel(discovery, best_tsv) if best_tsv else ""))
        # The deposited profile is what a reader will actually download, so quote ITS length.
        dep = discovery / "PACKAGE" / "03_hmm_profile" / "profile.hmm"
        hmm = dep if dep.exists() else discovery / f"run{best}" / "benchmark" / "hmm" / "benchmark_profile.hmm"
        st = _hmm_stats(hmm)
        if st.get("LENG"):
            out.append(_row(S, N, "canonical_hmm_length", _int(st["LENG"]), "match states",
                            "length of the deposited profile.hmm", _rel(discovery, hmm)))
        if st.get("NSEQ"):
            out.append(_row(S, N, "canonical_hmm_nseq", _int(st["NSEQ"]), "sequences",
                            "sequences the deposited model was built from", _rel(discovery, hmm)))

    if man.get("iteration_stop_reason"):
        out.append(_row(S, N, "iteration_stop_reason", man["iteration_stop_reason"], "text",
                        "why iterating stopped (convergence, or the iteration cap)",
                        "run_manifest.json"))
    return out


# --------------------------------------------------------------------------------------
# stage 02 — search
# --------------------------------------------------------------------------------------
def _stage02(discovery: Path) -> list:
    S, N = "02", "search"
    f = discovery / "database_hit_summary.csv"
    db = _read_table(f)
    if db.empty or "database" not in db.columns:
        return []
    src = "database_hit_summary.csv"
    ALL = "ALL (deduplicated across databases)"
    per = db[db["database"] != ALL]
    tot = db[db["database"] == ALL]

    # Split the databases by engine status BEFORE counting anything. A failed or skipped
    # database also shows 0 hits, and rolling it in with the completed ones presented a
    # tooling failure as negative biological evidence ("searched, found nothing").
    status = {i: str(r.get("status", "") or "").strip().lower() for i, r in per.iterrows()} \
        if "status" in per.columns else {}
    done = [i for i, s in status.items() if s in _DB_OK]
    failed = [i for i, s in status.items() if s in _DB_FAILED]
    skipped = [i for i, s in status.items() if s in _DB_SKIPPED]
    unknown = [i for i, s in status.items() if s not in _DB_OK | _DB_FAILED | _DB_SKIPPED]

    out = [
        _row(S, N, "databases_attempted", len(per), "databases",
             "databases the engine was asked to search in the canonical iteration, whatever "
             "the outcome", src),
    ]
    if status:
        out.append(
            _row(S, N, "databases_searched", len(done), "databases",
                 "databases that COMPLETED, including those returning 0 hits — for these, and "
                 "only these, a 0 is a result (the gene is absent from the protein databases "
                 "because it is unannotated). Failed and skipped databases are counted "
                 "separately below and their 0 means nothing was searched", src))
        out.append(
            _row(S, N, "databases_failed", len(failed), "databases",
                 "the search errored out; a 0 hit count here is NOT evidence of absence", src))
        out.append(
            _row(S, N, "databases_skipped", len(skipped), "databases",
                 "never searched (database missing or disabled); a 0 hit count here is NOT "
                 "evidence of absence", src))
        if unknown:
            out.append(
                _row(S, N, "databases_status_unknown", len(unknown), "databases",
                     "status word the engine reported is not one this summary recognises: "
                     + "; ".join(sorted({status[i] or "(blank)" for i in unknown}))
                     + " — treat their hit counts as unverified", src))
    else:
        # No status column at all (an older run): report that instead of publishing a
        # confident databases_searched that would be indistinguishable from "all completed".
        out.append(
            _row(S, N, "databases_status_unknown", len(per), "databases",
                 "this run's database_hit_summary.csv carries no status column, so which "
                 "databases actually completed cannot be established — no 0 here may be read "
                 "as absence", src))

    # State the SEARCH SPACE THAT WAS NOT COVERED. A reviewer asking "did you look in the gut
    # phage catalogues, or in bacterial genomes for a prophage copy?" otherwise has no answer in
    # the package: the tables list only what WAS searched, so an unsearched catalogue is
    # indistinguishable from one that returned nothing. Naming them is what makes the coverage
    # claim falsifiable.
    try:
        import sys as _sys
        _eng = str(Path(__file__).resolve().parents[1] / "engine")
        if _eng not in _sys.path:
            _sys.path.insert(0, _eng)
        from databases.builtin import BUILTIN_DATABASES        # a plain dict, no I/O
        attempted = {str(x) for x in per.get("database", [])}
        missing = [str(d.get("name")) for d in BUILTIN_DATABASES
                   if str(d.get("name")) not in attempted]
        if missing:
            out.append(_row(
                S, N, "catalog_databases_not_searched", len(missing), "databases",
                "in the catalog but NOT searched in this run, so this run says nothing either "
                "way about them: " + "; ".join(missing)
                + " — a homolog present only in one of these would not appear above",
                "engine/databases/builtin.py"))
    except Exception:
        pass

    out += [
        _row(S, N, "databases_with_hits", int(sum(1 for h in per.get("hits", []) if _int(h) > 0)),
             "databases", "", src),
        _row(S, N, "nucleotide_databases",
             int(sum(1 for t in per.get("type", []) if "nucleotide" in str(t))), "databases",
             "of the databases attempted; searched by six-frame translation", src),
        _row(S, N, "protein_databases",
             int(sum(1 for t in per.get("type", []) if str(t) == "protein")), "databases",
             "of the databases attempted", src),
        _row(S, N, "total_hits", sum(_int(h) for h in per.get("hits", [])), "hits",
             "hit rows summed over databases (the same genome copy in two databases counts twice)", src),
        _row(S, N, "search_runtime", round(sum(_num(r) or 0.0 for r in per.get("runtime_seconds", [])), 1),
             "seconds", "canonical iteration only", src),
    ]
    if not tot.empty:
        t = tot.iloc[0]
        out += [
            _row(S, N, "unique_homolog_proteins", _int(t.get("unique_sequences")), "proteins",
                 "deduplicated ACROSS databases, by locus then protein — not a sum of the "
                 "per-database counts", src),
            _row(S, N, "unique_organisms", _int(t.get("unique_organisms")), "organisms",
                 "distinct phages (canonical organism), immune to one phage being catalogued "
                 "in several databases", src),
        ]
    # per-database rollup: one row per database, so the table shows WHICH database found what
    for i, r in per.iterrows():
        name = str(r.get("database", ""))
        # A 0 from a database that did not complete must not be readable as absence, so the
        # caveat rides on the row itself and not only on the counts above.
        warn = "" if (i in done or (not status)) else \
            " — NOT searched successfully; this hit count is not evidence of absence"
        out.append(_row(S, N, f"hits_in_{_slug(name)}", _int(r.get("hits")), "hits",
                        f"{name} [{r.get('type','')}; status={r.get('status','')}; "
                        f"{_int(r.get('unique_sequences'))} unique proteins]{warn}", src))
    return out


# --------------------------------------------------------------------------------------
# stage 03 — validation
# --------------------------------------------------------------------------------------
def _stage03(discovery: Path) -> list:
    S, N = "03", "validation"
    best, tsv, d = _canonical_run(discovery)
    if d.empty:
        return []
    src = _rel(discovery, tsv) if tsv else ""
    n = len(d)
    out = [_row(S, N, "hits_examined", n, "hits", f"canonical run{best}", src)]
    # The whole ORF-filter block is emitted ONLY when the column exists. `.get` would hand
    # back the default for an absent column, and a hits.tsv written without
    # passes_orf_filter would then publish "pass_rate 0.0 %" — a confident statement that no
    # hit passed validation, from a check that never ran. Same treatment as the optional
    # ORF sanity flags below.
    if "passes_orf_filter" in d.columns:
        passed = int(sum(1 for v in d["passes_orf_filter"] if str(v) == "True"))
        out += [
            _row(S, N, "passes_orf_filter", passed, "hits",
                 "hit sits in a real stop-to-stop ORF of the expected length/coverage — the "
                 "six-frame evidence test, which does NOT require a gene caller to have "
                 "annotated the locus", src),
            _row(S, N, "fails_orf_filter", n - passed, "hits", "", src),
            _row(S, N, "pass_rate", round(100.0 * passed / n, 1) if n else "", "%", "", src),
        ]
    for col, unit, note in (("source_type", "hits", "how the hit was found"),
                            ("confidence_tier", "hits", "tier assigned from E-value/bit/coverage")):
        for val, cnt in _counts(d.get(col, pd.Series(dtype=str))):
            out.append(_row(S, N, f"{col}_{_slug(val)}", cnt, unit, f"{note}: {val}", src))
    # ORF sanity flags — these are what make a six-frame call believable on its own
    for col, note in (("has_start_M", "ORF begins at a methionine"),
                      ("ends_at_stop", "ORF runs to a stop codon"),
                      ("internal_stops", "hits with >0 internal stops (read-through candidates)")):
        if col not in d.columns:
            continue
        if col == "internal_stops":
            cnt = int(sum(1 for v in d[col] if (_num(v) or 0) > 0))
        else:
            cnt = int(sum(1 for v in d[col] if str(v) == "True"))
        out.append(_row(S, N, col, cnt, "hits", note, src))
    return out


# --------------------------------------------------------------------------------------
# stage 04 — homologs
# --------------------------------------------------------------------------------------
def _stage04(discovery: Path) -> list:
    S, N = "04", "homologs"
    f = discovery / "hits_deduplicated.csv"
    d = _read_table(f)
    if d.empty:
        return []
    src = "hits_deduplicated.csv"
    out = [
        _row(S, N, "homolog_proteins", len(d), "proteins",
             "one row per distinct homolog PROTEIN in the canonical run", src),
    ]
    if "n_loci" in d.columns:
        out.append(_row(S, N, "genomic_loci", sum(_int(v) for v in d["n_loci"]), "loci",
                        "genomic gene copies — several phages can carry an identical protein, "
                        "so loci >= proteins", src))
    if "n_genomes" in d.columns:
        out.append(_row(S, N, "genome_records", sum(_int(v) for v in d["n_genomes"]), "genomes",
                        "n_genomes summed over proteins — a genome carrying two different "
                        "homolog proteins is counted twice; use distinct_organisms for breadth", src))

    # Organism/genome breadth is recomputed from the canonical run's hits rather than summed
    # out of the deduplicated table: summing per-protein counts double-counts a phage that
    # carries two different homolog proteins.
    _best, tsv, h = _canonical_run(discovery)
    if not h.empty:
        orgs = {canonical_organism(r.get("organism", ""), r.get("genome_id", ""))
                for _, r in h.iterrows()}
        orgs.discard("")
        genomes = {_base_acc(g) for g in h.get("genome_id", pd.Series(dtype=str)) if str(g).strip()}
        rsrc = _rel(discovery, tsv) if tsv else ""
        out.append(_row(S, N, "distinct_organisms", len(orgs), "organisms",
                        "distinct phages carrying a homolog (canonical organism)", rsrc))
        out.append(_row(S, N, "distinct_genomes", len(genomes), "genomes",
                        "distinct genome accessions (version suffix stripped, so NC_023589.1 "
                        "and NC_023589 are one). Higher than distinct_organisms because the "
                        "SAME phage is catalogued under a RefSeq and a GenBank accession — "
                        "breadth claims must use distinct_organisms", rsrc))

    for col, unit, note in (("domain_aa_len", "aa", "HMM-matched domain length"),
                            ("full_length_aa", "aa", "surrounding ORF length"),
                            ("best_bit_score", "bits", "best HMMER bit score")):
        vals = [v for v in (_num(x) for x in d.get(col, [])) if v is not None]
        if not vals:
            continue
        vals.sort()
        # residue counts are integers; bit scores are not, so only the lengths are collapsed
        cast = int if unit == "aa" else (lambda x: x)
        out.append(_row(S, N, f"{col}_min", cast(vals[0]), unit, note, src))
        out.append(_row(S, N, f"{col}_median", cast(vals[len(vals) // 2]), unit, note, src))
        out.append(_row(S, N, f"{col}_max", cast(vals[-1]), unit, note, src))

    for val, cnt in _counts(d.get("confidence_tier", pd.Series(dtype=str))):
        out.append(_row(S, N, f"tier_{_slug(val)}", cnt, "proteins",
                        f"homolog proteins in tier '{val}'", src))
    # A protein re-trimmed shorter by a later model would otherwise be published short; this
    # says how many rows carry a longer envelope from an earlier round.
    if {"domain_aa_len", "max_domain_aa_len_any_round"} <= set(d.columns):
        longer = int(sum(1 for _, r in d.iterrows()
                         if (_num(r["max_domain_aa_len_any_round"]) or 0) > (_num(r["domain_aa_len"]) or 0)))
        out.append(_row(S, N, "longer_envelope_in_earlier_round", longer, "proteins",
                        "canonical round trimmed these shorter than an earlier round did — "
                        "max_domain_aa_len_any_round records the longest call", src))
    return out


# --------------------------------------------------------------------------------------
# stage 05 — controls
# --------------------------------------------------------------------------------------
def _stage05(discovery: Path) -> list:
    S, N = "05", "controls"
    out = []
    rep = _read_json(discovery / "controls" / "control_report.json")
    if rep:
        src = "controls/control_report.json"
        roc = rep.get("roc") or {}
        for key, unit, note in (
            ("sensitivity", "fraction", "positive controls (the seeds) detected"),
            ("specificity", "fraction", "negative controls correctly rejected"),
            ("false_positive_rate", "fraction", ""),
            ("total_positives", "sequences", ""),
            ("true_positives", "sequences", ""),
            ("total_negatives", "sequences", "fungal + mammalian + archaeal + shuffled seeds"),
            ("false_positives", "sequences", ""),
            ("n_controls", "control sets", ""),
            ("sensitivity_strict", "fraction", "at the strict bit-score threshold"),
            ("specificity_strict", "fraction", "at the strict bit-score threshold"),
        ):
            if key in rep:
                out.append(_row(S, N, key, rep[key], unit, note, src))
        for key, unit, note in (
            ("auc", "AUC", "ROC area under the curve over the control sets"),
            ("strict_threshold", "bits", "bit score used as the strict cutoff"),
            ("optimal_threshold", "bits", str(roc.get("threshold_method", ""))),
            ("n_positive", "sequences", ""),
            ("n_negative", "sequences", ""),
            ("separable", "flag", "positives and negatives do not overlap in score"),
        ):
            if key in roc:
                out.append(_row(S, N, f"roc_{key}", roc[key], unit, note, src))
        # per-control-set detail: the numbers a referee spot-checks
        for c in (rep.get("results") or []):
            if not isinstance(c, dict):
                continue
            nm = _slug(c.get("control_name", ""))
            out.append(_row(S, N, f"control_{nm}_hit_rate", c.get("hit_rate_pct", ""), "%",
                            f"{c.get('role','')}: {c.get('desc','')} "
                            f"({_int(c.get('n_hits'))}/{_int(c.get('n_seqs'))} sequences)", src))

    dec = _read_json(discovery / "controls" / "sixframe_decoy_control.json")
    if dec and str(dec.get("status", "")) == "ok":
        src = "controls/sixframe_decoy_control.json"
        note_decoy = ("reversed six-frame ORFs from the SAME genomes: same length and "
                      "composition, no homology — the only control that measures the false "
                      "discovery rate of a six-frame search")
        for key, unit, note in (
            ("n_decoy_sequences", "sequences", note_decoy),
            ("n_sixframe_orfs_total", "ORFs", "size of the real six-frame search space"),
            ("sampled_fraction", "fraction", "decoys sampled / total six-frame ORFs"),
            ("threshold", "bits", "strict bit-score cutoff the FDR is quoted at"),
            ("best_decoy_bit_score", "bits", "the single best-scoring decoy"),
            ("weakest_true_positive_bit_score", "bits", "the weakest accepted real hit"),
            ("gap_bits", "bits", "weakest true positive - best decoy; >0 = clean separation"),
            ("clean_separation", "flag", "no decoy scores as high as any accepted hit"),
            ("decoys_at_or_above_threshold", "sequences", ""),
            ("expected_decoys_in_full_search_space", "sequences",
             "decoys at/above threshold, scaled from the sample to the whole search space"),
            ("true_positives_at_or_above_threshold", "hits", ""),
            ("empirical_fdr", "fraction", "expected decoys / (expected decoys + true positives)"),
            ("hmmsearch_filters", "text",
             "reporting filters deliberately opened so weak decoys are not censored"),
        ):
            if key in dec:
                out.append(_row(S, N, f"decoy_{key}", dec[key], unit, note, src))
    return out


# --------------------------------------------------------------------------------------
# stage 06 — seeds
# --------------------------------------------------------------------------------------
def _stage06(discovery: Path) -> list:
    S, N = "06", "seeds"
    out = []
    f = discovery / "seed_qc" / "seed_recovery.csv"
    d = _read_table(f)
    if not d.empty:
        src = "seed_qc/seed_recovery.csv"
        n = len(d)
        out.append(_row(S, N, "seeds_scored", n, "sequences",
                        "one row per input FASTA record", src))
        # Emitted only when the flag column is really there: `.get` returns the default for
        # an absent column, which would publish "recovery_rate_final 0.0 %" — the model lost
        # every seed — out of a measurement that was never made.
        if "before_recovered" in d.columns:
            before = int(sum(1 for v in d["before_recovered"] if str(v) == "True"))
            out.append(_row(S, N, "recovered_by_initial_model", before, "sequences",
                            "seed scores above the strict threshold against the FIRST model", src))
        if "after_recovered" in d.columns:
            after = int(sum(1 for v in d["after_recovered"] if str(v) == "True"))
            out += [
                _row(S, N, "recovered_by_final_model", after, "sequences",
                     "and against the refined model — this measures whether iterating LOST any "
                     "seed, NOT whether the seed's genomic locus was re-found by the search", src),
                _row(S, N, "recovery_rate_final", round(100.0 * after / n, 1) if n else "", "%", "", src),
            ]
        for val, cnt in _counts(d.get("status", pd.Series(dtype=str))):
            out.append(_row(S, N, f"status_{_slug(val)}", cnt, "sequences",
                            f"seed QC status '{val}'", src))
        bits = [v for v in (_num(x) for x in d.get("after_bit", [])) if v is not None]
        if bits:
            out.append(_row(S, N, "weakest_seed_bit_score", min(bits), "bits",
                            "lowest score any seed gets against the final model", src))

    # Family census is written by family_census.py; absent runs simply skip these rows.
    cen = _read_table(discovery / "family_census.csv")
    if not cen.empty and "identity_threshold" in cen.columns:
        src = "family_census.csv"
        for _, r in cen.iterrows():
            thr = str(r.get("identity_threshold", "")).strip()
            head = str(r.get("headline", "")).strip()
            note = f"clustered at {thr} identity" + (" (HEADLINE)" if head in ("True", "1", "yes") else "")
            for col, unit in (("n_seed_proteins", "proteins"), ("n_homolog_proteins", "proteins"),
                              ("n_clusters_union", "clusters"), ("n_clusters_shared", "clusters"),
                              ("n_clusters_seed_only", "clusters"), ("n_clusters_new", "clusters"),
                              ("pct_seeds_refound", "%")):
                if col in cen.columns and str(r.get(col, "")).strip():
                    out.append(_row(S, N, f"{col}_at_{_slug(thr)}", r.get(col), unit, note, src))
    return out


# --------------------------------------------------------------------------------------
# stage 07 — overprinting
# --------------------------------------------------------------------------------------
def _stage07(discovery: Path) -> list:
    S, N = "07", "overprint"
    out = []

    # overprinting_summary.csv is already written in a metric/value/note shape, so it maps
    # straight onto this schema — no reinterpretation, no chance of the two disagreeing.
    summ = _read_table(discovery / "overprinting_summary.csv")
    if not summ.empty and "metric" in summ.columns:
        for _, r in summ.iterrows():
            out.append(_row(S, N, str(r.get("metric", "")), r.get("value", ""),
                            r.get("units", ""), r.get("note", ""), "overprinting_summary.csv"))

    loci = _read_table(discovery / "overprinted_loci.csv")
    src = "overprinted_loci.csv"
    if loci.empty:
        # Fallback to the raw read-through table so the stage is still populated on a run
        # made before overprint_report existed. Same loci, fewer columns (no host gene).
        loci = _read_table(discovery / "interrupted_homologs.tsv", sep="\t")
        src = "interrupted_homologs.tsv"
    if loci.empty:
        return out

    out.append(_row(S, N, "interrupted_loci", len(loci), "loci",
                    "homolog copies carrying a premature stop, found by read-through "
                    "translation — the stop-to-stop six-frame search misses these", src))
    for val, cnt in _counts(loci.get("overprinting_support", pd.Series(dtype=str))):
        out.append(_row(S, N, f"support_{_slug(val)}", cnt, "loci",
                        f"overprinting support '{val}' — the premature stops are synonymous "
                        f"in the overlapping antisense frame", src))
    if "antisense_open_stops" in loci.columns:
        # `_num(v) == 0` on purpose, NOT `_num(v) or ...`: zero is the interesting value
        # here (a fully open antisense frame) and would be swallowed as falsy.
        openf = int(sum(1 for v in loci["antisense_open_stops"] if _num(v) == 0))
        out.append(_row(S, N, "fully_open_antisense_frame", openf, "loci",
                        "0 stops across the whole domain in the antisense frame — a real "
                        "overlapping ORF, which is what overprinting requires", src))
    for col, unit, note in (
        ("host_product", "loci", "loci with a named antisense host gene resolved"),
        ("antisense_orf_matches_host_gene", "loci",
         "the computed antisense ORF coincides with the annotated host CDS — the single "
         "strongest line of evidence"),
        ("nested_fully", "loci", "domain lies entirely inside the host gene"),
    ):
        if col not in loci.columns:
            continue
        cnt = int(sum(1 for v in loci[col] if str(v).strip() and str(v) not in ("False", "0", "nan")))
        out.append(_row(S, N, col, cnt, unit, note, src))
    lens = [v for v in (_num(x) for x in loci.get("domain_aa_len", [])) if v is not None]
    if lens:
        out.append(_row(S, N, "domain_aa_len_median", int(sorted(lens)[len(lens) // 2]), "aa",
                        "median matched domain length of the interrupted copies", src))
    return out


# --------------------------------------------------------------------------------------
# stage 08 — downstream
# --------------------------------------------------------------------------------------
def _stage08(discovery: Path) -> list:
    S, N = "08", "downstream"
    out = []

    # The alignment/tree live under downstream/tree in the run and are mirrored into
    # PACKAGE; prefer the run copy, fall back to the package copy.
    aln_dirs = [discovery / "downstream" / "tree", discovery / "PACKAGE" / "04_alignment_phylogeny"]
    stats_p = next((p / "hits.aln.stats.json" for p in aln_dirs if (p / "hits.aln.stats.json").exists()), None)
    if stats_p:
        st = _read_json(stats_p)
        src = _rel(discovery, stats_p)
        for key, unit, note in (("n_sequences", "sequences", "sequences in the alignment (homologs + seeds)"),
                                ("aln_length", "columns", "aligned length"),
                                ("gap_pct", "%", "gap fraction of the alignment"),
                                ("conserved_columns", "columns", "columns conserved across the family"),
                                ("avg_pairwise_id", "%", "mean pairwise identity")):
            if key in st:
                out.append(_row(S, N, f"alignment_{key}", st[key], unit, note, src))
        flagged = st.get("flagged_sequences")
        if isinstance(flagged, list):
            out.append(_row(S, N, "alignment_flagged_sequences", len(flagged), "sequences",
                            "sequences flagged as poorly aligned", src))

    tree_p = next((p / "hits.treefile" for p in aln_dirs if (p / "hits.treefile").exists()), None)
    if tree_p:
        try:
            txt = tree_p.read_text(encoding="utf-8", errors="replace")
            # Newick leaf labels: a name always follows '(' or ','. Internal nodes yield an
            # empty capture (the next char is '(') and are dropped, so this counts TIPS.
            tips = [t.strip() for t in re.findall(r"[(,]\s*([^(),:;]*)", txt) if t.strip()]
            out.append(_row(S, N, "tree_tips", len(tips), "tips",
                            "tips in the ML tree — homologs plus the input seeds, which are "
                            "included for context and labelled with a '_seed' suffix",
                            _rel(discovery, tree_p)))
            seeds = sum(1 for t in tips if t.lower().endswith("_seed") or t.upper().startswith("SEED"))
            if seeds:
                out.append(_row(S, N, "tree_seed_tips", seeds, "tips",
                                "input seeds included in the tree for context",
                                _rel(discovery, tree_p)))
        except OSError:
            pass

    memb = discovery / "downstream" / "clinker" / "cluster_membership.tsv"
    m = _read_table(memb, sep="\t")
    if not m.empty and "cluster_id" in m.columns:
        out.append(_row(S, N, "synteny_clusters", len({str(c) for c in m["cluster_id"] if str(c).strip()}),
                        "clusters", "gene-neighbourhood clusters compared with clinker",
                        _rel(discovery, memb)))
        out.append(_row(S, N, "clustered_hits", len(m), "hits",
                        "hits assigned to a synteny cluster", _rel(discovery, memb)))

    gbk = discovery / "downstream" / "genbank_with_sequence"
    if gbk.is_dir():
        files = sorted(gbk.glob("*.gbk"))
        if files:
            out.append(_row(S, N, "genbank_neighbourhoods", len(files), "files",
                            "real-sequence GenBank neighbourhood records — one file per "
                            "genome ACCESSION, NOT one per phage: the same phage catalogued "
                            "under a RefSeq and a GenBank accession (and under a versioned "
                            "and an unversioned id) gets a file each, so this count is "
                            "inflated by cross-database aliases. Use "
                            "genbank_neighbourhood_organisms for a phage count",
                            _rel(discovery, gbk)))
            # The honest phage count: canonical organism of each record's ORGANISM line,
            # collapsed with the same rule every other organism count in the package uses.
            orgs = set()
            unnamed = 0
            for p in files:
                name = _gbk_organism(p)
                if name:
                    orgs.add(canonical_organism(name, p.stem))
                else:
                    unnamed += 1
            orgs.discard("")
            if orgs:
                out.append(_row(S, N, "genbank_neighbourhood_organisms", len(orgs), "organisms",
                                "distinct phages among those records (canonical organism of "
                                "each record's ORGANISM line)"
                                + (f"; {unnamed} record(s) carried no ORGANISM line and are "
                                   f"not counted" if unnamed else ""),
                                _rel(discovery, gbk)))

    nb = discovery / "downstream" / "synteny" / "neighbour_gene_annotations.csv"
    nbd = _read_table(nb)
    if not nbd.empty:
        out.append(_row(S, N, "annotated_neighbour_genes", len(nbd), "genes",
                        "ordered neighbouring genes annotated around every hit", _rel(discovery, nb)))
    return out


# --------------------------------------------------------------------------------------
# public entry point
# --------------------------------------------------------------------------------------
# (code, slug, builder). The slug is both the `stage_name` value and part of the filename,
# so a reader can map a row back to its own table without a lookup.
_STAGES = (
    ("00", "input", _stage00),
    ("01", "model", _stage01),
    ("02", "search", _stage02),
    ("03", "validation", _stage03),
    ("04", "homologs", _stage04),
    ("05", "controls", _stage05),
    ("06", "seeds", _stage06),
    ("07", "overprint", _stage07),
    ("08", "downstream", _stage08),
)


def stage_rows(discovery: Path) -> dict:
    """{(code, slug): [row, ...]} for every stage. Each stage is isolated: one that raises
    (missing file, malformed CSV, a column another module renamed) contributes [] and the
    rest still export. Never raises."""
    discovery = Path(discovery)
    out = {}
    for code, slug, fn in _STAGES:
        try:
            rows = fn(discovery) or []
        except Exception:
            rows = []               # a broken stage must never cost the reader the others
        out[(code, slug)] = [r for r in rows if isinstance(r, dict)]
    return out


def build(discovery: Path, log=None) -> list:
    """Write one `stage<NN>_<slug>_summary.csv` per non-empty stage plus the concatenated
    `pipeline_stage_summary.csv`, into the run directory. Returns the paths written.

    Pure filesystem — reads files this run already produced and writes CSVs. `export_csv`
    calls it on every export (twice per run), so it must stay cheap and must not raise."""
    discovery = Path(discovery)
    written: list = []
    try:
        by_stage = stage_rows(discovery)
    except Exception as e:                     # defensive: stage_rows already swallows
        if log:
            log(f"  (stage summaries skipped: {e})")
        return written

    allrows: list = []
    for (code, slug), rows in sorted(by_stage.items()):
        path = discovery / f"stage{code}_{slug}_summary.csv"
        if not rows:
            # Remove a stale table from an earlier export in the same run. Deleting it HERE
            # is only half the job — the PACKAGE mirror in export_csv is what actually ships,
            # and it deletes the mirrored copy of anything missing from the run root (it is
            # a copy loop, so it can never remove a file on its own).
            try:
                path.unlink()
            except OSError:
                pass
            continue
        try:
            pd.DataFrame(rows).reindex(columns=list(COLUMNS)).fillna("").to_csv(path, index=False)
            written.append(str(path))
            allrows += rows
        except Exception as e:
            if log:
                log(f"  (stage {code} table skipped: {e})")

    p = discovery / SUMMARY_NAME
    if allrows:
        try:
            pd.DataFrame(allrows).reindex(columns=list(COLUMNS)).fillna("").to_csv(p, index=False)
            written.append(str(p))
        except Exception as e:
            if log:
                log(f"  ({SUMMARY_NAME} skipped: {e})")
    else:
        # Every stage came back empty. Leaving the previous export's concatenation in place
        # would publish a whole run's numbers that this run no longer stands behind — the
        # same stale-table failure as the per-stage tables above.
        try:
            p.unlink()
        except OSError:
            pass
    if log and written:
        log(f"  wrote {len(written)} stage summary table(s) ({len(allrows)} metrics)")
    return written


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discovery-dir", type=Path, required=True)
    args = ap.parse_args()
    for f in build(args.discovery_dir):
        print(f"  wrote {f}")


if __name__ == "__main__":
    main()
