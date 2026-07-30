#!/usr/bin/env python3
"""coverage_report.py — one table stating exactly what search space was covered, and what was not.

A discovery run's tables list what WAS searched, which makes an unsearched database
indistinguishable from one that returned nothing. That is the difference between "the gene is
absent from the gut phage catalogues" and "we never looked there", and only one of those is a
finding. This consolidates every coverage-closing scan into a single statement, including the
spaces still NOT covered — because a coverage claim is only worth anything if the gaps are named.

Reads whatever exists and skips the rest, so it is useful after a partial run:
  <coverage>/seed_sources/missed_seed_scan.csv        seed source genomes
  <coverage>/catalogue_{gpd,gvd}/stream_scan_summary.json + collection_hits.tsv
  <coverage>/host_genera/results/collection_hits.tsv  prophage in host genera
  <run>/database_hit_summary.csv                      the main run's databases

Writes <coverage>/coverage_summary.csv with the shared stage schema
(space, searched, n_units, units, n_hits, verdict, note) plus a plain-language summary.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

COLS = ["space", "searched", "n_units", "units", "n_hits", "verdict", "note"]

# Spaces that remain uncovered even after a full-coverage run, with the reason. Naming these is
# the point: a reviewer asking "did you look?" gets a straight answer either way.
NOT_COVERED = [
    ("RefSeq bacterial genomes (full set)", "no", "",
     "~600 GB compressed / ~2 TB six-frame — days-to-weeks on a server. The host-genera scan is "
     "the bounded answer to the same question (a prophage copy in the phages' hosts)."),
    ("Pfam / PHROGs / VOGDB (annotation DBs)", "n/a", "",
     "profile libraries of ANNOTATED domains; an unannotated gene cannot be in them by "
     "construction, which is itself the finding for this family."),
    ("GenBank nt (whole nucleotide division)", "no", "",
     "not in the database catalog at all. Individual records ARE reachable — the seed-source "
     "scan fetches them by accession — but the division is not searched wholesale."),
]


def _rows(p: Path) -> list:
    try:
        with open(p, encoding="utf-8", errors="replace", newline="") as fh:
            return list(csv.DictReader(fh))
    except OSError:
        return []


def _tsv_n(p: Path) -> int:
    try:
        n = sum(1 for _ in open(p, encoding="utf-8", errors="replace")) - 1
        return max(0, n)
    except OSError:
        return 0


def build(run_dir: Path, cov: Path, log=print) -> dict:
    run_dir, cov = Path(run_dir), Path(cov)
    out = []

    # --- the main run's databases -------------------------------------------------------
    for r in _rows(run_dir / "database_hit_summary.csv"):
        db = r.get("database", "")
        if db.startswith("ALL"):
            continue
        hits = r.get("hits", "0")
        out.append({
            "space": db, "searched": "yes", "n_units": "", "units": "",
            "n_hits": hits, "verdict": "hits" if str(hits) not in ("0", "") else "no hits",
            "note": f"main discovery run; status={r.get('status','')}"
                    + ("" if str(hits) not in ("0", "") else
                       " — a 0 here IS a result: the gene is absent from this database"),
        })

    # --- seed source genomes -----------------------------------------------------------
    for cand in (cov / "seed_sources" / "missed_seed_scan.csv",
                 run_dir / "seed_qc" / "missed_seed_scan.csv"):
        rs = _rows(cand)
        if rs:
            found = [r for r in rs if str(r.get("verdict", "")).startswith("present")]
            nf = [r for r in rs if r.get("verdict") == "not_fetched"]
            out.append({
                "space": "seed source genomes (fetched individually)", "searched": "yes",
                "n_units": len(rs), "units": "seed genomes", "n_hits": len(found),
                "verdict": ("all present" if len(found) == len(rs)
                            else f"{len(found)}/{len(rs)} present"),
                "note": "each seed whose own genome the main search never touched, fetched by "
                        "accession and scanned. A hit here means the miss was database coverage, "
                        "not model sensitivity."
                        + (f" {len(nf)} could not be fetched." if nf else ""),
            })
            break

    # --- metagenome catalogues ---------------------------------------------------------
    for key, nice in (("gpd", "Gut Phage Database (GPD)"), ("gvd", "GVD-AVrC")):
        d = cov / f"catalogue_{key}"
        summ, prog = d / "stream_scan_summary.json", d / "_progress.json"
        st, have_summary = {}, False
        for p in (summ, prog):
            try:
                st = json.loads(p.read_text(encoding="utf-8"))
                have_summary = (p == summ)
                break
            except Exception:
                continue
        if not st:
            continue
        n = st.get("contigs_scanned", st.get("contigs_done", 0))
        h = _tsv_n(d / "collection_hits.tsv")
        # COMPLETION MUST BE POSITIVELY ASSERTED, never inferred from the absence of a flag.
        # `st.get("complete") is False` treated a MISSING key as complete — and the in-flight
        # checkpoint has no such key, it is only added by the final write. So a scan that was
        # still running, or had been killed, was published as fully searched; at zero hits that
        # became "the family is absent from this catalogue". stream_scan_summary.json is only
        # written when a scan ends, so its absence alone means the scan did not finish.
        partial = (bool(st.get("bounded_test"))
                   or st.get("complete") is not True
                   or not have_summary
                   or bool(st.get("ended_early"))
                   or int(st.get("failed_batches") or 0) > 0)
        # A hits table that is missing or empty while the checkpoint recorded hits means the
        # table is the unreliable witness — never silently prefer the 0.
        h_ckpt = int(st.get("hits") or 0)
        if h_ckpt and not h:
            h = h_ckpt
        out.append({
            "space": nice, "searched": "partial" if partial else "yes",
            "n_units": n, "units": "contigs", "n_hits": h,
            "verdict": ("hits (lower bound)" if (h and partial) else
                        "hits" if h else ("no hits — INCOMPLETE" if partial else "no hits")),
            "note": (("INCOMPLETE — only %s contig(s) were scanned" % f"{n:,}")
                     + (" (bounded by --max-contigs)" if st.get("bounded_test") else "")
                     + (" ; the stream ENDED EARLY" if st.get("ended_early") else "")
                     + (" ; %d batch scan(s) FAILED" % int(st.get("failed_batches") or 0)
                        if st.get("failed_batches") else "")
                     + (" ; no end-of-scan summary was written, so the scan did not finish"
                        if not have_summary else "")
                     + ". A hit count here is a LOWER BOUND and a 0 is NOT evidence of absence."
                     ) if partial else
                    ("streamed, batched and discarded; nothing cached. 0 hits = the family is "
                     "absent from this catalogue" if not h else
                     "streamed, batched and discarded; nothing cached"),
        })

    # --- host genera (prophage) --------------------------------------------------------
    for cand in (cov / "host_genera" / "results" / "collection_hits.tsv",
                 run_dir / "host_genera_scan" / "results" / "collection_hits.tsv"):
        if cand.exists():
            h = _tsv_n(cand)
            urls = cand.parent.parent / "host_genome_urls.txt"
            ng = 0
            try:
                ng = sum(1 for ln in open(urls, encoding="utf-8") if ln.strip())
            except OSError:
                pass
            out.append({
                "space": "host genera — RefSeq representative genomes", "searched": "yes",
                "n_units": ng or "", "units": "bacterial genomes", "n_hits": h,
                "verdict": "hits" if h else "no hits",
                "note": "where a phage gene would sit as a prophage. 0 hits = the family is "
                        "phage-specific in these hosts" if not h else "prophage copies found",
            })
            break

    for space, searched, n, note in NOT_COVERED:
        out.append({"space": space, "searched": searched, "n_units": n, "units": "",
                    "n_hits": "", "verdict": "not searched", "note": note})

    cov.mkdir(parents=True, exist_ok=True)
    dst = cov / "coverage_summary.csv"
    with dst.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(out)

    # ALSO write it into the RUN directory. This table is the whole point of the exercise — it
    # is what makes a coverage claim falsifiable — and it was landing only in the scan's own
    # output folder, outside the run and outside the shareable package. export_csv's
    # TABLE_EXPORTS then mirrors it into PACKAGE/01_summary_tables.
    try:
        run_copy = Path(run_dir) / "coverage_summary.csv"
        if run_copy.resolve() != dst.resolve():
            run_copy.write_text(dst.read_text(encoding="utf-8"), encoding="utf-8")
            log(f"  -> {run_copy}")
    except OSError as e:
        log(f"  (could not copy the coverage summary into the run dir: {e})")

    n_yes = sum(1 for r in out if r["searched"] == "yes")
    n_part = sum(1 for r in out if r["searched"] == "partial")
    n_no = sum(1 for r in out if r["searched"] == "no")
    log("")
    log("  COVERAGE SUMMARY")
    for r in out:
        flag = {"yes": "[searched]", "partial": "[PARTIAL] ", "no": "[not done]",
                "n/a": "[n/a]     "}.get(r["searched"], "[?]")
        u = f"{r['n_units']} {r['units']}".strip()
        log(f"    {flag} {str(r['space'])[:44]:<45} {u:<22} hits={r['n_hits']}")
    log(f"  -> {dst}")
    log(f"     {n_yes} space(s) fully searched, {n_part} partial, {n_no} not searched "
        f"(the not-searched ones are named in the table, with why)")
    return {"csv": str(dst), "rows": len(out), "searched": n_yes,
            "partial": n_part, "not_searched": n_no}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--coverage-dir", required=True, type=Path)
    a = ap.parse_args()
    build(a.run_dir, a.coverage_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
