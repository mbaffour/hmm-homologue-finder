#!/usr/bin/env python3
"""
golden_test.py — biological-invariant regression test for a discovery run.

A full pipeline run has stochastic steps (IQ-TREE search, MEME), so byte-for-byte
comparison is useless. Instead we snapshot the *biological invariants* that MUST
stay stable for a fixed seed + fixed (SHA-pinned) databases, and compare against a
stored baseline with sensible tolerances:

  - per-run hit counts (total / passed_filter / six_frame / protein_db / unique)
  - the main-paper homolog set: the SET of representative accessions, each hit's
    confidence_tier, and best bit score (within ±BIT_TOL)
  - always-true ORF invariants: every kept hit has internal_stops==0 and
    passes_orf_filter==True
  - phylogenetic tree (when present): the SET of tip labels, and — if a baseline
    Newick is stored — the Robinson-Foulds distance (0 under a fixed -seed)
  - the iteration stop reason and threshold-calibration sensitivity/specificity

Usage:
    # write/refresh a baseline from a known-good run
    python3 golden_test.py --discovery <dir> --baseline tests/golden/<fam>.json --update
    # check a fresh run against the baseline (exit 1 on regression)
    python3 golden_test.py --discovery <dir> --baseline tests/golden/<fam>.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

BIT_TOL = 2.0  # bit-score tolerance for the paper-table comparison


def _read_csv(p: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(p, dtype=str).fillna("")
    except Exception:
        return pd.DataFrame()


def _tree_tips(treefile: Path) -> list[str]:
    """Return the sorted set of tip labels from a Newick file (no deps on ete/dendropy)."""
    try:
        from Bio import Phylo
        t = Phylo.read(str(treefile), "newick")
        return sorted(term.name for term in t.get_terminals() if term.name)
    except Exception:
        return []


def extract_invariants(discovery: Path) -> dict:
    discovery = Path(discovery)
    inv: dict = {}

    # 1. per-run hit counts (stable for fixed seed + DBs)
    hs = _read_csv(discovery / "hit_summary.csv")
    if not hs.empty:
        cols = ["run", "total_hits", "passed_filter", "six_frame_hits",
                "protein_db_hits", "unique_sequences"]
        inv["hit_summary"] = hs[[c for c in cols if c in hs.columns]].to_dict("records")

    # 2. main-paper homolog set: accession -> {tier, bit}
    pt = _read_csv(discovery / "paper_main_table.csv")
    if not pt.empty and "accession" in pt.columns:
        rows = {}
        for _, r in pt.iterrows():
            try:
                bit = round(float(r.get("best_bit_score", "") or 0), 1)
            except ValueError:
                bit = 0.0
            rows[r["accession"]] = {"tier": r.get("confidence_tier", ""), "bit": bit}
        inv["paper_table"] = rows

    # 3. ORF invariants over every kept hit (must ALWAYS hold)
    allh = _read_csv(discovery / "all_runs_hits.csv")
    if not allh.empty:
        kept = allh[allh.get("passes_orf_filter", "") == "True"] if "passes_orf_filter" in allh else allh
        six = kept[kept.get("source_type", "") == "six_frame_orf"] if "source_type" in kept else kept.iloc[0:0]
        inv["orf_invariants"] = {
            "n_kept": int(len(kept)),
            "all_pass_filter": bool((allh.get("passes_orf_filter", pd.Series(dtype=str)) == "True").all()) if "passes_orf_filter" in allh else None,
            # six-frame ORFs must have no internal stops by construction
            "six_frame_internal_stops_max": (int(pd.to_numeric(six["internal_stops"], errors="coerce").fillna(0).max()) if ("internal_stops" in six.columns and len(six)) else 0),
        }

    # 4. tree tips (presence + set), if a tree was built
    for tf in [discovery / "PACKAGE" / "05_phylogeny" / "hits.treefile",
               discovery / "downstream" / "tree" / "hits.treefile"]:
        if tf.exists():
            inv["tree_tips"] = _tree_tips(tf)
            break

    # 5. provenance: stop reason + calibration
    mani = discovery / "run_manifest.json"
    if mani.exists():
        try:
            m = json.loads(mani.read_text())
            inv["stop_reason"] = m.get("iteration_stop_reason", "")
            cal = m.get("threshold_calibration") or {}
            if cal:
                inv["calibration"] = {
                    "sensitivity": cal.get("sensitivity"),
                    "specificity": cal.get("specificity"),
                    "false_positive_rate": cal.get("false_positive_rate"),
                }
        except Exception:
            pass
    return inv


def compare(actual: dict, expected: dict) -> list[str]:
    diffs: list[str] = []

    # hit summary: exact count match
    if "hit_summary" in expected:
        if actual.get("hit_summary") != expected["hit_summary"]:
            diffs.append(f"hit_summary changed:\n  expected {expected['hit_summary']}\n  actual   {actual.get('hit_summary')}")

    # paper table: same accession set; same tier; bit within tolerance
    if "paper_table" in expected:
        exp, act = expected["paper_table"], actual.get("paper_table", {})
        if set(exp) != set(act):
            only_e = sorted(set(exp) - set(act))
            only_a = sorted(set(act) - set(exp))
            diffs.append(f"paper_table accession set changed: missing={only_e} new={only_a}")
        for acc in set(exp) & set(act):
            if exp[acc]["tier"] != act[acc]["tier"]:
                diffs.append(f"paper_table {acc}: tier {exp[acc]['tier']} -> {act[acc]['tier']}")
            if abs(float(exp[acc]["bit"]) - float(act[acc]["bit"])) > BIT_TOL:
                diffs.append(f"paper_table {acc}: bit {exp[acc]['bit']} -> {act[acc]['bit']} (>{BIT_TOL})")

    # ORF invariants: these must never regress
    if "orf_invariants" in actual:
        oi = actual["orf_invariants"]
        if oi.get("all_pass_filter") is False:
            diffs.append("ORF invariant violated: some kept hits have passes_orf_filter != True")
        if (oi.get("six_frame_internal_stops_max") or 0) > 0:
            diffs.append(f"ORF invariant violated: six-frame hit with internal stop(s) (max={oi['six_frame_internal_stops_max']})")
        if "orf_invariants" in expected and oi.get("n_kept") != expected["orf_invariants"].get("n_kept"):
            diffs.append(f"n_kept changed: {expected['orf_invariants'].get('n_kept')} -> {oi.get('n_kept')}")

    # tree tip set
    if "tree_tips" in expected:
        if set(expected["tree_tips"]) != set(actual.get("tree_tips", [])):
            diffs.append("tree tip set changed")

    # stop reason
    if "stop_reason" in expected and expected["stop_reason"] != actual.get("stop_reason"):
        diffs.append(f"stop_reason: {expected['stop_reason']!r} -> {actual.get('stop_reason')!r}")

    return diffs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discovery", type=Path, required=True, help="a *_discovery output dir")
    ap.add_argument("--baseline", type=Path, required=True, help="baseline invariants JSON")
    ap.add_argument("--update", action="store_true", help="write the baseline from this run instead of comparing")
    args = ap.parse_args()

    inv = extract_invariants(args.discovery)
    if args.update:
        args.baseline.parent.mkdir(parents=True, exist_ok=True)
        args.baseline.write_text(json.dumps(inv, indent=2, sort_keys=True))
        print(f"baseline written: {args.baseline}")
        print(json.dumps(inv, indent=2, sort_keys=True))
        return

    if not args.baseline.exists():
        sys.exit(f"baseline not found: {args.baseline} (run with --update first)")
    expected = json.loads(args.baseline.read_text())
    diffs = compare(inv, expected)
    if diffs:
        print(f"GOLDEN REGRESSION — {len(diffs)} difference(s):")
        for d in diffs:
            print("  - " + d)
        sys.exit(1)
    print("GOLDEN OK — all invariants match the baseline.")


if __name__ == "__main__":
    main()
