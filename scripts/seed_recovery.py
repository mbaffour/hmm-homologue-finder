#!/usr/bin/env python3
"""seed_recovery.py — per-seed QC: which input seeds does the model actually recover?

The controls report an aggregate *sensitivity* (e.g. "96/101 recovered") but not
*which* seeds missed. This QC scores every input seed, by name, against:

  * the INITIAL model  — built directly from the seeds (run1 HMM)   ["before"]
  * the FINAL model    — the refined, most-complete run's HMM       ["after"]

so an outlier seed the model barely recognises — or one lost / gained as the
profile broadens over iterations — is visible by name. Writes
`seed_recovery.csv` and returns a summary. Never raises.
"""
from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path

STRICT_BITS = 45.0  # same strict bit threshold the pipeline uses to tier hits


def _seed_ids(faa: Path) -> list:
    """Seed identifiers (first token of each FASTA header), in file order."""
    ids = []
    for ln in Path(faa).read_text(errors="replace").splitlines():
        if ln.startswith(">"):
            tok = ln[1:].split()
            if tok:
                ids.append(tok[0])
    return ids


def parse_tblout_best_bits(text: str) -> dict:
    """Parse `hmmsearch --tblout` text -> {target: best full-sequence bit score}.

    Pure and testable. Full-sequence score is column 6 (0-indexed 5) of the tabular
    output; a target may appear more than once, so keep the maximum.
    """
    best: dict = {}
    for ln in text.splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        p = ln.split()
        if len(p) < 6:
            continue
        tgt = p[0]
        try:
            sc = float(p[5])
        except ValueError:
            continue
        if tgt not in best or sc > best[tgt]:
            best[tgt] = sc
    return best


def per_sequence_best_bits(hmm, faa, cpu=1) -> dict:
    """{seq_id: best full-sequence bit score} from hmmsearch of `hmm` vs `faa`.

    Returns {} on any failure (missing files, hmmsearch absent, etc.)."""
    hmm, faa = Path(hmm), Path(faa)
    if not hmm.exists() or not faa.exists():
        return {}
    with tempfile.TemporaryDirectory() as td:
        tbl = Path(td) / "seed_recovery.tbl"
        try:
            subprocess.run(
                ["hmmsearch", "--noali", "--cpu", str(cpu),
                 "--tblout", str(tbl), str(hmm), str(faa)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception:
            return {}
        try:
            return parse_tblout_best_bits(tbl.read_text(errors="replace"))
        except Exception:
            return {}


def classify(before_recovered: bool, after_recovered: bool) -> str:
    if before_recovered and after_recovered:
        return "recovered"
    if before_recovered and not after_recovered:
        return "lost_after_refinement"
    if not before_recovered and after_recovered:
        return "gained_after_refinement"
    return "never_recovered"


def seed_recovery_report(seeds_faa, before_hmm, after_hmm, out_dir,
                         cpu=1, log=print, strict: float = STRICT_BITS) -> dict:
    """Score each seed vs the initial and final models; write seed_recovery.csv.

    Returns a summary dict (empty on failure). Never raises."""
    try:
        seeds_faa = Path(seeds_faa)
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        ids = _seed_ids(seeds_faa)
        if not ids:
            return {}
        before = per_sequence_best_bits(before_hmm, seeds_faa, cpu)
        # Single-iteration runs reuse one model — don't pay for a second search.
        same = bool(before_hmm and after_hmm
                    and Path(before_hmm).resolve() == Path(after_hmm).resolve())
        after = before if same else per_sequence_best_bits(after_hmm, seeds_faa, cpu)

        rows, missing_after = [], []
        for sid in ids:
            bb, ab = before.get(sid, 0.0), after.get(sid, 0.0)
            br, ar = bb >= strict, ab >= strict
            if not ar:
                missing_after.append(sid)
            rows.append({"seed_id": sid,
                         "before_bit": round(bb, 1), "before_recovered": br,
                         "after_bit": round(ab, 1), "after_recovered": ar,
                         "status": classify(br, ar)})

        csv_path = out_dir / "seed_recovery.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["seed_id", "before_bit", "before_recovered",
                                               "after_bit", "after_recovered", "status"])
            w.writeheader()
            w.writerows(rows)

        n = len(ids)
        nb = sum(1 for r in rows if r["before_recovered"])
        na = sum(1 for r in rows if r["after_recovered"])
        summary = {"strict_bits": strict, "n_seeds": n,
                   "recovered_before": nb, "recovered_after": na,
                   "not_recovered_after": missing_after, "single_model": same,
                   "csv": str(csv_path)}
        log(f"Seed-recovery QC: {nb}/{n} seeds recovered by the initial model, "
            f"{na}/{n} by the final model (strict bit≥{strict:g}). -> {csv_path.name}")
        if missing_after:
            log("  WARNING: " + str(len(missing_after)) + " seed(s) NOT recovered by the final "
                "model: " + ", ".join(missing_after[:8])
                + (" …" if len(missing_after) > 8 else "")
                + " (likely divergent outliers; see seed_recovery.csv)")
        return summary
    except Exception as e:
        try:
            log(f"  (seed-recovery QC skipped: {e})")
        except Exception:
            pass
        return {}
