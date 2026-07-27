#!/usr/bin/env python3
"""Empirical false-discovery control against the search space that actually produced the
hits: the six-frame ORFs of the genome databases, reversed.

Why this exists. The standing control suite scores the model against curated protein sets
(fungal, mammalian, archaeal) and shuffled seeds. None of those is a six-frame translation
of genomic DNA — yet for an unannotated gene essentially every hit comes from a six-frame
ORF. "Specificity 1.0 against unrelated proteomes" therefore only says those negatives were
easy; it does not bound the false-positive rate of the search that was actually run. A
reviewer asking "how do you know these aren't spurious six-frame ORFs?" is not answered by
it. This is.

The decoy reverses each ORF's amino-acid string. Reversal preserves length and composition
*exactly*, and unlike a shuffle it also preserves local composition runs — it destroys only
the order of the motifs the model recognises. A score reached on this set is what the model
attains on realistic non-homologous sequence of the same shape.

Reporting filters are opened (``--max -E 100000``): HMMER's defaults censor precisely the
weak decoy scores this measurement needs, and with them the best decoy is unmeasurable and
the gap cannot be stated. E-values use ``-Z`` = the full ORF count, so they stay on the
scale of the real search even when the decoy is a sample of it.

Outputs ``sixframe_decoy_control.json``: best decoy score, weakest true positive, the gap
between them, and the empirical FDR at the operating threshold (decoy hits scaled back up
to the full search space).

Standalone:
    python3 sixframe_decoy_control.py --hmm profile.hmm --hits-tsv validated/hits.tsv \\
        --cache-dir ~/.cache/hmm-homologue-finder --out controls/
"""

from __future__ import annotations

import argparse
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

MIN_KEEP_AA = 30          # matches the engine's six-frame min ORF length
DEFAULT_SAMPLE = 200_000  # ORFs; enough to bound the tail, ~1 min of hmmsearch


def find_sixframe_files(cache_dir: Path, db_names=None) -> list:
    """Cached six-frame ORF FASTAs, i.e. the exact sequences that were searched.

    The engine writes ``<source>.sixframe.min<N>.faa`` beside each downloaded nucleotide
    database, so the decoy can be built from the real search space instead of a stand-in.
    """
    cache_dir = Path(cache_dir)
    roots = [cache_dir / "cache", cache_dir]
    files = []
    for root in roots:
        if root.is_dir():
            files += sorted(root.glob("*/*.sixframe.min*.faa"))
    if db_names:
        want = {str(d).lower().replace(" ", "_") for d in db_names}
        sel = [f for f in files if any(w in f.parent.name.lower() for w in want)]
        files = sel or files
    return files


def _iter_fasta(path: Path):
    """Yield amino-acid sequences from a FASTA, streaming (files here are many GB)."""
    seq = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if seq:
                    yield "".join(seq)
                    seq = []
            else:
                seq.append(ln.strip())
    if seq:
        yield "".join(seq)


def build_decoy(sixframe_files, out_faa: Path, n: int = DEFAULT_SAMPLE,
                seed: int = 42, log=print) -> tuple:
    """Reservoir-sample ORFs from the real six-frame set and write them reversed.

    Returns (out_faa, n_sampled, n_total). One pass gives both the sample and the true
    search-space size needed for ``-Z`` and for scaling the FDR back up.
    """
    rng = random.Random(seed)
    reservoir, total = [], 0
    for f in sixframe_files:
        for s in _iter_fasta(f):
            if len(s) < MIN_KEEP_AA:
                continue
            total += 1
            if len(reservoir) < n:
                reservoir.append(s)
            else:
                j = rng.randrange(total)
                if j < n:
                    reservoir[j] = s
    out_faa = Path(out_faa)
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    with out_faa.open("w", encoding="utf-8") as fh:
        for i, s in enumerate(reservoir, 1):
            # reversal: identical length and composition, motif order destroyed
            fh.write(f">decoy_{i:07d} reversed_sixframe_orf len={len(s)}\n{s[::-1]}\n")
    log(f"  decoy: {len(reservoir):,} reversed ORFs sampled from {total:,} real six-frame ORFs")
    return out_faa, len(reservoir), total


def _hmmsearch_exe() -> str:
    """Locate hmmsearch, falling back to the interpreter's own env/bin.

    Not finding it must never be mistaken for "the decoys scored nothing" — that would
    turn a broken measurement into a clean bill of health.
    """
    exe = shutil.which("hmmsearch")
    if exe:
        return exe
    cand = Path(sys.executable).parent / "hmmsearch"
    return str(cand) if cand.exists() else ""


def _scan(hmm: Path, faa: Path, tbl: Path, z_total: int, cpu: int = 4):
    """hmmsearch with reporting filters OPEN. Returns the decoy bit scores, or None if
    the scan could not be performed.

    ``--max -E 100000`` matters: with HMMER's defaults the weak scores that define the
    decoy ceiling are never reported, so the gap would be unmeasurable.

    The None-vs-[] distinction is load-bearing. An earlier version returned [] when
    hmmsearch was missing, so `max(scores)` fell back to 0.0 and the control reported
    "best decoy 0.0 bits, clean separation, FDR 0.0" — a confident result from a scan
    that never ran.
    """
    exe = _hmmsearch_exe()
    if not exe:
        return None
    cmd = [exe, "--max", "-E", "100000", "--cpu", str(cpu), "--tblout", str(tbl)]
    if z_total > 0:
        cmd += ["-Z", str(z_total)]
    cmd += [str(hmm), str(faa)]
    try:
        proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=False)
    except OSError:
        return None
    if proc.returncode != 0 or not tbl.exists():
        return None
    scores = []
    for ln in tbl.read_text(errors="replace").splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        if len(f) > 5:
            try:
                scores.append(float(f[5]))
            except ValueError:
                pass
    return scores


def _true_scores(hits_tsv: Path) -> list:
    """Bit scores of the validated hits (the true positives being defended)."""
    try:
        lines = Path(hits_tsv).read_text(errors="replace").splitlines()
    except OSError:
        return []
    if not lines:
        return []
    hdr = lines[0].split("\t")
    try:
        bi = hdr.index("bit_score")
    except ValueError:
        return []
    out = []
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) > bi:
            try:
                out.append(float(f[bi]))
            except ValueError:
                pass
    return out


def run(hmm: Path, hits_tsv: Path, cache_dir: Path, out_dir: Path, threshold: float = 45.0,
        n: int = DEFAULT_SAMPLE, cpu: int = 4, db_names=None, log=print) -> dict:
    """Build the decoy, scan it, and report the gap + empirical FDR. Never raises."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    res = {"status": "skipped", "threshold": threshold}
    try:
        files = find_sixframe_files(cache_dir, db_names)
        if not files:
            res["reason"] = "no cached six-frame ORF files found"
            log("  (six-frame decoy control skipped: no cached six-frame ORFs)")
            return res
        faa, n_sampled, n_total = build_decoy(files, out_dir / "sixframe_decoy.faa",
                                              n=n, log=log)
        decoy = _scan(hmm, faa, out_dir / "sixframe_decoy.tbl", n_total, cpu)
        if decoy is None:                       # scan could not run — report that, do NOT
            res.update({"status": "error",      # let a missing tool look like a clean gap
                        "reason": "hmmsearch unavailable or failed; no FDR measured"})
            (out_dir / "sixframe_decoy_control.json").write_text(json.dumps(res, indent=2))
            log("  (six-frame decoy control: hmmsearch unavailable — FDR NOT measured)")
            return res
        true = [s for s in _true_scores(hits_tsv)]
        frac = (n_sampled / n_total) if n_total else 1.0
        best_decoy = max(decoy) if decoy else None    # None = scan ran, nothing scored
        tp = [s for s in true if s >= threshold]
        weakest_tp = min(tp) if tp else None
        n_decoy_over = sum(1 for s in decoy if s >= threshold)
        expected_full = (n_decoy_over / frac) if frac else 0.0
        res = {
            "status": "ok",
            "decoy": "reversed six-frame ORFs from the searched genome databases",
            "n_decoy_sequences": n_sampled,
            "n_sixframe_orfs_total": n_total,
            "sampled_fraction": round(frac, 6),
            "threshold": threshold,
            "n_decoys_scored": len(decoy),
            "best_decoy_bit_score": (None if best_decoy is None else round(best_decoy, 1)),
            "weakest_true_positive_bit_score": (None if weakest_tp is None
                                                else round(weakest_tp, 1)),
            "gap_bits": (None if (weakest_tp is None or best_decoy is None)
                         else round(weakest_tp - best_decoy, 1)),
            "clean_separation": (weakest_tp is not None and best_decoy is not None
                                 and weakest_tp > best_decoy),
            "decoys_at_or_above_threshold": n_decoy_over,
            "expected_decoys_in_full_search_space": round(expected_full, 1),
            "true_positives_at_or_above_threshold": len(tp),
            "empirical_fdr": (round(expected_full / (expected_full + len(tp)), 6)
                              if (expected_full + len(tp)) > 0 else 0.0),
            "hmmsearch_filters": "--max -E 100000 (reporting filters open so weak decoys are not censored)",
        }
        (out_dir / "sixframe_decoy_control.json").write_text(json.dumps(res, indent=2))
        if best_decoy is None:
            log(f"  Six-frame decoy FDR: scanned {n_sampled:,} reversed ORFs, none scored at all "
                f"(FDR {res['empirical_fdr']})")
        elif res["clean_separation"]:
            log(f"  Six-frame decoy FDR: best decoy {res['best_decoy_bit_score']} bits vs "
                f"weakest true positive {res['weakest_true_positive_bit_score']} bits "
                f"(gap {res['gap_bits']}); empirical FDR {res['empirical_fdr']}")
        else:
            log(f"  Six-frame decoy FDR: NO clean separation — best decoy "
                f"{res['best_decoy_bit_score']} bits reaches the true-positive range; "
                f"empirical FDR {res['empirical_fdr']}")
    except Exception as e:                                    # never fail the run
        res = {"status": "error", "error": str(e)[:200], "threshold": threshold}
        log(f"  (six-frame decoy control failed: {e})")
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--hits-tsv", required=True)
    ap.add_argument("--cache-dir", default=str(Path.home() / ".cache" / "hmm-homologue-finder"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--threshold", type=float, default=45.0)
    ap.add_argument("--n", type=int, default=DEFAULT_SAMPLE)
    ap.add_argument("--cpu", type=int, default=4)
    a = ap.parse_args()
    r = run(Path(a.hmm), Path(a.hits_tsv), Path(a.cache_dir), Path(a.out),
            threshold=a.threshold, n=a.n, cpu=a.cpu)
    print(json.dumps(r, indent=2))
    return 0 if r.get("status") in ("ok", "skipped") else 1


if __name__ == "__main__":
    sys.exit(main())
