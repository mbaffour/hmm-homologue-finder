#!/usr/bin/env python3
"""find_interrupted.py — find homologs interrupted by a premature stop codon.

WHY THIS EXISTS
---------------
The standard search builds **stop-to-stop** six-frame ORFs, so a homolog whose
gene carries an internal stop codon gets its ORF truncated *at* the stop and is
split or missed entirely. That is exactly what happens to an **overprinted** gene:
when a small gene is encoded antisense to (or in a different frame from) a large
essential gene — e.g. gp75 overprinted on a virion DNA-directed RNA polymerase —
a point mutation can be a **premature stop in the small-gene frame yet silent
(synonymous) in the RNA-polymerase frame**. Selection keeps the polymerase intact
and tolerates the truncation, so the small gene is interrupted but invisible to a
stop-to-stop search.

WHAT IT DOES
------------
Re-translates every reading frame with **read-through** (stop codons are kept, not
broken on), so the full domain — including any internal stops — is searchable with
the family HMM. Every HMM match whose domain envelope contains >= 1 internal stop
is reported as a candidate **interrupted / overprinted** homolog, with the stop
positions, so you can study the truncation (and, with the overlapping frame,
whether each stop is silent there).

  python3 find_interrupted.py --genomes genomes.fa[.gz] --hmm profile.hmm \
      --out interrupted_homologs.tsv [--min-bit 25] [--cpu 8]
"""
from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

STOP_SEARCH = "X"   # HMMER-friendly placeholder for a stop in the SEARCH sequence
MIN_AA = 30         # ignore frames shorter than this (matches the six-frame floor)


def open_maybe_gz(p: Path):
    p = Path(p)
    return gzip.open(p, "rt") if p.suffix == ".gz" else open(p)


def read_through_aa(nt: str, strand: str, frame: int) -> tuple[str, str]:
    """Translate one reading frame with read-through. Returns (search_aa, marker_aa):
    search_aa has 'X' at stops (so the HMM matches across them); marker_aa keeps '*'
    at the same positions so the stops can be located afterwards."""
    s = nt.upper().replace("U", "T")
    if strand == "-":
        s = str(Seq(s).reverse_complement())
    s = s[frame:]
    s = s[: len(s) // 3 * 3]
    marker = str(Seq(s).translate()) if s else ""   # '*' at stops
    return marker.replace("*", STOP_SEARCH), marker


def count_envelope_stops(marker_aa: str, env_from: int, env_to: int) -> tuple[int, list]:
    """Internal stops within the 1-based inclusive envelope [env_from, env_to].
    Returns (count, [1-based positions]). The terminal position is excluded so a
    stop that merely ends the envelope is not counted as 'internal'."""
    seg = marker_aa[env_from - 1:env_to]
    if seg.endswith("*"):
        seg = seg[:-1]
    pos = [env_from + i for i, c in enumerate(seg) if c == "*"]
    return len(pos), pos


def _frames(seq: str):
    """Yield (strand, frame, search_aa, marker_aa) for all six frames."""
    for strand in ("+", "-"):
        for frame in (0, 1, 2):
            search, marker = read_through_aa(seq, strand, frame)
            if len(search) >= MIN_AA:
                yield strand, frame, search, marker


BATCH_CONTIGS = 500   # search this many contigs per hmmsearch call (bounds temp + memory)
WIN_AA = 5000         # window long read-through frames: HMMER aborts (SIGABRT) on whole-
WIN_STEP = 4400       # genome-length sequences (10k-100k+ aa). Overlap = WIN_AA - WIN_STEP
                      # (>> any domain), so every domain lies fully inside >=1 window.
ROW_COLS = ["contig", "strand", "frame", "env_from_aa", "env_to_aa", "domain_aa_len",
            "internal_stops", "stop_aa_positions", "domain_bit_score", "i_evalue",
            "domain_aa_with_stops"]


def _windows(n: int):
    """Yield (offset, length) windows covering a length-n frame with WIN overlap."""
    if n <= WIN_AA:
        yield 0, n
        return
    off = 0
    while off < n:
        yield off, min(WIN_AA, n - off)
        if off + WIN_AA >= n:
            break
        off += WIN_STEP


def _search_batch(frames: list, markers: dict, hmm: Path, workdir: Path,
                  min_bit: float, cpu: int) -> tuple[int, list]:
    """hmmsearch one batch of (name, search_aa); return (n_scored, interrupted_rows).
    markers maps name -> marker_aa (with '*') for the SAME batch (kept small)."""
    import shutil
    sfa = workdir / "batch.faa"
    with open(sfa, "w") as f:
        for name, search in frames:
            f.write(f">{name}\n{search}\n")
    dt = workdir / "batch.domtbl"
    subprocess.run(["hmmsearch", "--noali", "--cpu", str(cpu), "--domtblout", str(dt),
                    str(hmm), str(sfa)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    n_scored, rows = 0, []
    for ln in dt.read_text(errors="replace").splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        if len(p) < 23:
            continue
        try:
            name, dom_bits, i_eval = p[0], float(p[13]), float(p[12])
            env_from, env_to = int(p[19]), int(p[20])
        except (ValueError, IndexError):
            continue
        if dom_bits < min_bit:
            continue
        n_scored += 1
        marker = markers.get(name, "")
        n_stops, positions = count_envelope_stops(marker, env_from, env_to)
        if n_stops < 1:
            continue
        contig, sf_tag, off_s = name.rsplit("__", 2)   # name = contig__SF__windowOffset
        offset = int(off_s)
        rows.append({
            "contig": contig, "strand": sf_tag[0], "frame": sf_tag[1],
            "env_from_aa": offset + env_from, "env_to_aa": offset + env_to,
            "domain_aa_len": env_to - env_from + 1,
            "internal_stops": n_stops,
            "stop_aa_positions": ";".join(str(offset + x) for x in positions),
            "domain_bit_score": round(dom_bits, 1),
            "i_evalue": f"{i_eval:.2g}",
            "domain_aa_with_stops": marker[env_from - 1:env_to],
        })
    sfa.unlink(missing_ok=True)
    dt.unlink(missing_ok=True)
    return n_scored, rows


def _run(genomes: Path, hmm: Path, out: Path, min_bit: float, cpu: int, log=print) -> dict:
    """Batched read-through scan — streams the DB in BATCH_CONTIGS-sized chunks so a
    huge nucleotide DB never produces one giant temp file (which abort hmmsearch).
    Temp lives next to the output (spacious), not /tmp."""
    import csv
    import shutil
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    base = out.parent if out.parent.exists() else Path(".")
    workdir = Path(tempfile.mkdtemp(prefix="interrupted_", dir=str(base)))
    all_rows, n_scored, n_batches, n_frames = [], 0, 0, 0
    frames, markers, nctg = [], {}, 0
    try:
        with open_maybe_gz(genomes) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                for strand, frame, search, marker in _frames(str(rec.seq)):
                    for off, wlen in _windows(len(search)):
                        sw = search[off:off + wlen]
                        if len(sw) < MIN_AA:
                            continue
                        name = f"{rec.id}__{strand}{frame}__{off}"
                        frames.append((name, sw))
                        markers[name] = marker[off:off + wlen]
                        n_frames += 1
                nctg += 1
                if nctg >= BATCH_CONTIGS:
                    try:
                        ns, rows = _search_batch(frames, markers, hmm, workdir, min_bit, cpu)
                        n_scored += ns
                        all_rows += rows
                    except Exception as e:
                        log(f"  (find-interrupted: a batch failed, skipped: {e})")
                    n_batches += 1
                    frames, markers, nctg = [], {}, 0
            if frames:
                try:
                    ns, rows = _search_batch(frames, markers, hmm, workdir, min_bit, cpu)
                    n_scored += ns
                    all_rows += rows
                    n_batches += 1
                except Exception as e:
                    log(f"  (find-interrupted: final batch failed, skipped: {e})")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    log(f"  read-through: {n_frames} windowed frames in {n_batches} batch(es); "
        f"{n_scored} matches ≥{min_bit} bits scored")
    # Overlapping windows can report the same interrupted domain twice — dedup by
    # (contig, strand, frame, stop positions), keeping the best-scoring copy.
    best: dict = {}
    for r in all_rows:
        key = (r["contig"], r["strand"], r["frame"], r["stop_aa_positions"])
        if key not in best or r["domain_bit_score"] > best[key]["domain_bit_score"]:
            best[key] = r
    all_rows = list(best.values())
    all_rows.sort(key=lambda r: -r["domain_bit_score"])
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ROW_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)
    log(f"  find-interrupted: {len(all_rows)} match(es) carry an internal stop "
        f"(interrupted candidates) -> {out.name}")
    return {"matches_scored": n_scored, "interrupted_candidates": len(all_rows), "out": str(out)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", type=Path, required=True,
                    help="nucleotide FASTA (.fa/.fna, optionally .gz) to scan")
    ap.add_argument("--hmm", type=Path, required=True, help="the family profile HMM")
    ap.add_argument("--out", type=Path, default=Path("interrupted_homologs.tsv"))
    ap.add_argument("--min-bit", type=float, default=25.0,
                    help="minimum domain bit score to report (default 25)")
    ap.add_argument("--cpu", type=int, default=8)
    args = ap.parse_args()
    if not args.genomes.exists():
        sys.exit(f"genomes FASTA not found: {args.genomes}")
    if not args.hmm.exists():
        sys.exit(f"HMM not found: {args.hmm}")
    s = _run(args.genomes, args.hmm, args.out, args.min_bit, args.cpu)
    if "error" in s:
        sys.exit(s["error"])
    print(f"  matches scored: {s['matches_scored']}; interrupted candidates: "
          f"{s['interrupted_candidates']}; written to {s['out']}")


if __name__ == "__main__":
    main()
