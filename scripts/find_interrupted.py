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


def _run(genomes: Path, hmm: Path, out: Path, min_bit: float, cpu: int, log=print) -> dict:
    tmp = Path(tempfile.mkdtemp(prefix="interrupted_"))
    search_faa = tmp / "readthrough_search.faa"
    marker_faa = tmp / "readthrough_marker.faa"
    n_frames = 0
    # Pass 1: stream contigs -> read-through frames (search + marker), low memory.
    with open_maybe_gz(genomes) as fh, open(search_faa, "w") as sf, open(marker_faa, "w") as mf:
        for rec in SeqIO.parse(fh, "fasta"):
            cid = rec.id
            for strand, frame, search, marker in _frames(str(rec.seq)):
                name = f"{cid}__{strand}{frame}"
                sf.write(f">{name}\n{search}\n")
                mf.write(f">{name}\n{marker}\n")
                n_frames += 1
    log(f"  read-through: {n_frames} frames written; searching with the family HMM…")
    # hmmsearch the family HMM against the read-through frames.
    domtbl = tmp / "hits.domtbl"
    try:
        subprocess.run(["hmmsearch", "--noali", "--cpu", str(cpu), "--domtblout", str(domtbl),
                        str(hmm), str(search_faa)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        log(f"  (find-interrupted: hmmsearch failed: {e})")
        return {"error": str(e)}
    # Index markers so we only load the (few) sequences that actually hit.
    marker_idx = SeqIO.index(str(marker_faa), "fasta")
    rows, n_total, n_interrupted = [], 0, 0
    for ln in domtbl.read_text(errors="replace").splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        if len(p) < 23:
            continue
        try:
            name = p[0]
            dom_bits = float(p[13])
            i_eval = float(p[12])
            env_from, env_to = int(p[19]), int(p[20])
        except (ValueError, IndexError):
            continue
        if dom_bits < min_bit:
            continue
        n_total += 1
        marker = str(marker_idx[name].seq) if name in marker_idx else ""
        n_stops, positions = count_envelope_stops(marker, env_from, env_to)
        if n_stops < 1:
            continue
        n_interrupted += 1
        contig, sf_tag = name.rsplit("__", 1)
        strand, frame = sf_tag[0], sf_tag[1]
        dom_aa = marker[env_from - 1:env_to]
        rows.append({
            "contig": contig, "strand": strand, "frame": frame,
            "env_from_aa": env_from, "env_to_aa": env_to,
            "domain_aa_len": env_to - env_from + 1,
            "internal_stops": n_stops,
            "stop_aa_positions": ";".join(str(x) for x in positions),
            "domain_bit_score": round(dom_bits, 1),
            "i_evalue": f"{i_eval:.2g}",
            "domain_aa_with_stops": dom_aa,
        })
    marker_idx.close()
    rows.sort(key=lambda r: -r["domain_bit_score"])
    import csv
    cols = ["contig", "strand", "frame", "env_from_aa", "env_to_aa", "domain_aa_len",
            "internal_stops", "stop_aa_positions", "domain_bit_score", "i_evalue",
            "domain_aa_with_stops"]
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
    summary = {"matches_scored": n_total, "interrupted_candidates": n_interrupted, "out": str(out)}
    log(f"  find-interrupted: {n_total} read-through matches ≥{min_bit} bits; "
        f"{n_interrupted} carry an internal stop (interrupted candidates) -> {out.name}")
    return summary


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
