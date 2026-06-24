#!/usr/bin/env python3
"""scan_genome.py — does ONE genome contain my gene?

A focused, single-genome counterpart to the full discovery pipeline. Build (or
accept) a profile HMM for your gene of interest and scan **one** genome for it,
reporting a per-genome result: present / not detected, the hit's genome
coordinates, ORF validation (start M, ends at a stop, internal stops), the HMM
score, and the hit sequence (protein + DNA). Because the scan is read-through
(stop codons kept, not broken on), a clean copy is found with 0 internal stops
**and** a stop-interrupted / overprinted copy is found too (reported only with
``--find-interrupted``), with the same overprinting (silent-stop) test the
discovery pipeline uses.

    # from seed sequences (builds the HMM for you; nucleotide seeds are translated):
    python3 scan_genome.py --seeds gene_seeds.faa --genome genome.fna --out scan_out

    # from an existing profile HMM:
    python3 scan_genome.py --hmm gene.hmm --genome genome.fna --out scan_out

    # fetch the genome from NCBI by accession instead of a local file:
    python3 scan_genome.py --hmm gene.hmm --accession KX098390 --email you@inst.edu

    # also report stop-interrupted / overprinted copies (with the overprinting test):
    python3 scan_genome.py --seeds gene_seeds.faa --genome genome.fna --find-interrupted
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

import find_interrupted as FI   # read-through translation, coord mapping, overprinting

ROW_COLS = ["contig", "strand", "frame", "nt_start", "nt_end", "domain_aa_len",
            "internal_stops", "status", "domain_bit_score", "i_evalue",
            "orf_nt_start", "orf_nt_end", "orf_aa_len", "has_start_M", "ends_at_stop",
            "overprinting_support", "antisense_open_stops", "stop_nt_positions",
            "domain_aa", "orf_aa", "orf_nt"]


def fetch_genome(accessions: str, email: str | None, out_dir: Path, log) -> Path:
    """Pull one or more NCBI **nucleotide** accessions (comma-separated) into a single
    genome FASTA via Entrez efetch, so you can scan a genome by accession instead of a
    local file. NCBI requires an email. Retries transient failures. Returns the path."""
    from Bio import Entrez
    email = email or os.environ.get("NCBI_EMAIL")
    if not email:
        sys.exit("--accession needs an email for NCBI Entrez: pass --email you@inst.edu "
                 "(or set $NCBI_EMAIL).")
    Entrez.email = email
    socket.setdefaulttimeout(60)
    ids = [a.strip() for a in accessions.replace(",", " ").split() if a.strip()]
    if not ids:
        sys.exit("no accession given to --accession")
    out_fa = out_dir / ((ids[0].replace("/", "_") + ".fna") if len(ids) == 1
                        else "fetched_genome.fna")
    log(f"Fetching {len(ids)} accession(s) from NCBI nucleotide: {', '.join(ids)} …")
    data = ""
    for attempt in (1, 2, 3):
        try:
            with Entrez.efetch(db="nucleotide", id=",".join(ids),
                               rettype="fasta", retmode="text") as h:
                data = h.read()
            if data.strip().startswith(">"):
                break
        except Exception as e:
            log(f"  (NCBI fetch attempt {attempt} failed: {e})")
        time.sleep(3 * attempt)
    if not data.strip().startswith(">"):
        sys.exit(f"NCBI returned no FASTA for: {accessions} "
                 "(check the accession is a nucleotide record; assembly GCF_/GCA_ ids "
                 "aren't fetched directly — use their nucleotide/contig accessions).")
    out_fa.write_text(data)
    log(f"  fetched {data.count('>')} sequence(s) -> {out_fa.name}")
    return out_fa


def build_hmm_from_seeds(seeds: Path, table: int, out_dir: Path, cpu: int, log) -> Path:
    """Build a profile HMM from seed sequences (translating a nucleotide seed first).
    MAFFT-align (accuracy-first L-INS-i for ≤500 seqs, else --auto) then hmmbuild;
    a single seed is built directly. Same toolchain as the discovery pipeline."""
    import hmm_finder as H   # sibling: reuse the nucleotide-seed detection + translation
    seeds = Path(seeds)
    if H._looks_like_nucleotide(seeds):
        log(f"Seed looks like nucleotide CDS; translating with genetic code {table}…")
        seeds = H.translate_seed(seeds, table, out_dir, log)
    n = sum(1 for _ in SeqIO.parse(str(seeds), "fasta"))
    if n == 0:
        sys.exit("No sequences found in the seed FASTA.")
    aln = out_dir / "seeds.aln.faa"
    if n == 1:
        shutil.copy(str(seeds), str(aln))           # hmmbuild accepts a single sequence
    else:
        mode = (["--localpair", "--maxiterate", "1000"] if n <= 500 else ["--auto"])
        with open(aln, "w") as f:
            subprocess.run(["mafft", *mode, "--thread", str(cpu), str(seeds)],
                           check=True, stdout=f, stderr=subprocess.DEVNULL)
    hmm = out_dir / "profile.hmm"
    subprocess.run(["hmmbuild", "--amino", str(hmm), str(aln)],
                   check=True, capture_output=True, text=True)
    log(f"Built HMM from {n} seed(s) -> {hmm.name}")
    return hmm


def scan(genome: Path, hmm: Path, out_dir: Path, min_bit: float, find_interrupted: bool,
         cpu: int, log) -> dict:
    """Read-through six-frame scan of one genome with the HMM. Returns a summary dict
    and writes scan_hits.tsv + scan_hits_{aa.faa,nt.fna} + scan_report.txt."""
    out_dir.mkdir(parents=True, exist_ok=True)
    workdir = Path(tempfile.mkdtemp(prefix="scan_", dir=str(out_dir)))
    frames, markers, contig_nt = [], {}, {}
    n_contigs = 0
    try:
        with FI.open_maybe_gz(genome) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                seq = str(rec.seq).upper().replace("U", "T")
                contig_nt[rec.id] = seq
                n_contigs += 1
                for strand, frame, search, marker in FI._frames(seq):
                    for off, wlen in FI._windows(len(search)):
                        sw = search[off:off + wlen]
                        if len(sw) < FI.MIN_AA:
                            continue
                        name = f"{rec.id}__{strand}{frame}__{off}"
                        frames.append((name, sw))
                        markers[name] = marker[off:off + wlen]
        if not frames:
            log("  (genome too short to scan in any frame)")
            return _finish(out_dir, [], min_bit, find_interrupted, n_contigs, log)
        sfa = workdir / "frames.faa"
        with open(sfa, "w") as f:
            for name, search in frames:
                f.write(f">{name}\n{search}\n")
        dt = workdir / "scan.domtbl"
        log(f"hmmsearch: {len(frames)} read-through frame window(s) over {n_contigs} contig(s)…")
        subprocess.run(["hmmsearch", "--noali", "--cpu", str(cpu), "--domtblout", str(dt),
                        str(hmm), str(sfa)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rows = _parse(dt, markers, contig_nt, min_bit, find_interrupted)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return _finish(out_dir, rows, min_bit, find_interrupted, n_contigs, log)


def _parse(dt: Path, markers: dict, contig_nt: dict, min_bit: float,
           find_interrupted: bool) -> list:
    rows = []
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
        marker = markers.get(name, "")
        n_stops, positions = FI.count_envelope_stops(marker, env_from, env_to)
        if n_stops > 0 and not find_interrupted:
            continue   # interrupted copy, but the user asked only for clean detection
        contig, sf_tag, off_s = name.rsplit("__", 2)
        strand, frame, offset = sf_tag[0], int(sf_tag[1]), int(off_s)
        nt = contig_nt.get(contig, "")
        nt_start, nt_end, dom_nt = FI.aa_to_nt(nt, strand, frame, offset + env_from, offset + env_to)
        orf_from, orf_to, orf_aa = FI.extend_orf(marker, env_from, env_to)
        orf_nt_start, orf_nt_end, orf_nt = FI.aa_to_nt(nt, strand, frame, offset + orf_from, offset + orf_to)
        terminal_stop = orf_aa.endswith("*")
        # ORF validation: a clean gene starts with M (after the upstream stop) and ends at a stop.
        has_M = orf_aa[:1] == "M"
        status = "clean" if n_stops == 0 else "interrupted"
        overpr, anti_open = "", ""
        stop_pos = ""
        if n_stops > 0:
            stop_fwds = [FI.stop_nt(nt, strand, frame, offset + x) for x in positions]
            stop_pos = ";".join(str(s) for s in stop_fwds)
            opa = FI.analyze_overprinting(nt, strand, nt_start, nt_end, stop_fwds)
            overpr, anti_open = opa["support"], opa["open_stops"]
        rows.append({
            "contig": contig, "strand": strand, "frame": frame,
            "nt_start": nt_start, "nt_end": nt_end,
            "domain_aa_len": env_to - env_from + 1, "internal_stops": n_stops,
            "status": status, "domain_bit_score": round(dom_bits, 1),
            "i_evalue": f"{i_eval:.2g}",
            "orf_nt_start": orf_nt_start, "orf_nt_end": orf_nt_end,
            "orf_aa_len": orf_to - orf_from + 1,
            "has_start_M": int(has_M), "ends_at_stop": int(terminal_stop),
            "overprinting_support": overpr, "antisense_open_stops": anti_open,
            "stop_nt_positions": stop_pos,
            "domain_aa": marker[env_from - 1:env_to], "orf_aa": orf_aa, "orf_nt": orf_nt,
        })
    # one read-through frame can report a domain twice across overlapping windows; keep best
    best = {}
    for r in rows:
        key = (r["contig"], r["strand"], r["frame"], r["nt_start"], r["nt_end"])
        if key not in best or r["domain_bit_score"] > best[key]["domain_bit_score"]:
            best[key] = r
    out = sorted(best.values(), key=lambda r: -r["domain_bit_score"])
    return out


def _finish(out_dir: Path, rows: list, min_bit: float, find_interrupted: bool,
            n_contigs: int, log) -> dict:
    tsv = out_dir / "scan_hits.tsv"
    with open(tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ROW_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(out_dir / "scan_hits_aa.faa", "w") as fa, \
         open(out_dir / "scan_hits_nt.fna", "w") as fn:
        for i, r in enumerate(rows, 1):
            hdr = (f"{r['contig']}|{r['strand']}{r['frame']}|{r['nt_start']}-{r['nt_end']}|"
                   f"{r['status']}|bit={r['domain_bit_score']}")
            if r["orf_aa"]:
                fa.write(f">{hdr}\n{r['orf_aa']}\n")
            if r["orf_nt"]:
                fn.write(f">{hdr}\n{r['orf_nt']}\n")
    clean = [r for r in rows if r["status"] == "clean"]
    interrupted = [r for r in rows if r["status"] == "interrupted"]
    if clean:
        verdict = f"GENE PRESENT — {len(clean)} clean hit(s)"
    elif interrupted:
        verdict = (f"PRESENT but INTERRUPTED — {len(interrupted)} stop-interrupted/"
                   f"overprinted copy(ies), no clean ORF")
    else:
        verdict = "GENE NOT DETECTED"
    lines = [f"Single-genome scan — {verdict}",
             f"  contigs scanned: {n_contigs};  reporting threshold: {min_bit:g} bits"
             + ("  (read-through: clean + interrupted)" if find_interrupted else "  (clean ORFs only)"),
             ""]
    for r in (rows[:10] or []):
        line = (f"  • {r['contig']} {r['strand']}{r['frame']} {r['nt_start']}-{r['nt_end']} "
                f"| {r['status']} | bit {r['domain_bit_score']} E {r['i_evalue']} "
                f"| ORF {r['orf_aa_len']}aa (M={r['has_start_M']}, stop={r['ends_at_stop']})")
        if r["status"] == "interrupted":
            line += f" | overprinting={r['overprinting_support']} (antisense_open_stops={r['antisense_open_stops']})"
        lines.append(line)
    if len(rows) > 10:
        lines.append(f"  … and {len(rows) - 10} more in scan_hits.tsv")
    report = "\n".join(lines) + "\n"
    (out_dir / "scan_report.txt").write_text(report)
    log("\n" + report)
    log(f"  -> {tsv.name}, scan_hits_aa.faa, scan_hits_nt.fna, scan_report.txt")
    return {"verdict": verdict, "n_clean": len(clean), "n_interrupted": len(interrupted),
            "n_contigs": n_contigs, "hits": len(rows), "tsv": str(tsv)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--seeds", type=Path, help="seed FASTA (protein or nucleotide CDS) to BUILD the HMM from")
    src.add_argument("--hmm", type=Path, help="an existing profile HMM to scan with")
    gsrc = ap.add_mutually_exclusive_group(required=True)
    gsrc.add_argument("--genome", type=Path,
                      help="a local single nucleotide genome to scan (.fna/.fa, optionally .gz)")
    gsrc.add_argument("--accession",
                      help="NCBI nucleotide accession(s) to fetch & scan, comma-separated "
                           "(e.g. KX098390 or NC_031062). Needs --email.")
    ap.add_argument("--email", default=None,
                    help="NCBI email (required with --accession; or set $NCBI_EMAIL)")
    ap.add_argument("--out", type=Path, default=Path("genome_scan"), help="output directory")
    ap.add_argument("--min-bit", type=float, default=25.0,
                    help="minimum domain bit score to report (default 25)")
    ap.add_argument("--find-interrupted", action="store_true",
                    help="also report stop-interrupted / overprinted copies (with the overprinting test)")
    ap.add_argument("--trans-table", type=int, default=11,
                    help="genetic code for translating a nucleotide SEED (default 11)")
    ap.add_argument("--cpu", type=int, default=4)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    log = print
    # Resolve the genome: a local FASTA, or fetch it from NCBI by accession.
    if args.accession:
        genome = fetch_genome(args.accession, args.email, args.out, log)
    else:
        if not args.genome.exists():
            sys.exit(f"genome not found: {args.genome}")
        genome = args.genome
    if args.hmm:
        if not args.hmm.exists():
            sys.exit(f"HMM not found: {args.hmm}")
        hmm = args.hmm
    else:
        if not args.seeds.exists():
            sys.exit(f"seeds not found: {args.seeds}")
        hmm = build_hmm_from_seeds(args.seeds, args.trans_table, args.out, args.cpu, log)
    s = scan(genome, hmm, args.out, args.min_bit, args.find_interrupted, args.cpu, log)
    # exit code 0 if the gene was detected (clean or interrupted), 1 if absent — handy in scripts
    sys.exit(0 if (s["n_clean"] or s["n_interrupted"]) else 1)


if __name__ == "__main__":
    main()
