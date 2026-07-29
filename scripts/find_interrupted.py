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
import functools
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

STOP_SEARCH = "X"   # HMMER-friendly placeholder for a stop in the SEARCH sequence
MIN_AA = 30         # ignore frames shorter than this (matches the six-frame floor)
MET_MARGIN = 60     # a start codon may sit at most this many aa upstream of the domain envelope


def open_maybe_gz(p: Path):
    p = Path(p)
    return gzip.open(p, "rt") if p.suffix == ".gz" else open(p)


def read_through_aa(nt: str, strand: str, frame: int, table: int = 11) -> tuple[str, str]:
    """Translate one reading frame with read-through, under NCBI genetic code `table`
    (11 = bacterial/archaeal/phage default; e.g. 4 = Mycoplasma TGA->W, 25 = TGA->G).
    Returns (search_aa, marker_aa): search_aa has 'X' at stops (so the HMM matches
    across them); marker_aa keeps '*' at the same positions so stops can be located."""
    s = nt.upper().replace("U", "T")
    if strand == "-":
        s = str(Seq(s).reverse_complement())
    s = s[frame:]
    s = s[: len(s) // 3 * 3]
    marker = str(Seq(s).translate(table=table)) if s else ""   # '*' at stops
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


def extend_orf(marker_aa: str, env_from: int, env_to: int) -> tuple[int, int, str]:
    """Extend the domain envelope to the flanking stops in the read-through frame —
    the full ORF the gene occupies, read THROUGH any premature internal stops — so
    the sequence AFTER the premature stop (the rest of the gene, to its natural
    stop) is captured. Returns (orf_from, orf_to, full_aa): 1-based inclusive bounds
    and the full amino-acid sequence with '*' kept at every stop (premature ones in
    the middle; the terminal one is the natural gene end)."""
    n = len(marker_aa)
    left = marker_aa.rfind("*", 0, max(0, env_from - 1))     # upstream stop (0-based) or -1
    after_stop = left + 1 if left >= 0 else 0                # 0-based, first residue past it
    # Start the gene at its ATG/Met START CODON, not at the upstream stop, otherwise we prepend
    # the residues between the upstream in-frame stop and the real start (a 145-aa stop-to-stop
    # ORF for a 138-aa gene). Use the Met CLOSEST to and within MET_MARGIN aa upstream of the
    # domain envelope start (and after the upstream stop) — NOT the earliest Met, which on a long
    # stop-free (antisense/overprint) frame would be hundreds of residues upstream and is not the
    # gene start. Fall back to the domain start (no upstream extension) when no Met is near.
    lo = max(after_stop, env_from - 1 - MET_MARGIN)
    met = marker_aa.rfind("M", max(0, lo), max(0, env_from))   # closest Met at/before env_from
    orf_from = (met + 1) if met >= 0 else env_from              # 1-based
    right = marker_aa.find("*", env_to)                        # first stop at/after the domain
    orf_to = (right + 1) if right >= 0 else n                  # include the natural terminal stop
    return orf_from, orf_to, marker_aa[orf_from - 1:orf_to]


def aa_to_nt(contig_seq: str, strand: str, frame: int, aa_from: int, aa_to: int):
    """Map a frame-relative aa span [aa_from, aa_to] (1-based inclusive) back to the
    genome: returns (fwd_start, fwd_end, coding_nt) where fwd_start/fwd_end are
    FORWARD-strand contig coordinates (1-based inclusive) and coding_nt is the
    coding DNA 5'->3' (reverse-complemented for '-' hits). Lets you see exactly
    where in the genome the gene sits and pull its DNA."""
    L = len(contig_seq)
    if strand == "-":
        rc_a = frame + (aa_from - 1) * 3
        rc_b = frame + aa_to * 3
        fwd_start, fwd_end = L - rc_b + 1, L - rc_a
        coding = str(Seq(contig_seq[fwd_start - 1:fwd_end]).reverse_complement())
    else:
        fwd_start = frame + (aa_from - 1) * 3 + 1
        fwd_end = frame + aa_to * 3
        coding = contig_seq[fwd_start - 1:fwd_end]
    return fwd_start, fwd_end, coding


def stop_nt(contig_seq: str, strand: str, frame: int, aa_pos: int) -> int:
    """Forward-strand 1-based coordinate of the first base of the stop codon at the
    given frame-relative aa position (lowest of the codon's three forward coords)."""
    if strand == "-":
        return len(contig_seq) - (frame + (aa_pos - 1) * 3) - 2
    return frame + (aa_pos - 1) * 3 + 1


# ---------------------------------------------------------------------------
# Overprinting / silent-stop analysis — the PROOF (not just the location) of an
# overprinted homolog. A premature stop is the overprinting signature when it (1)
# sits inside an OPEN overlapping antisense reading frame (a real candidate gene,
# e.g. gp75's RNA polymerase) and (2) is SYNONYMOUS in that frame — i.e. the
# nonsense mutation can be reverted to sense in the small gene without changing the
# antisense protein, so selection on the antisense gene tolerates (or favours) it.
# ---------------------------------------------------------------------------
_COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _comp(b: str) -> str:
    return b.translate(_COMP)


@functools.lru_cache(maxsize=8192)
def _aa(codon: str, table: int = 11) -> str:
    """Translate one codon (genetic code 11 = bacterial/phage by default). Ambiguous
    or short codons -> 'X'. Cached: codons repeat heavily across the scan."""
    if len(codon) != 3 or any(c not in "ACGT" for c in codon):
        return "X"
    return str(Seq(codon).translate(table=table))


def _codon_covering(contig: str, strand: str, frame: int, fwd_i0: int):
    """The codon (read 5'->3' in the gene) and the 0-based index within it that cover
    forward 0-based position fwd_i0 in reading frame (strand, frame). Returns
    (codon_str, idx) or None if the position falls outside a complete codon. Works on
    a local 3-base slice — no whole-contig reverse-complement."""
    L = len(contig)
    if strand == "+":
        k = (fwd_i0 - frame) // 3
        s0 = frame + 3 * k
        if s0 < 0 or s0 + 3 > L:
            return None
        return contig[s0:s0 + 3], fwd_i0 - s0
    r = L - 1 - fwd_i0                       # position in the reverse-complement strand
    k = (r - frame) // 3
    s0 = frame + 3 * k
    if s0 < 0 or s0 + 3 > L:
        return None
    codon = _revcomp(contig[L - s0 - 3:L - s0])   # R[s0:s0+3] expressed from forward bases
    return codon, r - s0


def _frame_stop_count(contig: str, strand: str, frame: int, lo1: int, hi1: int) -> int:
    """Number of stop codons in reading frame (strand, frame) whose codon lies fully
    within forward window [lo1, hi1] (1-based inclusive). 0 == the frame is OPEN across
    the window (a candidate overlapping ORF)."""
    L = len(contig)
    s = contig if strand == "+" else _revcomp(contig)
    cnt, cs = 0, frame
    while cs + 3 <= len(s):
        if strand == "+":
            flo, fhi = cs + 1, cs + 3
        else:
            flo, fhi = L - cs - 2, L - cs
        if flo >= lo1 and fhi <= hi1 and _aa(s[cs:cs + 3]) == "*":
            cnt += 1
        cs += 3
    return cnt


def _codon_at_pos(contig: str, strand: str, frame: int, fwd_i0: int):
    """(codon, fwd_lo0, fwd_hi0) for the frame codon covering forward 0-based `fwd_i0`, or
    None off the contig. A thin wrapper over `_codon_covering` that additionally reports
    WHERE the codon sits, so a caller can step codon-by-codon without re-deriving the frame
    arithmetic (one source of truth). Note the codon is read 5'->3' in the gene, so on '-'
    codon[0] is the HIGHEST forward base — hence the two branches below."""
    cov = _codon_covering(contig, strand, frame, fwd_i0)
    if not cov:
        return None
    codon, idx = cov
    lo0 = (fwd_i0 - idx) if strand == "+" else (fwd_i0 + idx - 2)
    if lo0 < 0 or lo0 + 3 > len(contig):
        return None
    return codon, lo0, lo0 + 2


def antisense_open_extent(contig: str, anti_strand: str, frame: int,
                          dom_lo1: int, dom_hi1: int, table: int = 11):
    """Full extent of the OPEN reading frame that the domain is nested inside, measured in
    the overlapping (antisense) frame (`anti_strand`, `frame`).

    WHY: `_frame_stop_count` only asks whether the DOMAIN WINDOW is stop-free, so the true
    length of the overlapping ORF is never computed — yet that length is the discriminating,
    length-dependent half of the overprinting argument (an open frame across L codons is
    improbable by chance, and the improbability grows with L; docs/METHODOLOGY.md quotes
    ~0.7% for the 137-aa gp75 length). Reporting a REAL L instead of an illustrative one is
    the difference between an anecdote and a measurement.

    Walks outward from the domain window, codon by codon in that frame, to the first stop in
    each direction — reusing `_codon_covering` via `_codon_at_pos`, so no whole-contig
    reverse-complement is ever built (these contigs are whole phage genomes).

    Returns (orf_lo1, orf_hi1, orf_aa_len): FORWARD-strand 1-based inclusive nucleotide
    bounds of the stop-free stretch — the flanking stop codons themselves are EXCLUDED — and
    its length in codons. Returns (0, 0, 0) when the frame or coordinates cannot be resolved
    (house style: never raise; the caller leaves the column empty).

    SELF-GATING (was a real bug): the returned span is only ever a stop-free ORF, never an
    "envelope". The starting envelope — the two antisense codons COVERING the domain edges —
    is stop-tested BEFORE the outward walk, and so is every codon between them; if any is a
    stop this returns (0, 0, 0). Without that test a locus whose antisense frame is closed
    (`antisense_open_stops` > 0) still got a span back, and a caller could publish a
    stop-riddled interval as "the computed open antisense ORF". A non-zero length from this
    function now MEANS the frame is open across the whole domain window.

    STILL NOT EVIDENCE OF EXPRESSION: an open frame is a necessary sequence signature of
    overprinting, not proof the antisense ORF is transcribed or translated.
    """
    try:
        L = len(contig or "")
        if L < 3 or anti_strand not in ("+", "-") or int(frame) not in (0, 1, 2):
            return 0, 0, 0
        frame = int(frame)
        d_lo, d_hi = int(dom_lo1), int(dom_hi1)
        if d_hi < d_lo:
            d_lo, d_hi = d_hi, d_lo
        if d_hi < 1 or d_lo > L:
            return 0, 0, 0        # window entirely off the contig: clamping it to the edge
                                  # would report an ORF for a locus that isn't on this contig
        lo0 = max(0, min(L - 1, d_lo - 1))     # a 1-2 bp overhang at a contig end is normal
        hi0 = max(0, min(L - 1, d_hi - 1))     # for a windowed frame, so clamp rather than bail
        # The domain's edges are codon boundaries of the SMALL gene, not of this frame, so
        # take the antisense codons that COVER those edges as the starting envelope.
        a = _codon_at_pos(contig, anti_strand, frame, lo0)
        b = _codon_at_pos(contig, anti_strand, frame, hi0)
        if not a or not b:
            return 0, 0, 0
        # The envelope codons are part of the span this function returns, so they must be
        # stop-tested like any other codon of the ORF. (They are also exactly the codons
        # `_frame_stop_count` cannot see: it counts only codons lying FULLY inside the window,
        # so an edge-straddling stop is invisible to `antisense_open_stops`.)
        if _aa(a[0], table) == "*" or _aa(b[0], table) == "*":
            return 0, 0, 0
        left, right = min(a[1], b[1]), max(a[2], b[2])
        # ...and so is everything between them: the domain window is NOT taken on trust. A
        # closed antisense frame has no ORF containing the domain, and reporting its envelope
        # as one is how a stop-riddled interval ends up in a CSV labelled "open ORF".
        cur = left
        while cur + 2 < right:          # `left` is codon a; stop once cur IS the last codon
            c = _codon_at_pos(contig, anti_strand, frame, cur + 3)
            if not c:
                return 0, 0, 0
            if _aa(c[0], table) == "*":
                return 0, 0, 0          # frame closed across the domain -> no ORF to report
            cur = c[1]
        # Codons of one frame tile the contig in steps of 3 on BOTH strands, so stepping the
        # forward coordinate by 3 lands exactly on the neighbouring codon either way.
        cur = left
        while True:
            c = _codon_at_pos(contig, anti_strand, frame, cur - 3)
            if not c or _aa(c[0], table) == "*":
                break                                  # contig edge, or the flanking stop
            cur = c[1]
        orf_lo0 = cur
        cur = right
        while True:
            c = _codon_at_pos(contig, anti_strand, frame, cur + 3)
            if not c or _aa(c[0], table) == "*":
                break
            cur = c[2]
        orf_hi0 = cur
        return orf_lo0 + 1, orf_hi0 + 1, (orf_hi0 - orf_lo0 + 1) // 3
    except Exception:
        return 0, 0, 0


def _stop_silent_in_frame(contig: str, small_strand: str, stop_fwd1: int,
                          anti: str, frame: int, table: int = 11) -> bool:
    """Is the premature stop (small-gene codon at forward 1-based stop_fwd1..+2)
    SYNONYMOUS in antisense reading frame `frame`? True iff some single-base
    substitution that REMOVES the stop in the small gene leaves the antisense frame's
    amino acid unchanged (and that antisense codon actually encodes a residue, not a
    stop)."""
    L = len(contig)
    i0 = stop_fwd1 - 1
    if i0 < 0 or i0 + 3 > L:
        return False
    tri = contig[i0:i0 + 3]
    if any(c not in "ACGT" for c in tri):
        return False
    small = tri if small_strand == "+" else _revcomp(tri)
    if _aa(small, table) != "*":
        return False                         # not a stop in the small frame -> nothing to test
    for j in range(3):
        for b in "ACGT":
            if b == tri[j]:
                continue
            mtri = tri[:j] + b + tri[j + 1:]
            mc = mtri if small_strand == "+" else _revcomp(mtri)
            if _aa(mc, table) == "*":
                continue                     # still a stop -> not a stop-removing change
            cov = _codon_covering(contig, anti, frame, i0 + j)
            if not cov:
                continue
            codon, idx = cov
            oaa = _aa(codon, table)
            if oaa in ("", "X", "*"):
                continue                     # antisense frame not coding here
            nb = b if anti == "+" else _comp(b)
            if _aa(codon[:idx] + nb + codon[idx + 1:], table) == oaa:
                return True
    return False


def analyze_overprinting(contig: str, small_strand: str, dom_lo1: int, dom_hi1: int,
                         stop_fwds: list, table: int = 11) -> dict:
    """Per interrupted homolog: pick the antisense frame with the FEWEST stops across
    the domain (the candidate overprinted ORF), then test whether each premature stop
    is synonymous in that frame. Returns {open_frame, open_stops, per_stop_silent,
    support}. support: 'strong' = the antisense frame is fully open (0 stops) AND every
    premature stop is silent in it; 'partial' = some stops silent; 'none' = no evidence."""
    anti = "-" if small_strand == "+" else "+"
    counts = [(f, _frame_stop_count(contig, anti, f, dom_lo1, dom_hi1)) for f in (0, 1, 2)]
    open_frame, open_stops = min(counts, key=lambda c: c[1])
    per_stop = [_stop_silent_in_frame(contig, small_strand, sp, anti, open_frame, table)
                for sp in stop_fwds]
    n_sil = sum(per_stop)
    if per_stop and open_stops == 0 and n_sil == len(per_stop):
        support = "strong"
    elif n_sil:
        support = "partial"
    else:
        support = "none"
    return {"open_frame": open_frame, "open_stops": open_stops,
            "per_stop_silent": per_stop, "support": support}


def _frames(seq: str, table: int = 11):
    """Yield (strand, frame, search_aa, marker_aa) for all six frames."""
    for strand in ("+", "-"):
        for frame in (0, 1, 2):
            search, marker = read_through_aa(seq, strand, frame, table)
            if len(search) >= MIN_AA:
                yield strand, frame, search, marker


BATCH_CONTIGS = 500   # search this many contigs per hmmsearch call (bounds temp + memory)
WIN_AA = 5000         # window long read-through frames: HMMER aborts (SIGABRT) on whole-
WIN_STEP = 4400       # genome-length sequences (10k-100k+ aa). Overlap = WIN_AA - WIN_STEP
                      # (>> any domain), so every domain lies fully inside >=1 window.
ROW_COLS = ["contig", "strand", "frame", "domain_nt_start", "domain_nt_end",
            "domain_aa_len", "internal_stops", "stop_nt_positions", "stop_aa_positions",
            "overprinting_support", "antisense_open_frame", "antisense_open_stops",
            "stop_silent_antisense",
            "domain_bit_score", "i_evalue", "orf_aa_len",
            "aa_before_first_stop", "aa_after_last_stop",
            "orf_nt_start", "orf_nt_end", "natural_stop_nt",
            "domain_nt", "domain_aa_with_stops", "full_orf_aa", "full_orf_nt"]


def _fasta_header(r: dict) -> str:
    """Descriptive, unique FASTA header for one interrupted-homolog row."""
    return (f"{r['contig']}__{r['strand']}{r['frame']}__"
            f"{r.get('domain_nt_start', '?')}-{r.get('domain_nt_end', '?')} "
            f"stops={r.get('internal_stops', '?')}@{r.get('stop_aa_positions', '')} "
            f"bit={r.get('domain_bit_score', '?')} ievalue={r.get('i_evalue', '')}")


def write_aa_fastas(rows: list, out_tsv: Path) -> list:
    """Write the PROTEIN sequences of the interrupted homologs beside the TSV, with
    the internal stop(s) shown as '*':
      <stem>_domain_aa.faa   — the HMM-matched domain ('*' at every internal stop)
      <stem>_full_orf_aa.faa — the full read-through ORF (premature stops kept '*',
                               terminal stop is the natural gene end)
    Returns the two paths (always written, even if empty, so links never dangle)."""
    stem = Path(out_tsv).with_suffix("")          # strip .tsv
    dom = Path(f"{stem}_domain_aa.faa")
    orf = Path(f"{stem}_full_orf_aa.faa")
    with open(dom, "w") as dh, open(orf, "w") as oh:
        for r in rows:
            hdr = _fasta_header(r)
            if r.get("domain_aa_with_stops"):
                dh.write(f">{hdr}\n{r['domain_aa_with_stops']}\n")
            if r.get("full_orf_aa"):
                oh.write(f">{hdr}\n{r['full_orf_aa']}\n")
    return [dom, orf]


def write_orf_nt_fasta(rows: list, out_tsv: Path) -> Path:
    """Write the NUCLEOTIDE of each interrupted homolog's full read-through ORF
    (coding 5'->3', ending in the actual stop codon triplet) -> <stem>_full_orf_nt.fna.
    Translates back to full_orf_aa. Always written so report/PACKAGE links never dangle."""
    nt = Path(f"{Path(out_tsv).with_suffix('')}_full_orf_nt.fna")
    with open(nt, "w") as nh:
        for r in rows:
            if r.get("full_orf_nt"):
                nh.write(f">{_fasta_header(r)}\n{r['full_orf_nt']}\n")
    return nt


def add_organism_column(tsv_path: Path, meta_csv: Path) -> bool:
    """Insert an 'organism' column (right after 'contig') into an interrupted-homologs
    TSV, resolved offline from genome_metadata.csv (genome_id -> organism). Falls back to
    the contig id when an organism is unknown. Idempotent (no-op if 'organism' already
    present). Returns True if the file was rewritten."""
    import csv
    tsv_path, meta_csv = Path(tsv_path), Path(meta_csv)
    if not tsv_path.exists():
        return False
    with open(tsv_path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or "organism" in rows[0]:
        return False
    import re
    def _base(a):   # strip a trailing version (NC_049456.2 -> NC_049456) so a versioned
        return re.sub(r"\.\d+$", "", (a or "").strip())   # contig matches a base-accession key
    org = {}
    if meta_csv.exists():
        with open(meta_csv, newline="") as fh:
            for m in csv.DictReader(fh):
                o = m.get("organism", "")
                if m.get("genome_id"):
                    org[_base(m["genome_id"])] = o
                for a in (m.get("accessions", "") or "").split(";"):   # all known aliases
                    if a.strip():
                        org[_base(a)] = o
    for r in rows:
        c = r.get("contig", "")
        r["organism"] = org.get(_base(c), "") or c
    cols = []
    for c in list(rows[0].keys()):
        if c == "organism":
            continue
        cols.append(c)
        if c == "contig":
            cols.append("organism")
    if "organism" not in cols:
        cols.insert(0, "organism")
    with open(tsv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return True


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


def _search_batch(frames: list, markers: dict, contig_nt: dict, hmm: Path, workdir: Path,
                  min_bit: float, cpu: int, table: int = 11) -> tuple[int, list]:
    """hmmsearch one batch of (name, search_aa); return (n_scored, interrupted_rows).
    markers maps name -> marker_aa; contig_nt maps contig -> forward DNA (both for
    the SAME batch, kept small)."""
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
        strand, frame, offset = sf_tag[0], int(sf_tag[1]), int(off_s)
        # full read-through ORF (through the premature stop to the natural stop).
        orf_from, orf_to, full_aa = extend_orf(marker, env_from, env_to)
        terminal_stop = 1 if full_aa.endswith("*") else 0
        aa_before = positions[0] - orf_from              # intact N-terminal residues
        aa_after = max(0, orf_to - positions[-1] - terminal_stop)  # C-terminal continuation
        # Map the DOMAIN back to genome DNA coordinates (frame-relative -> contig).
        nt = contig_nt.get(contig, "")
        dfrom, dto = offset + env_from, offset + env_to       # frame-relative aa
        nt_start, nt_end, dom_nt = (aa_to_nt(nt, strand, frame, dfrom, dto)
                                    if nt else (0, 0, ""))
        stop_fwds = [stop_nt(nt, strand, frame, offset + x) for x in positions] if nt else []
        stop_nt_pos = ";".join(str(s) for s in stop_fwds)
        # Overprinting (silent-stop) test: is the premature stop synonymous in an OPEN
        # overlapping antisense frame? (the proof of overprinting, not just location).
        if nt and stop_fwds:
            opa = analyze_overprinting(nt, strand, nt_start, nt_end, stop_fwds, table)
            anti_open_frame, anti_open_stops = opa["open_frame"], opa["open_stops"]
            stop_silent = ";".join("1" if s else "0" for s in opa["per_stop_silent"])
            overpr = opa["support"]
        else:
            anti_open_frame, anti_open_stops, stop_silent, overpr = -1, -1, "", "none"
        # Map the FULL read-through ORF back to genome DNA — its coding sequence
        # 5'->3' includes the actual stop codon triplet at the end; and locate the
        # natural/terminal stop codon's forward coordinate (0 if the ORF ran to the
        # contig/window edge without one).
        ofrom, oto = offset + orf_from, offset + orf_to
        orf_nt_start, orf_nt_end, orf_nt = (aa_to_nt(nt, strand, frame, ofrom, oto)
                                            if nt else (0, 0, ""))
        natural_stop_nt = (stop_nt(nt, strand, frame, oto) if (nt and terminal_stop) else 0)
        rows.append({
            "contig": contig, "strand": strand, "frame": frame,
            "domain_nt_start": nt_start, "domain_nt_end": nt_end,
            "domain_aa_len": env_to - env_from + 1,
            "internal_stops": n_stops,
            "stop_nt_positions": stop_nt_pos,
            "stop_aa_positions": ";".join(str(offset + x) for x in positions),
            "overprinting_support": overpr,
            "antisense_open_frame": anti_open_frame,
            "antisense_open_stops": anti_open_stops,
            "stop_silent_antisense": stop_silent,
            "domain_bit_score": round(dom_bits, 1),
            "i_evalue": f"{i_eval:.2g}",
            # protein length EXCLUDING the terminal stop codon (a '*' is not a residue), so this
            # column matches scan_genome.py and extract_validated_hits.py for the same gene.
            "orf_aa_len": (orf_to - orf_from + 1) - int(terminal_stop),
            "aa_before_first_stop": aa_before,
            "aa_after_last_stop": aa_after,
            "orf_nt_start": orf_nt_start,
            "orf_nt_end": orf_nt_end,
            "natural_stop_nt": natural_stop_nt,
            "domain_nt": dom_nt,
            "domain_aa_with_stops": marker[env_from - 1:env_to],
            "full_orf_aa": full_aa,
            "full_orf_nt": orf_nt,
        })
    sfa.unlink(missing_ok=True)
    dt.unlink(missing_ok=True)
    return n_scored, rows


def _run(genomes: Path, hmm: Path, out: Path, min_bit: float, cpu: int, log=print,
         emit_fasta: bool = True, table: int = 11) -> dict:
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
    frames, markers, contig_nt, nctg = [], {}, {}, 0
    try:
        with open_maybe_gz(genomes) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                seq = str(rec.seq)
                contig_nt[rec.id] = seq.upper().replace("U", "T")   # for genome coords + DNA
                for strand, frame, search, marker in _frames(seq, table):
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
                        ns, rows = _search_batch(frames, markers, contig_nt, hmm, workdir, min_bit, cpu, table)
                        n_scored += ns
                        all_rows += rows
                    except Exception as e:
                        log(f"  (find-interrupted: a batch failed, skipped: {e})")
                    n_batches += 1
                    frames, markers, contig_nt, nctg = [], {}, {}, 0
            if frames:
                try:
                    ns, rows = _search_batch(frames, markers, contig_nt, hmm, workdir, min_bit, cpu, table)
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
    summary = {"matches_scored": n_scored, "interrupted_candidates": len(all_rows),
               "out": str(out)}
    if emit_fasta:
        dom, orf = write_aa_fastas(all_rows, out)
        nt = write_orf_nt_fasta(all_rows, out)
        summary["domain_faa"], summary["orf_faa"] = str(dom), str(orf)
        summary["orf_fna"] = str(nt)
        log(f"  find-interrupted: sequences -> {dom.name}, {orf.name}, {nt.name}")
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
    ap.add_argument("--trans-table", type=int, default=11,
                    help="NCBI genetic code for the read-through translation (default 11 = "
                         "bacterial/archaeal/phage; e.g. 4 = Mycoplasma, 25 = candidate SR1)")
    args = ap.parse_args()
    if not args.genomes.exists():
        sys.exit(f"genomes FASTA not found: {args.genomes}")
    if not args.hmm.exists():
        sys.exit(f"HMM not found: {args.hmm}")
    s = _run(args.genomes, args.hmm, args.out, args.min_bit, args.cpu, table=args.trans_table)
    if "error" in s:
        sys.exit(s["error"])
    print(f"  matches scored: {s['matches_scored']}; interrupted candidates: "
          f"{s['interrupted_candidates']}; written to {s['out']}")


if __name__ == "__main__":
    main()
