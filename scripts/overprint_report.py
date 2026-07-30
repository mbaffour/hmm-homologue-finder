#!/usr/bin/env python3
"""overprint_report.py — NAME the host gene each interrupted homolog is overprinted on.

WHY THIS EXISTS
---------------
find_interrupted.py already detects the overprinting signature: a premature stop in the
small gene that is synonymous inside an OPEN overlapping antisense frame. But it stores
that frame only as a number (`antisense_open_frame` = 0/1/2) plus a stop count. Nothing in
the pipeline ever says WHICH GENE the homolog is nested inside — which, for this project,
is the headline biology: gp75 is unannotated precisely because it sits antisense inside a
virion RNA polymerase, and "antisense frame 1 is open" is a far weaker sentence than
"nested inside the 10 kb virion RNA polymerase of Erwinia phage vB_EamP_Rexella".

This stage closes that gap. For every interrupted locus it:
  * collapses the SAME locus reported under several database accessions (RefSeq / GenBank /
    versioned RefSeq of one genome) to ONE row keyed on canonical.canonical_organism + strand
    + domain coordinates, keeping the best-scoring record and listing the rest in
    `accession_aliases` — so an OP id names a genomic locus, not a database record;
  * resolves the contig to a nucleotide accession (engine.pipeline.synteny.nt_accession_in);
  * pulls the genome's REAL CDS annotation from NCBI (build_real_genbanks.fetch_ncbi);
  * picks the host gene — the greatest-overlap annotated CDS on the ANTISENSE strand;
  * measures the TRUE extent of the open antisense ORF
    (find_interrupted.antisense_open_extent) and asks whether that computed ORF and the
    annotated host CDS are the same gene (`antisense_orf_matches_host_gene`) — the
    strongest single line of evidence this pipeline can produce offline, because it ties a
    frame computed from raw sequence to a gene named by an independent annotator;
  * writes overprinted_loci.csv + overprinting_summary.csv and renders a family overview
    plus a per-locus antisense diagram.

SCIENTIFIC HEDGE (kept everywhere, including every figure caption)
-----------------------------------------------------------------
An open overlapping frame carrying a synonymous premature stop is a NECESSARY SEQUENCE
SIGNATURE of antisense overprinting. It is NOT proof that the antisense ORF is transcribed,
translated or selected — that needs orthogonal data (RNA-seq / ribo-seq, dN/dS in the
overlapping frame). Nothing written here should be quoted as evidence of expression.

OFFLINE
-------
Without --email nothing is fetched: every table is still written, host columns are empty,
`host_annotation_source` is "none", and the family overview is still rendered. The feature
degrades, it never disappears.

  python3 overprint_report.py --discovery-dir DIR [--email you@example.com] \
      [--cpu 8] [--out DIR]
"""
from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(HERE.parent / "engine")):      # sibling scripts/ + engine package
    if _p not in sys.path:
        sys.path.insert(0, _p)

STOP_COLOUR = "#c0392b"
SUPPORT_COLOURS = {"strong": "#1f6f4a", "partial": "#d98c00", "none": "#9aa3ad"}
SUPPORT_ORDER = {"strong": 0, "partial": 1, "none": 2}
MAX_LOCUS_FIGURES = 40      # per-locus diagrams are ~1 s each; a 500-locus family must not
                            # turn a 2-minute stage into a 10-minute one
FLANKS = 3                  # CDS drawn each side of the overprinted window, for context
P_STOP_PER_CODON = 3 / 64   # 3 of the 64 codons are stops under a uniform-composition null
N_ANTISENSE_FRAMES = 3      # analyze_overprinting picks the best of the 3 antisense frames

HEDGE = ("necessary sequence signature of antisense overprinting; NOT evidence that the "
         "antisense ORF is transcribed or translated")
# Pre-capitalised: str.capitalize() would lower-case "ORF"/"NOT" and quietly soften the caveat.
HEDGE_SENTENCE = "N" + HEDGE[1:] + "."

LOCI_COLS = [
    "locus_id", "contig", "organism", "accession", "accession_aliases", "strand", "frame",
    "domain_nt_start", "domain_nt_end", "domain_aa_len",
    "internal_stops", "stop_nt_positions",
    "overprinting_support", "antisense_open_frame", "antisense_open_stops",
    "stop_silent_antisense",
    "host_gene", "host_product", "host_locus_tag", "host_protein_id",
    "host_gene_start", "host_gene_end", "host_gene_strand", "host_gene_aa_len",
    "host_annotation_source", "host_category", "nested_fully",
    "overlap_bp", "overlap_pct_of_domain",
    "antisense_orf_is_open",
    "antisense_orf_nt_start", "antisense_orf_nt_end", "antisense_orf_aa_len",
    "antisense_orf_matches_host_gene", "figure",
]

SUMMARY_COLS = ["metric", "value", "note"]


# ---------------------------------------------------------------------------
# Guarded imports. Every one of these is optional at runtime: a missing engine
# module, a missing matplotlib or an absent DNA Features Viewer must degrade the
# report, never abort the pipeline stage that calls us.
# ---------------------------------------------------------------------------
def _accession_of(text: str) -> str:
    """Nucleotide accession embedded in a contig id / FASTA header, or ''. Delegates to
    engine.pipeline.synteny.nt_accession_in, which (unlike a \\b-anchored regex) resolves the
    underscore-joined seed headers this project's inputs use."""
    try:
        from pipeline.synteny import nt_accession_in
        return nt_accession_in(text) or ""
    except Exception:
        return ""


def _categorize(product: str, gene: str = "") -> str:
    """Broad functional category for a host product, using the SAME categorizer as the
    synteny figures and genome maps so a host gene is coloured/labelled identically there."""
    try:
        from synteny_figure import categorize
        return categorize(product or gene or "")
    except Exception:
        return ""


def _safe(name: str, limit: int = 60) -> str:
    """Filesystem-safe fragment for a figure name."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name or "")).strip("_")[:limit] or "locus"


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def read_interrupted(tsv_path) -> list:
    """Rows of an interrupted_homologs.tsv as dicts (empty list if absent/unreadable).
    The heavy sequence columns (domain_nt, full_orf_aa/nt) are dropped: this stage only
    needs coordinates, and carrying whole ORFs through would bloat memory pointlessly."""
    drop = {"domain_nt", "domain_aa_with_stops", "full_orf_aa", "full_orf_nt"}
    try:
        p = Path(tsv_path)
        if not p.exists():
            return []
        # explicit utf-8: organism strings carry accented phage/host names, and the default
        # locale codec on a Windows host would raise on them rather than read the file
        with open(p, newline="", encoding="utf-8", errors="replace") as fh:
            return [{k: v for k, v in r.items() if k not in drop}
                    for r in csv.DictReader(fh, delimiter="\t")]
    except Exception:
        return []


def _int(v, default=0) -> int:
    try:
        return int(float(str(v).strip()))
    except Exception:
        return default


def _float(v, default=0.0) -> float:
    try:
        return float(str(v).strip())
    except Exception:
        return default


def _primary_accession(row: dict) -> str:
    """The nucleotide accession a row's contig id resolves to (the contig id itself if it
    carries no recognisable accession)."""
    contig = str(row.get("contig", "") or "")
    return _accession_of(contig) or contig


def _canonical(row: dict) -> str:
    """`canonical.canonical_organism` for a row — the project-wide definition of "the same
    phage", used BOTH for the dedup key and for counting distinct organisms so those two
    numbers can never disagree. Works on a raw interrupted row (contig) or a built record
    (accession). Falls back to a lower-cased organism/accession if canonical.py is missing,
    which degrades the collapsing but never crashes the stage."""
    acc = str(row.get("accession", "") or "") or _primary_accession(row)
    name = str(row.get("organism", "") or "")
    try:
        from canonical import canonical_organism
        return canonical_organism(name, acc)
    except Exception:
        return (name.strip() or acc).lower()


# ---------------------------------------------------------------------------
# Cross-accession deduplication
# ---------------------------------------------------------------------------
def dedupe_rows(rows: list) -> list:
    """Collapse the SAME genomic locus reported under several database accessions.

    THE BUG THIS FIXES: the interrupted table carries one row per contig searched, and the
    same phage genome is catalogued under a RefSeq accession, a GenBank accession and a
    versioned RefSeq accession (NC_031062 / KX098389 / NC_031062.2 for Erwinia phage
    vB_EamP_Frozen). Numbering those three rows OP002 / OP004 / OP008 published 27 "loci"
    for 13 real ones, under a note claiming "one genomic locus each", in a table that
    simultaneously reported 13 distinct organisms — and contradicting this project's
    established one-locus-per-phage census.

    KEY: (canonical.canonical_organism(organism, accession), strand, domain_nt_start,
    domain_nt_end) — the same identity function every other table in this project counts
    organisms with, so "the same phage" means one thing pipeline-wide. Coordinates are part
    of the key on purpose: a genome carrying two genuinely different interrupted loci keeps
    both rows.

    REPRESENTATIVE: the best-scoring row of the group (highest domain_bit_score), ties broken
    by input order — which is already descending bit score, so the choice is deterministic
    and re-running the stage cannot reshuffle the OP ids.

    NOTHING IS LOST: every other accession of the group is preserved in `accession_aliases`
    (';'-joined, '' for a locus seen under a single accession), so `accession` +
    `accession_aliases` is the complete set of records this locus was found in.

    Idempotent — existing `accession_aliases` values are folded back in, so re-running this
    over its own output is a no-op rather than a silent loss of the alias list. Returns NEW
    dicts; the input rows are not mutated.
    """
    groups: dict = {}
    order: list = []
    for i, r in enumerate(rows or []):
        key = (_canonical(r), str(r.get("strand", "") or "").strip(),
               _int(r.get("domain_nt_start")), _int(r.get("domain_nt_end")))
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append((i, r))
    out = []
    for key in order:
        members = groups[key]
        idx, best = min(members, key=lambda ir: (-_float(ir[1].get("domain_bit_score")), ir[0]))
        accs: list = []
        for _i, r in members:
            for a in [_primary_accession(r)] + str(r.get("accession_aliases", "") or "").split(";"):
                a = a.strip()
                if a and a not in accs:
                    accs.append(a)
        rep = _primary_accession(best)
        rec = dict(best)
        rec["accession_aliases"] = ";".join(a for a in accs if a != rep)
        out.append((idx, rec))
    out.sort(key=lambda t: t[0])          # keep input order = descending domain bit score
    return [r for _i, r in out]


def _antisense_strand(small_strand: str) -> int:
    """+1/-1 of the frame OPPOSITE the small gene — the strand a host gene must lie on to be
    an antisense overprinting partner."""
    return -1 if str(small_strand).strip() == "+" else 1


# ---------------------------------------------------------------------------
# Host-gene assignment
# ---------------------------------------------------------------------------
def pick_host_gene(domain: dict, cds_list: list) -> dict:
    """The annotated CDS the interrupted domain is overprinted on: the GREATEST-OVERLAP CDS
    on the ANTISENSE strand whose interval intersects the domain. `{}` when there is none.

    `domain` = {"nt_start", "nt_end", "strand"} (1-based inclusive forward coordinates,
    strand '+'/'-' of the SMALL gene). `cds_list` = the tuples build_real_genbanks.fetch_ncbi
    returns: (start, end, strand, product, gene, translation), 1-based inclusive; longer
    tuples are tolerated so the two extra identifier fields are picked up automatically if
    that function ever grows them.

    The intersection test is scan_genome._select_neighbours' rule verbatim
    (`g_start <= dom_end and g_end >= dom_start`) so "overlapping partner" means the same
    thing in the neighbourhood tables and here. Greatest overlap (not nearest, not first)
    because a long host gene may be flanked by short CDS that clip the domain by a few bp.

    Returns {start, end, strand, product, gene, locus_tag, protein_id, aa_len, overlap_bp,
    nested_fully, category}.
    """
    try:
        d_lo, d_hi = _int(domain.get("nt_start")), _int(domain.get("nt_end"))
        if d_hi < d_lo:
            d_lo, d_hi = d_hi, d_lo
        if d_lo <= 0 or d_hi <= 0:
            return {}
        want = _antisense_strand(domain.get("strand", "+"))
        best, best_ov = None, 0
        for c in (cds_list or []):
            try:
                g_lo, g_hi, g_st = _int(c[0]), _int(c[1]), _int(c[2], 1)
            except Exception:
                continue
            if g_hi < g_lo:
                g_lo, g_hi = g_hi, g_lo
            if g_st != want:
                continue                                  # same-strand neighbour, not a partner
            if not (g_lo <= d_hi and g_hi >= d_lo):
                continue
            ov = min(g_hi, d_hi) - max(g_lo, d_lo) + 1
            if ov > best_ov:
                best, best_ov = c, ov
        if best is None:
            return {}
        g_lo, g_hi, g_st = _int(best[0]), _int(best[1]), _int(best[2], 1)
        if g_hi < g_lo:
            g_lo, g_hi = g_hi, g_lo
        product = str(best[3] or "") if len(best) > 3 else ""
        gene = str(best[4] or "") if len(best) > 4 else ""
        transl = str(best[5] or "") if len(best) > 5 else ""
        # fetch_ncbi currently returns 6-tuples; read positions 6/7 only if a richer tuple
        # ever arrives, so this stays correct either way (see _fetch_feature_table_ids).
        locus_tag = str(best[6] or "") if len(best) > 6 else ""
        protein_id = str(best[7] or "") if len(best) > 7 else ""
        # prefer the real translation length; fall back to the coding span minus the stop codon
        aa_len = len(transl) if transl else max(0, (g_hi - g_lo + 1) // 3 - 1)
        return {"start": g_lo, "end": g_hi, "strand": g_st, "product": product, "gene": gene,
                "locus_tag": locus_tag, "protein_id": protein_id, "aa_len": aa_len,
                "overlap_bp": best_ov,
                "nested_fully": 1 if (g_lo <= d_lo and g_hi >= d_hi) else 0,
                "category": _categorize(product, gene)}
    except Exception:
        return {}


def orf_matches_host(orf_lo: int, orf_hi: int, host: dict, tol_frac: float = 0.10):
    """Are the COMPUTED open antisense ORF and the ANNOTATED host CDS the same gene?
    Returns 1/0, or '' when either side is missing (never guess from half the evidence).

    Two conditions, both required:
      1. SAME READING FRAME. A CDS is a whole number of codons and so is the walked ORF, so
         an in-frame pair shares a codon boundary; either end may be tested (which one lines
         up exactly depends on whether the annotation includes the terminal stop codon).
      2. ESSENTIALLY THE SAME INTERVAL — the overlap covers >= (1-tol) of the CDS *and*
         >= (1-tol) of the ORF. The two-sided test matters: a stop-to-stop ORF legitimately
         starts a little upstream of the annotated ATG, but an ORF twice the length of the
         CDS is a different (or mis-annotated) gene and must not be scored as a match.

    A 1 here means an ORF derived purely from raw sequence coincides with a gene an
    independent annotator named — which is why it is the strongest column in the table. It
    still says nothing about whether the small antisense ORF is expressed.
    """
    try:
        if not host or not orf_lo or not orf_hi:
            return ""
        h_lo, h_hi = _int(host.get("start")), _int(host.get("end"))
        if h_lo <= 0 or h_hi <= 0:
            return ""
        in_frame = ((h_lo - int(orf_lo)) % 3 == 0) or ((h_hi - int(orf_hi)) % 3 == 0)
        ov = max(0, min(h_hi, int(orf_hi)) - max(h_lo, int(orf_lo)) + 1)
        h_len = max(1, h_hi - h_lo + 1)
        o_len = max(1, int(orf_hi) - int(orf_lo) + 1)
        same = ov >= (1 - tol_frac) * h_len and ov >= (1 - tol_frac) * o_len
        return 1 if (in_frame and same) else 0
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# NCBI
# ---------------------------------------------------------------------------
def _fetch_genomes(accessions: list, email: str, log=print):
    """(seqs, names, feats) from build_real_genbanks.fetch_ncbi, or three empty dicts.
    Each is keyed by accession with AND without the version suffix."""
    if not accessions or not email:
        return {}, {}, {}
    try:
        from build_real_genbanks import fetch_ncbi
        return fetch_ncbi(sorted(set(accessions)), email)
    except Exception as e:
        log(f"  (overprint: NCBI annotation fetch failed, host columns left empty: {e})")
        return {}, {}, {}


_FT_ACC = re.compile(r">Feature\s+(?:\w+\|)?([A-Za-z0-9_.]+)")


def _fetch_feature_table_ids(accessions: list, email: str, log=print) -> dict:
    """Best-effort {accession: {(lo, hi, strand): {"locus_tag": .., "protein_id": ..}}}.

    WHY A SECOND CALL: fetch_ncbi returns 6-tuples that drop /locus_tag and /protein_id, and
    that file is owned elsewhere — but those two identifiers are what let a reader pull the
    host protein straight out of NCBI, so they are recovered here. This is NOT a second
    genome download: rettype='ft' returns the bare feature table (a few kB per genome, no
    sequence) and is restricted to the handful of accessions that actually carry a host gene.
    Returns {} on any failure; the columns then stay empty and nothing else changes.
    """
    if not accessions or not email:
        return {}
    out: dict = {}
    try:
        from Bio import Entrez
        Entrez.email = email
        ids = sorted(set(a for a in accessions if a))
        h = Entrez.efetch(db="nuccore", id=",".join(ids), rettype="ft", retmode="text")
        text = h.read()
    except Exception as e:
        log(f"  (overprint: feature-table fetch failed, locus_tag/protein_id left empty: {e})")
        return {}
    try:
        acc, key, coord = "", "", None
        for ln in text.splitlines():
            if ln.startswith(">Feature"):
                m = _FT_ACC.search(ln)
                acc = m.group(1) if m else ""
                out.setdefault(acc, {})
                out.setdefault(acc.split(".")[0], out[acc])   # version-insensitive alias
                key, coord = "", None
                continue
            if not acc or not ln.strip():
                continue
            p = ln.split("\t")
            if p and p[0].strip():                       # "start<TAB>end[<TAB>key]" line
                try:
                    a = _int(re.sub(r"[<>]", "", p[0]))
                    b = _int(re.sub(r"[<>]", "", p[1]))
                except Exception:
                    continue
                if len(p) > 2 and p[2].strip():
                    key = p[2].strip()
                    coord = (min(a, b), max(a, b), -1 if a > b else 1)
                continue
            if key in ("CDS", "gene") and coord and len(p) >= 5:
                qual, val = p[3].strip(), p[4].strip()
                if qual in ("locus_tag", "protein_id") and val:
                    # protein_id arrives as gb|ANZ50912.1| — keep just the accession
                    val = val.split("|")[1] if val.count("|") >= 2 else val
                    slot = out[acc].setdefault(coord, {})
                    slot.setdefault(qual, val)           # a CDS beats a later gene feature
    except Exception as e:
        log(f"  (overprint: feature table unparsed, locus_tag/protein_id left empty: {e})")
        return {}
    return out


def _lookup_ids(ft_map: dict, acc: str, host: dict) -> dict:
    """locus_tag/protein_id for a chosen host CDS, matched on its exact interval+strand.
    Exact match only — a fuzzy match here would attach the wrong protein accession to the
    headline gene, which is worse than an empty cell."""
    try:
        table = ft_map.get(acc) or ft_map.get(str(acc).split(".")[0]) or {}
        return table.get((_int(host.get("start")), _int(host.get("end")),
                          _int(host.get("strand"), 1)), {}) or {}
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# Row assembly
# ---------------------------------------------------------------------------
def build_rows(rows: list, seqs: dict, feats: dict, log=print) -> list:
    """One record per DISTINCT interrupted locus, with locus_id OP001, OP002, … in input
    order (= descending domain bit score, so OP001 is the best-scoring locus).

    Rows are deduped FIRST (see `dedupe_rows`) so an OP id names a genomic locus and not a
    database record: the same locus reported under a RefSeq and a GenBank accession is one
    OP id with the other accessions listed in `accession_aliases`. Assigning the ids before
    deduping is what previously produced 27 ids for 13 loci.

    Host columns stay empty when no annotation is available — the row is never dropped."""
    try:
        from find_interrupted import antisense_open_extent
    except Exception as e:
        log(f"  (overprint: antisense ORF extent unavailable: {e})")
        antisense_open_extent = None
    rows = dedupe_rows(rows)          # idempotent; safe when the caller already deduped
    recs = []
    for i, r in enumerate(rows, 1):
        contig = str(r.get("contig", "") or "")
        acc = _accession_of(contig) or contig
        d_lo, d_hi = _int(r.get("domain_nt_start")), _int(r.get("domain_nt_end"))
        strand = str(r.get("strand", "") or "")
        rec = {c: "" for c in LOCI_COLS}
        rec.update({
            "locus_id": f"OP{i:03d}",
            "contig": contig,
            "organism": str(r.get("organism", "") or "") or contig,
            "accession": acc,
            # every OTHER database accession this same genome/locus was found under
            "accession_aliases": str(r.get("accession_aliases", "") or ""),
            "strand": strand,
            "frame": r.get("frame", ""),
            "domain_nt_start": d_lo, "domain_nt_end": d_hi,
            "domain_aa_len": _int(r.get("domain_aa_len")),
            "internal_stops": _int(r.get("internal_stops")),
            "stop_nt_positions": str(r.get("stop_nt_positions", "") or ""),
            "overprinting_support": str(r.get("overprinting_support", "") or "none"),
            "antisense_open_frame": r.get("antisense_open_frame", ""),
            "antisense_open_stops": r.get("antisense_open_stops", ""),
            "stop_silent_antisense": str(r.get("stop_silent_antisense", "") or ""),
            "host_annotation_source": "none",
        })
        contig_seq = seqs.get(acc) or seqs.get(acc.split(".")[0]) or ""
        anti_frame = _int(r.get("antisense_open_frame"), -1)
        anti_st = "+" if _antisense_strand(strand) > 0 else "-"
        # THE GATE. antisense_open_extent measures the extent of an ORF the domain sits
        # inside; on a locus whose antisense frame is NOT open (antisense_open_stops > 0)
        # there is no such ORF, and writing its coordinates into columns named
        # antisense_orf_* published a stop-riddled interval as "the computed open antisense
        # ORF" (one row with support=none / 7 stops emitted a 56-aa "ORF" translating to
        # GYFFIPVSDI*APPGV*HRRGHYIESMSCSL...). Belt and braces: the function now self-gates
        # too, and antisense_orf_is_open states the verdict explicitly so a reader can never
        # mistake an envelope for an ORF.
        #   1  = frame open across the domain and the ORF extent was measured
        #   0  = frame NOT open across the domain; there is no ORF and no coordinates
        #   '' = not determined (offline: no contig sequence, or no antisense frame recorded)
        # A missing or unparseable antisense_open_stops is "undetermined" (-1 -> ''), NOT
        # "closed": asserting a frame is closed on the strength of a value we could not read
        # would be the same class of mistake in the other direction.
        n_stops = _int(r.get("antisense_open_stops"), -1)
        if n_stops > 0:
            rec["antisense_orf_is_open"] = 0
        elif n_stops == 0 and contig_seq and antisense_open_extent and anti_frame in (0, 1, 2) \
                and d_lo and d_hi:
            o_lo, o_hi, o_aa = antisense_open_extent(contig_seq, anti_st, anti_frame, d_lo, d_hi)
            # o_aa == 0 means the function's own stop test failed (an edge-straddling stop
            # `antisense_open_stops` cannot see) — record that as "not open", not as blank.
            rec["antisense_orf_is_open"] = 1 if o_aa else 0
            if o_aa:
                rec["antisense_orf_nt_start"], rec["antisense_orf_nt_end"] = o_lo, o_hi
                rec["antisense_orf_aa_len"] = o_aa
        cds = feats.get(acc) or feats.get(acc.split(".")[0]) or []
        host = pick_host_gene({"nt_start": d_lo, "nt_end": d_hi, "strand": strand}, cds)
        if host:
            rec.update({
                "host_gene": host["gene"], "host_product": host["product"],
                "host_locus_tag": host["locus_tag"], "host_protein_id": host["protein_id"],
                "host_gene_start": host["start"], "host_gene_end": host["end"],
                "host_gene_strand": "+" if host["strand"] >= 0 else "-",
                "host_gene_aa_len": host["aa_len"],
                "host_annotation_source": "NCBI GenBank CDS",
                "host_category": host["category"],
                "nested_fully": host["nested_fully"],
                "overlap_bp": host["overlap_bp"],
                "overlap_pct_of_domain": (round(100.0 * host["overlap_bp"] / (d_hi - d_lo + 1), 1)
                                          if d_hi >= d_lo else ""),
            })
            rec["antisense_orf_matches_host_gene"] = orf_matches_host(
                rec["antisense_orf_nt_start"], rec["antisense_orf_nt_end"], host)
        rec["_host"] = host          # private, stripped before writing (figures need it)
        recs.append(rec)
    return recs


def write_loci_csv(recs: list, out_dir: Path, log=print) -> str:
    """overprinted_loci.csv — one row per interrupted locus. '' on failure."""
    try:
        p = Path(out_dir) / "overprinted_loci.csv"
        with open(p, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=LOCI_COLS, extrasaction="ignore")
            w.writeheader()
            w.writerows(recs)
        log(f"  overprinting: {len(recs)} locus row(s) -> {p.name}")
        return str(p)
    except Exception as e:
        log(f"  (overprint: overprinted_loci.csv not written: {e})")
        return ""


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
def _p_open_by_chance(aa_len: int) -> float:
    """P(at least one of the 3 antisense frames is stop-free across `aa_len` codons) under a
    uniform-codon null: 1 - (1 - (1 - 3/64)^L)^3.

    Deliberately computed over the DOMAIN length, not the full open-ORF length. The ORF's
    bounds were found BY walking to the flanking stops, so feeding that length back in would
    be circular and would manufacture an absurdly small p-value. The domain length is the
    span over which analyze_overprinting actually requires openness, so it is the honest L.
    """
    try:
        L = int(aa_len)
        if L <= 0:
            return 0.0
        per_frame = (1.0 - P_STOP_PER_CODON) ** L
        return 1.0 - (1.0 - per_frame) ** N_ANTISENSE_FRAMES
    except Exception:
        return 0.0


def _median(values):
    vals = [v for v in values if isinstance(v, (int, float)) and v > 0]
    return round(statistics.median(vals), 1) if vals else ""


def summarise(recs: list) -> list:
    """The long metric/value/note summary. Every note carries either the definition or the
    hedge, so the CSV is readable without the methods section.

    `recs` is the DEDUPED locus list (see `dedupe_rows`), so every count here is per genomic
    locus, never per database record — which is what makes n_interrupted_loci and
    n_distinct_organisms tell the same story instead of contradicting each other two rows
    apart."""
    from collections import Counter
    n = len(recs)
    sup = Counter(str(r.get("overprinting_support") or "none") for r in recs)
    products = [str(r.get("host_product") or "").strip() for r in recs]
    named = [p for p in products if p]
    prod_counts = Counter(named)
    top_prod, top_n = (prod_counts.most_common(1)[0] if prod_counts else ("", 0))
    dom_lens = [_int(r.get("domain_aa_len")) for r in recs]
    orf_lens = [_int(r.get("antisense_orf_aa_len")) for r in recs]
    med_dom = _median(dom_lens)
    med_orf = _median(orf_lens)
    n_match = sum(1 for r in recs if r.get("antisense_orf_matches_host_gene") == 1)
    n_open_orf = sum(1 for r in recs if r.get("antisense_orf_is_open") == 1)
    # accession + aliases = every record this locus was found in; the extra records collapsed
    n_alias = sum(len([a for a in str(r.get("accession_aliases") or "").split(";") if a.strip()])
                  for r in recs)
    orgs = {_canonical(r) for r in recs}
    p_chance = _p_open_by_chance(med_dom) if med_dom else 0.0
    orf_note = (f"median over the {n_open_orf} locus/loci with an OPEN antisense frame "
                f"(antisense_open_stops = 0); those ORFs run to a median of {med_orf} aa, far "
                f"past the domain. That extent is NOT folded into p_open_frame_by_chance (its "
                f"bounds were found by walking to stops, so doing so would be circular)"
                if med_orf else "no open antisense ORF measured (offline run, or no locus with "
                                "a fully open antisense frame)")
    return [
        {"metric": "n_interrupted_loci", "value": n,
         "note": "DISTINCT genomic loci: interrupted_homologs.tsv rows collapsed on "
                 "(canonical organism, strand, domain start, domain end), so the same locus "
                 "catalogued under a RefSeq and a GenBank accession counts once; the other "
                 "accessions are kept in the accession_aliases column"},
        {"metric": "n_accession_records_collapsed", "value": n_alias,
         "note": f"duplicate database records folded into those {n} loci "
                 f"({n + n_alias} rows in interrupted_homologs.tsv -> {n} loci); nothing is "
                 f"discarded — see accession_aliases"},
        {"metric": "n_strong", "value": sup.get("strong", 0),
         "note": f"antisense frame fully open across the domain AND every premature stop "
                 f"synonymous in it — {HEDGE}"},
        {"metric": "n_partial", "value": sup.get("partial", 0),
         "note": "some premature stops synonymous in the antisense frame, frame not fully open"},
        {"metric": "n_none", "value": sup.get("none", 0),
         "note": "no overprinting signature; the internal stop has some other explanation"},
        {"metric": "n_with_named_host_gene", "value": len(named),
         "note": "loci with an annotated antisense CDS overlapping the domain (needs --email)"},
        {"metric": "top_host_product", "value": top_prod,
         "note": "most frequent annotated host product across the interrupted loci"},
        {"metric": "n_loci_with_top_host_product", "value": top_n,
         "note": f"of {len(named)} loci with any named host gene"},
        {"metric": "n_distinct_host_products", "value": len(prod_counts),
         "note": "distinct /product strings; annotation wording varies between submitters"},
        {"metric": "n_fully_nested", "value": sum(1 for r in recs if r.get("nested_fully") == 1),
         "note": "host CDS entirely contains the domain (true nesting, not a partial overlap)"},
        {"metric": "median_domain_aa_len", "value": med_dom,
         "note": "HMM domain envelope length of the interrupted homolog"},
        {"metric": "n_with_open_antisense_orf", "value": n_open_orf,
         "note": "loci whose antisense frame is open across the whole domain AND whose ORF "
                 "extent was therefore measured; the antisense_orf_* columns are written for "
                 "these loci ONLY (antisense_orf_is_open = 1) — on a closed frame there is no "
                 "ORF, only an envelope full of stops, and it is not reported as one"},
        {"metric": "median_antisense_open_orf_aa_len", "value": med_orf,
         "note": orf_note},
        {"metric": "n_distinct_organisms", "value": len(orgs),
         "note": "phages carrying an interrupted locus, counted with canonical.canonical_organism "
                 "— the SAME identity used to collapse the loci above, so the two numbers are "
                 "consistent by construction"},
        {"metric": "p_open_frame_by_chance", "value": (f"{p_chance:.2e}" if p_chance else ""),
         "note": f"1-(1-(1-3/64)^L)^3 at the REAL median domain length L={med_dom} aa "
                 f"(previously quoted for an illustrative 137 aa); "
                 f"{n_match} locus/loci also have the computed open antisense ORF matching "
                 f"the annotated host CDS — {HEDGE}"},
    ]


def write_summary_csv(recs: list, out_dir: Path, log=print) -> str:
    """overprinting_summary.csv (metric, value, note). '' on failure."""
    try:
        p = Path(out_dir) / "overprinting_summary.csv"
        # utf-8: the notes carry em dashes, and losing the whole summary to a codec error on a
        # Windows host would be a silent, total loss of the stage's explanatory text
        with open(p, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=SUMMARY_COLS)
            w.writeheader()
            w.writerows(summarise(recs))
        log(f"  overprinting: summary -> {p.name}")
        return str(p)
    except Exception as e:
        log(f"  (overprint: overprinting_summary.csv not written: {e})")
        return ""


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _mpl():
    """The project's standard matplotlib preamble (identical to export_csv._db_barplot):
    Agg so it works headless, svg.fonttype='none' + pdf.fonttype=42 so text stays EDITABLE
    text in Illustrator/Inkscape rather than being converted to outlines. None if absent."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        matplotlib.rcParams["svg.fonttype"] = "none"
        matplotlib.rcParams["pdf.fonttype"] = 42
        import matplotlib.pyplot as plt
        return plt
    except Exception:
        return None


def overview_figure(recs: list, out_dir: Path, log=print) -> list:
    """overprinting_overview.{png,svg,pdf} — a horizontal lollipop, one line per locus,
    x = domain aa length, stem+dot coloured by overprinting support, annotated with the
    annotated host product. Reads as a family census: how many loci, how long, how well
    supported and what they are nested inside. Empty list if matplotlib is unavailable."""
    plt = _mpl()
    if plt is None or not recs:
        if plt is None:
            log("  (overprint: matplotlib unavailable, overview figure skipped)")
        return []
    try:
        order = sorted(recs, key=lambda r: (SUPPORT_ORDER.get(r.get("overprinting_support"), 3),
                                            -_int(r.get("domain_aa_len"))))
        n = len(order)
        lens = [_int(r.get("domain_aa_len")) for r in order]
        xmax = max(lens + [1])
        fig, ax = plt.subplots(figsize=(11.0, max(3.0, 0.30 * n + 2.2)))
        for i, r in enumerate(order):
            y = n - 1 - i                                   # first locus at the top
            col = SUPPORT_COLOURS.get(r.get("overprinting_support"), SUPPORT_COLOURS["none"])
            L = _int(r.get("domain_aa_len"))
            ax.plot([0, L], [y, y], color=col, lw=1.4, solid_capstyle="round", zorder=2)
            ax.plot([L], [y], marker="o", ms=6.5, color=col, zorder=3)
            prod = str(r.get("host_product") or "").strip()
            tag = prod if prod else "host gene not resolved (offline)"
            if r.get("antisense_orf_matches_host_gene") == 1:
                tag += "  ✓ ORF = CDS"          # computed ORF coincides with the annotation
            ax.text(L + xmax * 0.02, y, tag, va="center", ha="left", fontsize=6.8,
                    color=("#1a2230" if prod else "#8a939d"),
                    style=("normal" if prod else "italic"))
        ax.set_yticks(range(n))
        ax.set_yticklabels([f"{r['locus_id']}  {str(r.get('organism') or '')[:38]}"
                            for r in reversed(order)], fontsize=6.8)
        ax.set_xlim(0, xmax * 1.65)          # right-hand room for the host-product annotations
        ax.set_ylim(-1, n)
        ax.set_xlabel("interrupted domain length (aa)", fontsize=9)
        ax.set_title("Interrupted / overprinted loci and the host genes they are nested in",
                     fontsize=11, fontweight="bold")
        from matplotlib.lines import Line2D
        present = [s for s in ("strong", "partial", "none")
                   if any(r.get("overprinting_support") == s for r in order)]
        ax.legend(handles=[Line2D([0], [0], color=SUPPORT_COLOURS[s], lw=2.2, marker="o",
                                  ms=5, label=f"support: {s}") for s in present],
                  loc="lower right", fontsize=7.5, frameon=False)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="y", length=0)
        fig.text(0.01, 0.005, HEDGE_SENTENCE, fontsize=7, color="#666", ha="left")
        fig.tight_layout(rect=(0, 0.02, 1, 1))
        made = []
        for ext in ("png", "svg", "pdf"):
            p = Path(out_dir) / f"overprinting_overview.{ext}"
            fig.savefig(p, dpi=300 if ext == "png" else None, bbox_inches="tight")
            made.append(str(p))
        plt.close(fig)
        log("  overprinting: overview -> overprinting_overview.png / .svg / .pdf")
        return made
    except Exception as e:
        log(f"  (overprint: overview figure skipped: {e})")
        return []


def _stop_marks(rec: dict) -> list:
    """(nt_pos, label, colour) per premature stop, for genome_map.draw(marks=...). Only the
    FIRST stop is labelled: interrupted loci here carry stops ~12 bp apart, and on a 10 kb
    axis two labels at the same pixel are unreadable — so the count goes in the one label."""
    try:
        pos = [int(x) for x in str(rec.get("stop_nt_positions") or "").split(";") if x.strip()]
    except Exception:
        return []
    if not pos:
        return []
    head = "premature stop" + (f" (×{len(pos)})" if len(pos) > 1 else "")
    return [(p, head if i == 0 else "", STOP_COLOUR) for i, p in enumerate(pos)]


def _locus_genes(rec: dict, cds_list: list):
    """Gene list for one locus diagram, via genome_map.build_genes: the interrupted domain as
    the gold anchor, every CDS overlapping it (the host gene lands here with role 'overlap',
    which is what stacks it onto its own lane) and FLANKS CDS each side for context.
    Returns (genes, anchor) or (None, None)."""
    try:
        import genome_map as GM
        d_lo, d_hi = _int(rec.get("domain_nt_start")), _int(rec.get("domain_nt_end"))
        if d_lo <= 0 or d_hi <= 0:
            return None, None
        a_st = 1 if str(rec.get("strand")) == "+" else -1
        over = [c for c in cds_list if _int(c[0]) <= d_hi and _int(c[1]) >= d_lo]
        w_lo = min([_int(c[0]) for c in over] + [d_lo])
        w_hi = max([_int(c[1]) for c in over] + [d_hi])
        up = sorted([c for c in cds_list if _int(c[1]) < w_lo], key=lambda c: _int(c[1]))[-FLANKS:]
        down = sorted([c for c in cds_list if _int(c[0]) > w_hi], key=lambda c: _int(c[0]))[:FLANKS]
        called, flank_keys = [], set()
        for c in up + over + down:
            s, e, st = _int(c[0]), _int(c[1]), _int(c[2], 1)
            prod = str(c[3] or "") if len(c) > 3 else ""
            gene = str(c[4] or "") if len(c) > 4 else ""
            called.append((s, e, st, {"product": prod, "gene": gene,
                                      "category": _categorize(prod, gene)}))
            if c in up or c in down:
                flank_keys.add((s, e))
        anchor = (d_lo, d_hi, a_st)
        return GM.build_genes(anchor, called, flank_keys), anchor
    except Exception:
        return None, None


def _locus_caption(rec: dict) -> str:
    """One-line caption stating exactly what the reader is looking at — including whether the
    COMPUTED antisense ORF coincides with the ANNOTATED host CDS — and closing on the hedge."""
    bits = [f"{rec['locus_id']}: {rec.get('domain_aa_len')} aa domain interrupted by "
            f"{rec.get('internal_stops')} premature stop(s) (dashed red)"]
    if str(rec.get("antisense_open_stops")) != "":
        bits.append(f"antisense frame {rec.get('antisense_open_frame')} carries "
                    f"{rec.get('antisense_open_stops')} stop(s) across the domain")
    if rec.get("host_product"):
        bits.append(f"host CDS: {rec['host_product']} "
                    f"({rec.get('host_gene_aa_len')} aa, {rec.get('host_gene_strand')} strand)")
    if rec.get("antisense_orf_aa_len"):
        m = rec.get("antisense_orf_matches_host_gene")
        verdict = ("matches the annotated host CDS" if m == 1
                   else ("does NOT match the annotated host CDS" if m == 0 else "unmatched"))
        bits.append(f"computed open antisense ORF {rec['antisense_orf_nt_start']}–"
                    f"{rec['antisense_orf_nt_end']} ({rec['antisense_orf_aa_len']} aa) {verdict}")
    elif rec.get("antisense_orf_is_open") == 0:
        # Say it out loud rather than going quiet: a caption that simply omits the ORF reads
        # like a rendering gap, whereas "the frame is closed" is the actual finding.
        bits.append("the antisense frame is NOT open across the domain, so there is no "
                    "antisense ORF to report for this locus")
    return "; ".join(bits) + f".\n{HEDGE_SENTENCE}"


def locus_figures(recs: list, feats: dict, out_dir: Path, log=print,
                  cap: int = MAX_LOCUS_FIGURES) -> int:
    """Per-locus antisense diagram (DNA Features Viewer) for up to `cap` loci: the gold
    interrupted domain, the antisense host gene stacked on its own lane, the flanking genes,
    and a dashed red tick at every premature stop. Sets rec['figure'] to the path relative to
    `out_dir`. Returns how many were drawn; 0 (with a logged reason) offline."""
    try:
        import genome_map as GM
    except Exception as e:
        log(f"  (overprint: genome_map unavailable, per-locus diagrams skipped: {e})")
        return 0
    # genome_map only sets Agg; applying the project's full preamble here as well means these
    # locus SVG/PDFs open with LIVE, editable text in Illustrator/Inkscape like every other
    # publication figure the pipeline emits (the rcParams are the same values export_csv sets).
    _mpl()
    fig_dir = Path(out_dir) / "overprinting_loci"
    drawn, skipped = 0, 0
    for rec in recs:
        if drawn >= cap:
            skipped += 1
            continue
        acc = str(rec.get("accession") or "")
        cds = feats.get(acc) or feats.get(acc.split(".")[0]) or []
        if not cds:
            skipped += 1
            continue                              # offline, or a genome with no CDS annotation
        genes, anchor = _locus_genes(rec, cds)
        if not genes:
            skipped += 1
            continue
        try:
            fig_dir.mkdir(parents=True, exist_ok=True)
            base = fig_dir / f"{rec['locus_id']}_{_safe(rec.get('organism'))}_{_safe(acc, 20)}"
            got = GM.draw(genes, anchor, base, _locus_caption(rec), log=lambda *_a, **_k: None,
                          track_name=f"{rec.get('organism')}\n{acc}", tool="dfv",
                          labels=True, marks=_stop_marks(rec))
            if got:
                rec["figure"] = f"overprinting_loci/{Path(str(got)).name}.png"
                drawn += 1
            else:
                skipped += 1
        except Exception as e:
            skipped += 1
            log(f"  (overprint: diagram for {rec.get('locus_id')} skipped: {e})")
    if drawn:
        log(f"  overprinting: {drawn} per-locus diagram(s) -> {fig_dir.name}/"
            + (f" ({skipped} not drawn; cap {cap})" if skipped else ""))
    elif skipped:
        log(f"  (overprint: no per-locus diagrams — no CDS annotation available "
            f"for {skipped} locus/loci; re-run with --email)")
    return drawn


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def run(discovery_dir, out_dir=None, email: str = "", cpu: int = 1, log=print) -> dict:
    """Build the whole overprinting deliverable for a discovery directory. Always returns a
    dict (never raises); "loci" is 0 when there is nothing to report."""
    discovery = Path(discovery_dir)
    out = Path(out_dir) if out_dir else discovery
    try:
        out.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        log(f"  (overprint: cannot write to {out}: {e})")
        return {"loci": 0}
    rows = read_interrupted(discovery / "interrupted_homologs.tsv")
    if not rows:
        log("  (overprint: no interrupted_homologs.tsv rows — nothing to report)")
        return {"loci": 0}
    # Dedupe BEFORE anything downstream: the OP ids, the summary counts and the NCBI fetch
    # list must all be per-LOCUS. (build_rows dedupes again; it is idempotent.)
    n_raw = len(rows)
    rows = dedupe_rows(rows)
    if len(rows) != n_raw:
        log(f"  overprinting: {n_raw} interrupted row(s) -> {len(rows)} distinct locus/loci "
            f"({n_raw - len(rows)} cross-database accession alias(es) collapsed; "
            f"kept in accession_aliases)")
    accs = sorted({_primary_accession(r) for r in rows if r.get("contig")})
    seqs, _names, feats = _fetch_genomes(accs, email, log)
    if not email:
        log("  (overprint: no --email — running OFFLINE; host gene columns stay empty)")
    recs = build_rows(rows, seqs, feats, log)
    # locus_tag / protein_id come from a tiny feature-table pass over only the accessions that
    # actually produced a host gene (see _fetch_feature_table_ids for why it is a second call)
    host_accs = sorted({r["accession"] for r in recs if r.get("host_gene_start")})
    ft = _fetch_feature_table_ids(host_accs, email, log) if host_accs else {}
    if ft:
        for r in recs:
            ids = _lookup_ids(ft, r.get("accession", ""), r.get("_host") or {})
            r["host_locus_tag"] = r["host_locus_tag"] or ids.get("locus_tag", "")
            r["host_protein_id"] = r["host_protein_id"] or ids.get("protein_id", "")
    # FIGURES go to downstream/overprinting/, which is the directory assemble_package copies
    # into PACKAGE/10_overprinting. They were previously written beside the CSVs at the run
    # root, so that copy had no source and the antisense diagrams — the evidence for the most
    # interesting claim in the run — never reached the shareable bundle at all.
    # The two CSVs stay at the run root, where export_csv's TABLE_EXPORTS mirror picks them up
    # into 01_summary_tables. Both destinations are deliberate; neither is a fallback.
    fig_out = out / "downstream" / "overprinting" if out_dir is None else out
    try:
        fig_out.mkdir(parents=True, exist_ok=True)
    except OSError:
        fig_out = out
    n_fig = locus_figures(recs, feats, fig_out, log)
    figs = overview_figure(recs, fig_out, log)
    for r in recs:
        r.pop("_host", None)                      # private helper field, not a report column
    loci_csv = write_loci_csv(recs, out, log)
    summary_csv = write_summary_csv(recs, out, log)
    named = sum(1 for r in recs if str(r.get("host_product") or "").strip())
    matched = sum(1 for r in recs if r.get("antisense_orf_matches_host_gene") == 1)
    open_orf = sum(1 for r in recs if r.get("antisense_orf_is_open") == 1)
    log(f"  overprinting: {len(recs)} locus/loci, {named} with a named host gene, "
        f"{matched} with the computed antisense ORF matching that host CDS")
    return {"loci": len(recs), "named_host": named, "orf_matches_host": matched,
            "rows_in": n_raw, "aliases_collapsed": n_raw - len(rows),
            "open_antisense_orfs": open_orf,
            "locus_figures": n_fig, "loci_csv": loci_csv, "summary_csv": summary_csv,
            "overview": figs, "out_dir": str(out), "figure_dir": str(fig_out)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discovery-dir", type=Path, required=True,
                    help="a discovery run directory containing interrupted_homologs.tsv")
    ap.add_argument("--email", default="",
                    help="NCBI Entrez email; WITHOUT it the run is offline and the host-gene "
                         "columns are left empty (every table is still written)")
    ap.add_argument("--cpu", type=int, default=1,
                    help="accepted for interface parity with the other stage scripts; this "
                         "stage is network- and render-bound and runs single-threaded")
    ap.add_argument("--out", type=Path, default=None,
                    help="output directory (default: the discovery directory)")
    args = ap.parse_args()
    if not args.discovery_dir.exists():
        sys.exit(f"discovery directory not found: {args.discovery_dir}")
    s = run(args.discovery_dir, args.out, args.email, args.cpu)
    print(f"  overprinted loci: {s.get('loci', 0)}; named host gene: {s.get('named_host', 0)}; "
          f"ORF matches host CDS: {s.get('orf_matches_host', 0)}")


if __name__ == "__main__":
    main()
