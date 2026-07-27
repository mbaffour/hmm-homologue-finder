"""Homolog-locus identity and canonical-run selection — the single source of truth
shared by `hmm_finder.py` and `export_csv.py`.

Two rules live here, both of which were previously duplicated and had drifted apart:

1. **What makes two hits the same homolog** (`locus_ids`). Identity is the genomic
   locus, not the amino-acid string. `aa_sequence` is the HMM *envelope* slice, so the
   same gene is trimmed differently whenever the model is refined between iterations —
   one Erwinia gene was reported as a 63 aa, a 67 aa *and* a 118 aa "unique homolog".
   Coordinates do not drift, so the locus is the stable identity.

2. **Which iteration the package describes** (`best_run_index`). `hmm_finder` used to
   pick the run with the most hit ROWS (ties -> earliest) while `export_csv` picked the
   most unique SEQUENCES (ties -> latest). When two iterations tie on rows those rules
   disagree, and the package then shipped a `profile.hmm`, controls, tree and figures
   from one iteration alongside tables and a report headline from another — so the
   quoted sensitivity/specificity did not describe the model that produced the hits.
   One rule, applied by both, prevents that.
"""

from __future__ import annotations

import csv
import sys

from canonical import canonical_organism

__all__ = ["locus_ids", "n_unique_loci", "best_run_index", "read_hits_rows"]


def _interval(row):
    """(lo, hi) of the hit's ORF, preferring the full ORF over the domain envelope."""
    for s_col, e_col in (("orf_nt_start", "orf_nt_end"), ("nt_start", "nt_end")):
        try:
            a, b = int(float(row.get(s_col))), int(float(row.get(e_col)))
            return (min(a, b), max(a, b))
        except (TypeError, ValueError):
            continue
    return None


def _base_acc(gid: str) -> str:
    """Drop a version suffix so NC_023589.1 and NC_023589 are one genome."""
    g = str(gid or "")
    return g[: g.rfind(".")] if "." in g and g[g.rfind(".") + 1:].isdigit() else g


def locus_ids(rows) -> list:
    """Label each hit row with its physical genomic locus.

    Two rows are the same locus when they share a canonical organism and a strand and
    their ORF intervals overlap:

    * canonical organism collapses a genome catalogued under several accessions —
      notably a GenBank original and its RefSeq ``NC_`` mirror. These are the SAME
      physical sequence, so counting them as two databases overstates corroboration;
    * overlap rather than exact equality absorbs the few-nt start drift between rounds;
    * two *different* genes on one strand cannot overlap, so overlap is sufficient. The
      family's antisense partner sits on the opposite strand and stays separate.

    Protein-database hits carry no genomic coordinates and fall back to exact sequence
    identity. Returns one label per row, in input order.
    """
    rows = list(rows)
    lo, hi, key = [], [], []
    for r in rows:
        iv = _interval(r)
        lo.append(iv[0] if iv else None)
        hi.append(iv[1] if iv else None)
        org = canonical_organism(r.get("organism", ""), r.get("genome_id", ""))
        key.append((org or _base_acc(r.get("genome_id", "")), str(r.get("strand", ""))))

    buckets = {}
    for i, k in enumerate(key):
        buckets.setdefault(k, []).append(i)

    ids, n = [None] * len(rows), 0
    for _k, idxs in buckets.items():
        cur = cur_end = None
        for i in sorted((j for j in idxs if lo[j] is not None), key=lambda j: (lo[j], hi[j])):
            if cur_end is None or lo[i] > cur_end:      # disjoint -> a different gene
                n += 1
                cur, cur_end = n, hi[i]
            else:                                        # overlaps the open locus -> same gene
                cur_end = max(cur_end, hi[i])
            ids[i] = f"L{cur:04d}"
        seqmap = {}
        for i in (j for j in idxs if lo[j] is None):     # protein-DB hit: no coordinates
            s = str(rows[i].get("aa_sequence", ""))
            if s not in seqmap:
                n += 1
                seqmap[s] = n
            ids[i] = f"L{seqmap[s]:04d}"
    return ids


def n_unique_loci(rows) -> int:
    """How many distinct homolog loci a set of hit rows represents."""
    return len(set(locus_ids(rows)))


def read_hits_rows(tsv_path) -> list:
    """Read a validated `hits.tsv` into row dicts ([] if unreadable)."""
    try:
        csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    except (OverflowError, ValueError):
        pass
    try:
        with open(tsv_path, "r", encoding="utf-8", errors="replace", newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    except OSError:
        return []


def best_run_index(run_rows: dict) -> int:
    """Pick the iteration the whole package should describe.

    `run_rows` maps run index -> hit rows. The winner is the run recovering the most
    distinct homolog *loci*; ties resolve to the LATEST such run, i.e. the converged
    model, which is the one worth depositing. Falls back to 1.
    """
    best_i, best_n = None, -1
    for i in sorted(run_rows):
        n = n_unique_loci(run_rows[i])
        if n >= best_n:          # >= so a tie hands it to the later, converged round
            best_i, best_n = i, n
    return best_i if best_i is not None else 1
