#!/usr/bin/env python3
"""
canonical.py — single source of truth for collapsing the SAME phage that appears
under different database accessions (cross-database / cross-accession redundancy).

Used everywhere a count of "unique organisms" or a per-organism dedup matters, so
the number of distinct organisms is consistent across every table, figure, and the
tree — preventing the same phage from being counted several times just because it
is catalogued in INPHARED *and* RefSeq (etc.) under different accessions.
"""
from __future__ import annotations

import re


def canonical_organism(name: str, fallback: str = "") -> str:
    """Canonical phage identity for counting/deduplication.

    - Collapses host-genus aliases of the same phage:
      'Enterobacteria phage N4' and 'Escherichia phage N4'  ->  'n4'
    - Unnamed / metagenomic entries ('uncultured virus', blank) are NOT collapsed —
      they fall back to the genome accession so each distinct genome still counts once.
    - Returns '' only when neither a name nor a fallback accession is available.
    """
    s = re.sub(r"^(UNVERIFIED:?|MAG:?|TPA(?:_asm)?:?)\s*", "", str(name or "").strip(), flags=re.I).strip()
    if (not s) or re.search(r"uncultured|unclassified|metagenom|environmental", s, re.I):
        return (str(fallback).strip() or s).lower()
    m = re.search(r"\b(?:phage|virus)\b\s+(.+)$", s, flags=re.I)
    key = m.group(1).strip() if (m and m.group(1).strip()) else s
    return re.sub(r"\s+", " ", key).strip().lower()
