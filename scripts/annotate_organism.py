#!/usr/bin/env python3
"""
annotate_organism.py — add an `organism` column to a hits table
====================================================================
Looks up the source organism (phage name) for each hit's genome and inserts an
`organism` column immediately after `genome_id` in a hits.tsv.

  - NCBI-accession genomes (INPHARED / RefSeq): the phage name is fetched from
    NCBI (Entrez esummary), e.g. "Escherichia phage vB_EcoP_G7C".
  - Metagenomic genomes (GVD-AVrC / GPD): these are uncultured and unclassified,
    so the organism is recorded as "uncultured virus (<source database>)".

Names are cached so the same genome is only queried once. Safe to re-run.

USAGE
-----
  python3 annotate_organism.py --hits-tsv path/to/hits.tsv
  python3 annotate_organism.py --hits-tsv a.tsv b.tsv c.tsv   # multiple at once
"""
from __future__ import annotations

import argparse
import re
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez


def is_ncbi(s: str) -> bool:
    return bool(re.match(r"^[A-Z]{1,2}_?\d{5,8}", str(s)))


def clean_title(t: str) -> str:
    """'Escherichia phage X, complete genome' -> 'Escherichia phage X'."""
    return re.split(r",\s*(complete|partial|whole|genome assembly|DNA)\b", t, 1)[0].strip().rstrip(",").strip()


def fetch_names(accessions: list[str], email: str) -> dict[str, str]:
    """Return {accession: phage_name} from NCBI (with/without version key)."""
    Entrez.email = email
    names: dict[str, str] = {}
    for i in range(0, len(accessions), 40):
        chunk = accessions[i:i + 40]
        for attempt in range(3):
            try:
                h = Entrez.esummary(db="nuccore", id=",".join(chunk))
                for rec in Entrez.read(h):
                    acc = rec.get("AccessionVersion", "")
                    nm = clean_title(rec.get("Title", ""))
                    if acc and nm:
                        names[acc] = nm
                        names[acc.split(".")[0]] = nm
                break
            except Exception as e:
                print(f"  esummary retry {attempt+1}: {e}")
                time.sleep(4)
        time.sleep(0.4)
    return names


def annotate(tsv: Path, email: str) -> None:
    df = pd.read_csv(tsv, sep="\t")
    if df.empty or "genome_id" not in df.columns:
        # 0-hit table: just ensure an organism column exists, then return.
        if "organism" not in df.columns:
            df["organism"] = []
        df.to_csv(tsv, sep="\t", index=False)
        print(f"  {tsv.name}: 0 rows — organism column ensured.")
        return
    accs = sorted({str(g) for g in df["genome_id"] if is_ncbi(g)})
    names = fetch_names(accs, email) if accs else {}

    def org(row) -> str:
        gid = str(row["genome_id"])
        nm = names.get(gid) or names.get(gid.split(".")[0])
        if nm:
            return nm
        return f"uncultured virus ({row.get('db_name', 'metagenomic')})"

    df["organism"] = df.apply(org, axis=1)
    # place `organism` right after `genome_id`
    cols = list(df.columns)
    cols.remove("organism")
    idx = cols.index("genome_id") + 1
    cols = cols[:idx] + ["organism"] + cols[idx:]
    df = df[cols]
    df.to_csv(tsv, sep="\t", index=False)
    print(f"  {tsv.name}: organism added ({len(accs)} NCBI names; {len(df)} rows)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits-tsv", type=Path, nargs="+", required=True)
    ap.add_argument("--email", default="researcher@example.com")
    args = ap.parse_args()
    for tsv in args.hits_tsv:
        if tsv.exists():
            annotate(tsv, args.email)


if __name__ == "__main__":
    main()
