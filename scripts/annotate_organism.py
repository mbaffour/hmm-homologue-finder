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
import os
import re
import socket
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez

# Never let a stalled NCBI connection freeze an unattended run.
socket.setdefaulttimeout(60)


def is_ncbi(s: str) -> bool:
    return bool(re.match(r"^[A-Z]{1,2}_?\d{5,8}", str(s)))


def clean_title(t: str) -> str:
    """'Escherichia phage X, complete genome' -> 'Escherichia phage X'."""
    return re.split(r",\s*(complete|partial|whole|genome assembly|DNA)\b", t, 1)[0].strip().rstrip(",").strip()


_PROTEIN_PREFIXES = ("NP_", "YP_", "WP_", "XP_", "AP_", "QNP", "ADV", "AGT")


def _is_protein_acc(acc: str) -> bool:
    return any(str(acc).upper().startswith(p) for p in _PROTEIN_PREFIXES)


def _org_from_title(title: str) -> str:
    """Protein esummary Title is '<product> [<organism>]' -> the organism."""
    m = re.search(r"\[([^\]]+)\]\s*$", title or "")
    return m.group(1).strip() if m else ""


def _esummary_into(ids: list[str], db: str, namer, names: dict) -> None:
    """esummary `ids` against `db`; on batch error, resolve each id individually
    so one bad accession can't poison the whole batch. Quiet (no retry spam)."""
    for i in range(0, len(ids), 40):
        chunk = ids[i:i + 40]
        groups = [chunk]
        try:
            h = Entrez.esummary(db=db, id=",".join(chunk))
            for rec in Entrez.read(h):
                acc = rec.get("AccessionVersion") or rec.get("Caption", "")
                nm = namer(rec)
                if acc and nm:
                    names[acc] = nm
                    names[str(acc).split(".")[0]] = nm
            groups = []
        except Exception:
            groups = [[one] for one in chunk]  # fall back to per-accession
        for one in (g for grp in groups for g in grp):
            try:
                h = Entrez.esummary(db=db, id=one)
                for rec in Entrez.read(h):
                    acc = rec.get("AccessionVersion") or rec.get("Caption", "")
                    nm = namer(rec)
                    if acc and nm:
                        names[acc] = nm
                        names[str(acc).split(".")[0]] = nm
            except Exception:
                pass
            time.sleep(0.34)
        time.sleep(0.34)


def fetch_names(accessions: list[str], email: str) -> dict[str, str]:
    """Return {accession: name} from NCBI. Nucleotide accessions resolve to the
    genome title (phage name); protein accessions (YP_/NP_/WP_/…) are routed to
    the protein DB and resolve to their parent organism — so protein-DB hits and
    genome hits both get a proper organism, and a protein id never breaks a
    nuccore batch."""
    Entrez.email = email
    names: dict[str, str] = {}
    nuc = [a for a in accessions if not _is_protein_acc(a)]
    prot = [a for a in accessions if _is_protein_acc(a)]
    if nuc:
        _esummary_into(nuc, "nuccore", lambda r: clean_title(r.get("Title", "")), names)
    if prot:
        _esummary_into(prot, "protein", lambda r: _org_from_title(r.get("Title", "")), names)
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
    # Without a real NCBI email we must NOT contact Entrez (NCBI policy + project
    # rule: never send a placeholder address). Skip the lookup and fall back to the
    # generic organism label so the table is still well-formed offline.
    if accs and not email:
        print(f"  {tsv.name}: no --email/$NCBI_EMAIL — skipping NCBI organism lookup "
              f"({len(accs)} accession(s)); using generic labels.")
    names = fetch_names(accs, email) if (accs and email) else {}

    def org(row) -> str:
        gid = str(row["genome_id"])
        nm = names.get(gid) or names.get(gid.split(".")[0])
        if nm:
            return nm
        # No resolved name. Do NOT call a genome "uncultured virus" just because the lookup
        # failed or was skipped (offline): an NCBI accession (e.g. a cultured RefSeq NC_… record)
        # is a real cultured genome, so fall back to its accession. Reserve the metagenomic label
        # for genuinely non-NCBI ids.
        if is_ncbi(gid):
            return gid.split(".")[0]
        return f"uncultured virus ({row.get('db_name', 'metagenomic')})"

    df["organism"] = df.apply(org, axis=1)
    # place `organism` right after `genome_id`
    cols = list(df.columns)
    cols.remove("organism")
    idx = cols.index("genome_id") + 1
    cols = cols[:idx] + ["organism"] + cols[idx:]
    df = df[cols]
    df.to_csv(tsv, sep="\t", index=False)
    print(f"  {tsv.name}: organism added ({len(names)} NCBI names; {len(df)} rows)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits-tsv", type=Path, nargs="+", required=True)
    ap.add_argument("--email", default=None,
                    help="NCBI Entrez email for organism-name lookups. Never assumed: if "
                         "omitted (and $NCBI_EMAIL unset) the run stays offline and uses "
                         "generic organism labels — no address is ever sent to NCBI.")
    args = ap.parse_args()
    email = args.email or (os.environ.get("NCBI_EMAIL") or "").strip() or None
    for tsv in args.hits_tsv:
        if tsv.exists():
            annotate(tsv, email)


if __name__ == "__main__":
    main()
