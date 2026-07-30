#!/usr/bin/env python3
"""missed_seed_report.py — plan and report the hunt for seeds the search never re-found.

WHY THIS EXISTS. A discovery run reports the homologs it found in the databases it searched.
Some input seeds are not among them, and the run says nothing about why — which reads as a
sensitivity failure when it is usually a COVERAGE one. On the reference run every unrecovered
seed is either a GenBank-only accession newer than the searched INPHARED snapshot, or a
metagenomic/prophage record that a viral-genome database cannot contain. Either way the answer
is not "the model missed it" but "that genome was never searched".

So: take each missed seed, go and fetch ITS OWN source genome, scan the model against it, and
report a per-seed verdict. `explains_miss` is the column that answers the question.

NOTHING IS CACHED. Genomes are fetched in batches, scanned, and deleted by
`scan_genome_collection.sh`; peak disk is one batch. The database cache is never touched, and
this deliberately does NOT add a database — the user's constraint is to stream, not download.

Two modes, so the shell driver stays thin and this stays testable:

  --plan    read seed_qc/seed_status.csv, work out what to fetch, write the source list
  --report  join the collection scan's hits back onto the seeds, write missed_seed_scan.csv

THE ADDRESS IS NEVER WRITTEN TO THE LIST. NCBI wants an email in the URL, but the list file is
a shipped artefact, so the literal token ``__EMAIL__`` is written instead and substituted from
$NCBI_EMAIL by `_fetch` at request time.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "engine"))

EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&id={ids}&rettype=fasta&retmode=text"
          "&tool=hmm-homologue-finder&email=__EMAIL__")
IDS_PER_REQUEST = 50          # NCBI is fine with a few hundred; 50 keeps a failure small
EMAIL_PLACEHOLDER = "__EMAIL__"

REPORT_COLS = [
    "seed_id", "organism", "accession", "accession_prefix", "accession_class",
    "source", "fetched", "n_contigs", "best_bit", "n_clean", "n_interrupted",
    "strand", "frame", "domain_nt_start", "domain_nt_end", "overprinting_support",
    "verdict", "explains_miss",
]


def _rows(p: Path) -> list:
    try:
        with open(p, encoding="utf-8", errors="replace", newline="") as fh:
            return list(csv.DictReader(fh))
    except OSError:
        return []


def _tsv_rows(p: Path) -> list:
    try:
        with open(p, encoding="utf-8", errors="replace", newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    except OSError:
        return []


def _base(acc: str) -> str:
    a = str(acc or "")
    return a[: a.rfind(".")] if "." in a and a[a.rfind(".") + 1:].isdigit() else a


def missed_seeds(run_dir: Path) -> list:
    """The seeds to chase: status == 'missed' in seed_qc/seed_status.csv.

    That file is written by family_census and is the single source of truth for which seeds
    came back, so this cannot drift from the census the report quotes.
    """
    p = Path(run_dir) / "seed_qc" / "seed_status.csv"
    return [r for r in _rows(p) if str(r.get("status", "")).strip().lower() == "missed"]


def plan(run_dir: Path, out_dir: Path, log=print) -> dict:
    """Decide what to fetch and write the source list. Returns a summary dict.

    Splits the missed seeds by where their genome actually lives:
      * an NCBI nucleotide accession -> one efetch URL per batch of ids
      * a metagenomic catalogue id   -> stream the catalogue and keep ONLY those contigs
      * no resolvable accession      -> reported, never silently dropped
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seeds = missed_seeds(run_dir)
    if not seeds:
        log("  no missed seeds to chase (every distinct seed protein was re-found)")
        (out_dir / "missed_seed_sources.txt").write_text("")
        return {"n_missed": 0, "n_ncbi": 0, "n_catalogue": 0, "n_unresolved": 0}

    # Route by the FORM OF THE ACCESSION, not by accession_class. accession_class is derived
    # from the organism name, and a `MAG: Caudoviricetes sp. ... OR222882.1` seed is classed
    # metagenome while carrying an ordinary NCBI accession that efetch serves perfectly well.
    # Routing those to the GPD/GVD catalogues found no prefix match and reported them
    # "unfetchable" without ever trying — so only a genuine catalogue id goes down that path.
    try:
        from build_real_genbanks import CATALOGUES as _CATS
        cat_prefixes = tuple(_CATS)
    except Exception:
        cat_prefixes = ("uvig_", "GutCatV1_")
    ncbi, catalogue, unresolved = [], [], []
    for r in seeds:
        acc = (r.get("accession") or "").strip()
        if not acc:
            unresolved.append(r)
        elif acc.startswith(cat_prefixes):
            catalogue.append(r)
        else:
            ncbi.append(r)

    lines = []
    for i in range(0, len(ncbi), IDS_PER_REQUEST):
        ids = ",".join((r.get("accession") or "").strip() for r in ncbi[i:i + IDS_PER_REQUEST])
        lines.append(EFETCH.format(ids=ids))

    # Metagenomic contigs: stream the catalogue, keep only the wanted ids, never store the
    # catalogue itself. Written as a local path so the collection scanner's `cat` branch reads it.
    n_cat_pulled, got = 0, {}
    if catalogue:
        try:
            from build_real_genbanks import CATALOGUES, fetch_catalogue
            wanted = {(r.get("accession") or "").strip() for r in catalogue}
            wanted |= {_base(a) for a in wanted}
            for prefix, url in CATALOGUES.items():
                ids = {a for a in wanted if a.startswith(prefix)}
                if not ids:
                    continue
                log(f"  streaming {prefix}* catalogue for {len(ids)} contig(s) "
                    f"(the catalogue itself is NOT stored)")
                got.update(fetch_catalogue(url, ids) or {})
            if got:
                mag = out_dir / "_missed_seed_metagenome.fna"
                with mag.open("w", encoding="utf-8") as fh:
                    for acc, seq in got.items():
                        fh.write(f">{acc}\n{seq}\n")
                lines.append(str(mag))
                n_cat_pulled = len(got)
            else:
                log("  (no metagenomic contigs retrieved — their catalogue ids may not be in "
                    "GPD/GVD; they are still reported, with source=none)")
        except Exception as e:
            log(f"  (metagenomic catalogue pull skipped: {e})")

    lst = out_dir / "missed_seed_sources.txt"
    lst.write_text("\n".join(lines) + ("\n" if lines else ""))
    # Belt and braces: the address must never reach a shipped file.
    body = lst.read_text(encoding="utf-8")
    assert "@" not in body, "refusing to write an address into the source list"

    # Record WHICH accessions the plan actually put into the source list. Without this the
    # report cannot distinguish "we looked and the gene is not there" from "we never looked":
    # a metagenomic id whose catalogue pull returned nothing would otherwise be published as
    # gene_absent, which is a confident claim from a step that never ran.
    planned = {(r.get("accession") or "").strip() for r in ncbi} | set(got)
    try:
        import json
        (out_dir / "_missed_seed_planned.json").write_text(
            json.dumps({"planned_accessions": sorted(x for x in planned if x)}, indent=2))
    except Exception:
        pass
    log(f"  plan: {len(seeds)} missed seed(s) -> {len(ncbi)} via NCBI efetch "
        f"({len(lines) - (1 if n_cat_pulled else 0)} request(s)), "
        f"{n_cat_pulled} metagenomic contig(s), {len(unresolved)} with no accession")
    log(f"  -> {lst}")
    return {"n_missed": len(seeds), "n_ncbi": len(ncbi), "n_catalogue": n_cat_pulled,
            "n_unresolved": len(unresolved), "list": str(lst),
            "n_requests": len(lines)}


def _f(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def report(run_dir: Path, out_dir: Path, log=print) -> dict:
    """Join the collection scan's hits back onto the missed seeds and write the verdict table."""
    out_dir = Path(out_dir)
    seeds = missed_seeds(run_dir)
    hits = _tsv_rows(out_dir / "results" / "collection_hits.tsv")

    # index hits by base accession of the contig they landed on
    by_acc = {}
    for h in hits:
        acc = _base((h.get("contig") or "").strip())
        if acc:
            by_acc.setdefault(acc, []).append(h)

    # Which accessions the plan actually managed to put in front of the scanner. A seed whose
    # source was never materialised (e.g. a metagenomic id its catalogue did not contain) must
    # be reported as NOT FETCHED, never as gene_absent — "the gene is not there" is a claim we
    # have not earned if we never looked.
    planned = set()
    try:
        import json
        planned = set(json.loads(
            (out_dir / "_missed_seed_planned.json").read_text(encoding="utf-8")
        ).get("planned_accessions", []))
    except Exception:
        pass
    planned |= {_base(a) for a in planned}
    scan_ran = (out_dir / "results" / "collection_hits.tsv").exists()
    rows = []
    for r in seeds:
        acc = (r.get("accession") or "").strip()
        cls = (r.get("accession_class") or "").strip().lower()
        hs = by_acc.get(_base(acc), [])
        best = max(hs, key=lambda h: _f(h.get("domain_bit_score")), default=None)
        n_clean = sum(1 for h in hs if "interrupt" not in str(h.get("status", "")).lower())
        n_int = len(hs) - n_clean
        if not acc:
            verdict, why, src = "not_fetched", "unfetchable", "none"
        elif cls == "metagenome":
            src = "GPD/GVD-AVrC"
        else:
            src = "ncbi_efetch"
        was_fetched = bool(hs) or (_base(acc) in planned or acc in planned)
        if acc:
            if hs:
                verdict = "present_interrupted" if (n_int and not n_clean) else "present_clean"
                # The model DOES find it in its own genome -- so the miss was database
                # coverage, not sensitivity. This is the answer the user wanted.
                why = "interrupted_only" if verdict == "present_interrupted" \
                    else "genome_not_in_searched_dbs"
            elif was_fetched and scan_ran:
                verdict, why = "not_detected", "gene_absent"
            else:
                # never put in front of the scanner: its catalogue did not contain the id, or
                # the fetch could not be planned at all
                verdict, why = "not_fetched", "unfetchable"
        rows.append({
            "seed_id": r.get("seed_id", ""), "organism": r.get("organism", ""),
            "accession": acc, "accession_prefix": r.get("accession_prefix", ""),
            "accession_class": r.get("accession_class", ""),
            "source": (src if (acc and was_fetched) else "none"),
            "fetched": "True" if was_fetched else "False",
            "n_contigs": len({h.get("contig") for h in hs}),
            "best_bit": "" if best is None else round(_f(best.get("domain_bit_score")), 1),
            "n_clean": n_clean, "n_interrupted": n_int,
            "strand": (best or {}).get("strand", ""), "frame": (best or {}).get("frame", ""),
            "domain_nt_start": (best or {}).get("nt_start", ""),
            "domain_nt_end": (best or {}).get("nt_end", ""),
            "overprinting_support": (best or {}).get("overprinting_support", ""),
            "verdict": verdict, "explains_miss": why,
        })

    dst = out_dir / "missed_seed_scan.csv"
    with dst.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REPORT_COLS)
        w.writeheader()
        w.writerows(rows)
    # Also drop it beside the census inputs so it is packaged with seed_qc/.
    try:
        sq = Path(run_dir) / "seed_qc"
        if sq.is_dir():
            (sq / "missed_seed_scan.csv").write_text(dst.read_text(encoding="utf-8"),
                                                     encoding="utf-8")
    except OSError:
        pass

    tally = {}
    for r in rows:
        tally[r["explains_miss"]] = tally.get(r["explains_miss"], 0) + 1
    log(f"  -> {dst}  ({len(rows)} missed seed(s))")
    for k, v in sorted(tally.items(), key=lambda kv: -kv[1]):
        log(f"     {k}: {v}")
    return {"n_missed": len(rows), "explains_miss": tally, "csv": str(dst)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--plan", action="store_true")
    g.add_argument("--report", action="store_true")
    a = ap.parse_args()
    r = plan(a.run_dir, a.out) if a.plan else report(a.run_dir, a.out)
    # a plan with nothing to fetch is a valid outcome, not an error
    return 0 if r is not None else 1


if __name__ == "__main__":
    sys.exit(main())
