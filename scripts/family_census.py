#!/usr/bin/env python3
"""
family_census.py — did the protein family actually GROW?
========================================================

The pipeline's headline used to read "55 homologs from 101 seeds", which sounds like
a LOSS. It is not a loss; it is two different quantities counted two different ways
and then subtracted, which is meaningless:

  * **101** counts seed FASTA *records*. 22 of them are byte-identical duplicates, so
    the real input is **79 distinct proteins**. Clustering before this dedup makes the
    denominator 101 and every downstream percentage wrong.
  * **55** counts distinct homolog *proteins* recovered from raw genome sequence by
    the final profile HMM. It comes from `export_csv._dedup_hits` on the canonical
    run, where homolog identity is the genomic LOCUS, never the amino-acid string
    (the string is the HMM envelope slice and is re-trimmed between iterations —
    counting strings once inflated 55 into 71).
  * **Neither number is the family.** The family is the UNION of the two, clustered
    at a stated identity. That is the number a reader wants, and it is the only one
    that answers "did my family grow?".

On the gp75 reference run the answer is unambiguous:

    identity   family before   family after   new
    100 %      79              88             +9
     95 %      34              42             +8
     90 %      25              32             +7

i.e. at 95 % identity the family grows 34 -> 42, and separately 24 of the 34
already-known clusters were independently re-found from genomic DNA by the HMM.

A NOTE ON "100 % IDENTITY"
--------------------------
cd-hit computes identity over the ALIGNED region, so a short protein fully contained
in a longer one is merged even at ``-c 1.00``. That is why the 100 % union is 88 and
not 79 + 55 = 134: "100 % identity" here does NOT mean "identical string". This is
called out in the `note` column of every row so a reader cannot misread it.

OUTPUTS (written into the discovery dir)
----------------------------------------
  family_census.csv          one row per identity threshold — the headline table
  family_census_members.csv  one row per union member, with its cluster at each
                             threshold; also the side map from the synthetic
                             S####/H#### ids back to the real (pipe- and
                             space-bearing) labels
  seed_qc/seed_status.csv    one row per DISTINCT seed protein: was it re-found, at
                             what identity, by which homolog — the join point the
                             missed-seed feature consumes

USAGE
-----
  python3 scripts/family_census.py --discovery-dir DIR --seeds FASTA [--cpu N]

Every writer here is defensive: cd-hit missing, seeds unreadable, no hit table — each
degrades to a smaller but still-written census and logs a one-line reason. Nothing
in this module raises.

TWO RULES THAT MATTER MORE THAN THE NUMBERS
-------------------------------------------
1. **{} and no files, together.** The three CSVs are built in memory, staged as
   ``.tmp``, and os.replace()d into place only once ALL of them exist. Any failure —
   including one that happens after an earlier census succeeded — removes the .tmp
   files AND any stale ``family_census*.csv`` / ``seed_qc/seed_status.csv``. A census
   that returns {} leaves nothing behind that a reader could mistake for an answer.
2. **A degraded census says so in the sentence it is quoted by.** Without cd-hit the
   union is grouped by exact string, which on the gp75 reference reports 127 clusters
   (+48 new) instead of 88 (+9) and 7 re-found seeds instead of 46. Every rendering of
   that headline is prefixed ``APPROXIMATE — …``, and the summary carries
   ``census_complete: False`` plus ``clustering_method`` so a report can refuse to
   quote it without parsing prose.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
from pathlib import Path

import pandas as pd

# scripts/ and engine/ on sys.path so this works both as a CLI and as an import from
# hmm_finder / generate_report (which live in scripts/ themselves).
_HERE = Path(__file__).resolve().parent
for _p in (_HERE, _HERE.parent / "engine"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import export_csv                                          # noqa: E402
from run_selection import best_run_index                   # noqa: E402
# ONE cd-hit wrapper and ONE .clstr parser, shared with the clinker step, so the two
# clustering call sites cannot drift apart. Importing it also runs env_paths, which
# puts the conda env's cd-hit on PATH.
from cluster_and_clinker_corrected import cdhit, parse_clstr  # noqa: E402

try:
    from pipeline.synteny import nt_accession_in           # noqa: E402
except Exception:                                          # engine/ not importable
    def nt_accession_in(text: str) -> str:                 # type: ignore[misc]
        return ""


# ---------------------------------------------------------------------------
# Pinned cd-hit configuration
# ---------------------------------------------------------------------------
# Thresholds are ordered loosest-last so the CSV reads 100 % -> 90 %.
THRESHOLDS = (1.00, 0.95, 0.90)
HEADLINE_THRESHOLD = 0.95   # the row quoted in the report: strain-level, not species

# Word size: cd-hit REQUIRES -n 5 for -c >= 0.7 (it errors out otherwise), so this is
# fixed by the thresholds above rather than being a free choice.
CENSUS_WORD_SIZE = 5
# -aL 0 = cd-hit's default: no minimum alignment coverage of the longer sequence.
# Pinned EXPLICITLY rather than left implicit because it is the flag that makes a
# contained sequence merge at -c 1.00 (see the module docstring), which is the single
# most misreadable thing in this table. Changing it changes 88.
CENSUS_ALIGN_LONGER = 0.0
# Echoed verbatim into every CSV row so the numbers can be reproduced from the file
# alone. -T (threads) is deliberately excluded: it changes runtime, never the result,
# and including it would make two machines' censuses look different when they are not.
CDHIT_FLAGS = (f"-n {CENSUS_WORD_SIZE} -aL {CENSUS_ALIGN_LONGER:g} -M 0 -d 0")

SEED_PREFIX = "S"
HIT_PREFIX = "H"

# Prepended to EVERY rendering of the headline (summary dict, log line, CSV column) when
# the row it describes came from the exact-sequence fallback rather than cd-hit. Kept as
# one constant so the three renderings cannot drift apart, and worded so that the two
# things it gets wrong are named rather than hinted at.
APPROX_PREFIX = ("APPROXIMATE — cd-hit unavailable; exact-string grouping OVER-COUNTS "
                 "the union and UNDER-COUNTS re-found seeds: ")

# Everything this module writes. Listed once because it is both what a successful census
# commits and what a failed one has to remove.
_OUTPUT_NAMES = ("family_census.csv", "family_census_members.csv",
                 "seed_qc/seed_status.csv")
_OUTPUT_GLOBS = ("family_census*.csv", "family_census*.csv.tmp")

_METAGENOME_RE = re.compile(r"^(?:MAG|TPA)\b|uncultured|metagenom|environmental", re.I)


# ---------------------------------------------------------------------------
# Pure units
# ---------------------------------------------------------------------------
def read_fasta(path) -> list[tuple[str, str]]:
    """Read a FASTA into ``[(header_without_'>', sequence), ...]`` in file order.

    Deliberately hand-rolled rather than via Bio.SeqIO: seed headers are user-supplied
    and routinely carry pipes, colons and spaces, and we must keep the header EXACTLY
    as written so the census can be joined back to the user's own file. Returns [] on
    any read failure.
    """
    try:
        text = Path(path).read_text(errors="replace")
    except OSError:
        return []
    recs: list[tuple[str, str]] = []
    header, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                recs.append((header, "".join(buf)))
            header, buf = line[1:].strip(), []
        elif header is not None:
            buf.append(line.strip())
    if header is not None:
        recs.append((header, "".join(buf)))
    return [(h, s) for h, s in recs if s]


def normalize_aa(seq: str) -> str:
    """Comparable amino-acid string: upper-case, no whitespace, no trailing stop.

    Without this, `MTD…RT` and `MTD…RT*` are two 'distinct' seed proteins and the
    duplicate count silently drops.
    """
    return re.sub(r"[\s*]+", "", str(seq or "")).upper()


def dedup_by_sequence(records) -> tuple[list[tuple[str, str]], dict[str, str]]:
    """Collapse identical sequences -> ``(reps, header_to_rep)``.

    `reps` keeps the FIRST header seen for each distinct sequence, in file order;
    `header_to_rep` maps every input header (including the representative's own) to
    that representative's header.

    This must run BEFORE clustering. The gp75 seed file has 101 records but only 79
    distinct proteins; clustering the raw file makes the denominator 101, so "how many
    of my seeds came back?" is understated by 22 phantom seeds.
    """
    reps: list[tuple[str, str]] = []
    header_to_rep: dict[str, str] = {}
    first_by_seq: dict[str, str] = {}
    for header, seq in records:
        key = normalize_aa(seq)
        if not key:
            continue
        if key not in first_by_seq:
            first_by_seq[key] = header
            reps.append((header, key))
        header_to_rep[header] = first_by_seq[key]
    return reps, header_to_rep


def classify_clusters(clusters, seed_ids, hit_ids) -> dict:
    """Partition clusters into shared / seed-only / new.

    `clusters` is ``{cluster_id: [(member_id, is_representative), ...]}``.

      shared     — contains at least one seed AND at least one discovered homolog:
                   the HMM re-found something already in the input
      seed_only  — seeds only: an input protein nothing in the searched databases
                   matched at this identity
      new        — homologs only: a protein the family did not previously contain

    `shared + seed_only` is the family BEFORE discovery at this identity, and
    `shared + seed_only + new` is the family AFTER — that subtraction is the whole
    point of the file. Also returns `seeds_refound`, the seed member ids that landed
    in a shared cluster, so the per-threshold percentage and seed_status.csv cannot
    disagree about what "re-found" means.
    """
    seed_ids, hit_ids = set(seed_ids), set(hit_ids)
    out = {"n_union": 0, "n_shared": 0, "n_seed_only": 0, "n_new": 0,
           "n_unclassified": 0, "seeds_refound": set(), "member_to_cluster": {},
           "cluster_class": {}, "reps": set()}
    for cid, members in (clusters or {}).items():
        names = [m for m, _ in members]
        seeds = [m for m in names if m in seed_ids]
        hits = [m for m in names if m in hit_ids]
        if seeds and hits:
            klass = "shared"
            out["seeds_refound"].update(seeds)
        elif seeds:
            klass = "seed_only"
        elif hits:
            klass = "new"
        else:
            # A member neither list knows about: impossible unless a caller passed a
            # mismatched FASTA. Counted separately so the partition check below fails
            # loudly instead of the class counts quietly not summing to the union.
            klass = "unclassified"
        out[f"n_{klass}"] += 1
        out["n_union"] += 1
        out["cluster_class"][cid] = klass
        for name, is_rep in members:
            out["member_to_cluster"][name] = cid
            if is_rep:
                out["reps"].add(name)
    return out


def exact_clusters(members) -> dict:
    """cd-hit-free fallback: group union members by IDENTICAL sequence.

    `members` is ``[(member_id, sequence), ...]``. Same return shape as `cdhit`, so
    `classify_clusters` does not care which produced it. This is a STRICTLY weaker
    100 %-identity grouping than cd-hit's — it cannot merge a shorter protein into a
    longer one that contains it — so it over-counts the union relative to a cd-hit
    run. The census says so in the row's `note`.
    """
    by_seq: dict[str, list[str]] = {}
    for mid, seq in members:
        by_seq.setdefault(normalize_aa(seq), []).append(mid)
    return {i: [(m, j == 0) for j, m in enumerate(g)]
            for i, g in enumerate(by_seq.values())}


def cluster_at(faa, ident: float, workdir, cpu: int = 1) -> dict:
    """Cluster `faa` at `ident` -> ``{cluster_id: [(member_id, is_rep), ...]}``.

    Returns {} — never raises — when cd-hit is absent or fails, which is the whole
    reason the caller can still emit a (reduced) census on a machine without it.
    """
    faa, workdir = Path(faa), Path(workdir)
    if not faa.exists() or faa.stat().st_size == 0:
        return {}
    try:
        workdir.mkdir(parents=True, exist_ok=True)
        prefix = workdir / f"union_c{int(round(ident * 100)):03d}"
        return cdhit(faa, prefix, ident=ident, word=CENSUS_WORD_SIZE,
                     aL=CENSUS_ALIGN_LONGER, cpu=max(1, int(cpu or 1)))
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# Where the 55 homolog proteins come from
# ---------------------------------------------------------------------------
def homolog_table(discovery) -> pd.DataFrame:
    """The distinct homolog PROTEINS of a discovery run (55 on the gp75 reference).

    Rebuilt through `export_csv._dedup_hits` rather than read off a CSV so this file
    can never disagree with `hits_deduplicated.csv`, the paper table or the tree about
    what a homolog is. The canonical run is chosen by `run_selection.best_run_index`
    over the per-run frames of `all_runs_hits.csv`, exactly as `export_csv.export`
    does, and every round is passed as `all_rounds` so a locus keeps ONE identity
    across iterations.

    Falls back to a previously exported `hits_deduplicated.csv` if the raw hit table
    is missing, and to an empty frame if both are.
    """
    discovery = Path(discovery)
    allh_path = discovery / "all_runs_hits.csv"
    try:
        allh = pd.read_csv(allh_path, dtype=str).fillna("")
    except Exception:
        allh = pd.DataFrame()
    if not allh.empty and "run_label" in allh.columns:
        frames: dict[int, pd.DataFrame] = {}
        for label, g in allh.groupby("run_label", sort=True):
            try:
                frames[int(float(label))] = g.reset_index(drop=True)
            except (TypeError, ValueError):
                # Non-numeric run labels: fall back to first-seen order so the
                # selection still runs instead of dropping the round entirely.
                frames[len(frames) + 1] = g.reset_index(drop=True)
        if frames:
            best = best_run_index({i: f.to_dict("records") for i, f in frames.items()})
            dedup = export_csv._dedup_hits(frames.get(best, next(iter(frames.values()))),
                                           all_rounds=allh)
            if dedup is not None and not dedup.empty:
                return dedup
    try:
        return pd.read_csv(discovery / "hits_deduplicated.csv", dtype=str).fillna("")
    except Exception:
        return pd.DataFrame()


# ---------------------------------------------------------------------------
# Seed metadata
# ---------------------------------------------------------------------------
def seed_organism(header: str, accession: str = "") -> str:
    """Readable organism from an underscore-joined seed header.

    ``Escherichia_phage_vB_Eco_SPSP_OV049961.1`` -> ``Escherichia phage vB Eco SPSP``.
    Underscores become spaces because that is what `canonical.canonical_organism` (and
    a human) needs to see the ``phage <name>`` structure; the exact original header is
    preserved verbatim in the `label` column of family_census_members.csv, so nothing
    is lost by normalising here.
    """
    s = str(header or "").strip()
    if accession:
        s = s.replace(accession, " ")
    return re.sub(r"\s+", " ", s.replace("_", " ")).strip(" _.-")


def accession_class(label: str, accession: str) -> str:
    """`refseq` | `genbank` | `metagenome` | `unresolved` for a seed's accession.

    Precedence is deliberate and worth stating, because two of these describe
    different things: `metagenome` is about how the SAMPLE was obtained (a MAG or a
    TPA assembly, which has no isolate and often no host), while refseq/genbank is
    merely which database issued the id. Sample provenance is the more useful fact
    when a seed is not re-found — an assembled-from-reads contig is a much more likely
    false seed than a cultured isolate — so it wins.
    """
    if not accession:
        return "unresolved"
    if _METAGENOME_RE.search(str(label or "").replace("_", " ").strip()):
        return "metagenome"
    # RefSeq accessions are two letters, an UNDERSCORE, then digits (NC_, NZ_, NG_…);
    # primary GenBank/ENA/DDBJ accessions never carry that underscore.
    return "refseq" if "_" in accession.split(".")[0] else "genbank"


def accession_prefix(accession: str) -> str:
    """Leading letters (plus the RefSeq underscore) of an accession: ``NC_``, ``OV``."""
    m = re.match(r"^([A-Za-z]+_?)", str(accession or ""))
    return m.group(1) if m else ""


def _read_seed_recovery(discovery) -> dict:
    """{seed_id: row} from seed_qc/seed_recovery.csv ({} if it was never written).

    `seed_recovery.py` keys on the FIRST WHITESPACE TOKEN of each seed header, so the
    join key here must be built the same way or the whole column comes back empty.
    """
    path = Path(discovery) / "seed_qc" / "seed_recovery.csv"
    try:
        df = pd.read_csv(path, dtype=str).fillna("")
    except Exception:
        return {}
    if "seed_id" not in df.columns:
        return {}
    return {str(r["seed_id"]): dict(r) for _, r in df.iterrows()}


# ---------------------------------------------------------------------------
# The census
# ---------------------------------------------------------------------------
def _headline(t: float, c: dict, approximate: bool = False) -> str:
    """The one sentence a report quotes. `approximate` is NOT cosmetic.

    Without cd-hit the union is grouped by exact string, which cannot merge a short
    protein into the longer one containing it. On the gp75 reference that turns the
    true "88 (+9 new); 46 re-found" into "127 (+48 new); 7 re-found" — a 5x
    overstatement of novelty, which is the paper's actual claim. So the flag is
    threaded in from the caller and the sentence says so in its first words; a
    truncated quote of it still carries the caveat.
    """
    known = c["n_shared"] + c["n_seed_only"]
    pct = f"{t * 100:.0f}%"
    text = (f"at {pct} identity the family grows {known} -> {c['n_union']} protein "
            f"clusters (+{c['n_new']} new); {c['n_shared']} of the {known} known "
            f"clusters were re-found from genomic DNA")
    return (APPROX_PREFIX + text) if approximate else text


def _note(t: float, c: dict, cdhit_here: bool, complete: bool) -> str:
    parts = []
    if not cdhit_here:
        parts.append("cd-hit unavailable: clustered by EXACT sequence, which cannot "
                     "merge a short protein into a longer one containing it, so the "
                     "union is an over-count")
    elif t >= 1.0:
        parts.append("cd-hit identity is measured over the ALIGNED region, so a "
                     "shorter protein fully contained in a longer one clusters with "
                     "it even at -c 1.00; that is why the union is smaller than "
                     "n_seed_proteins + n_homolog_proteins. '100% identity' here does "
                     "NOT mean 'identical string'")
    parts.append("family before discovery = n_clusters_shared + n_clusters_seed_only; "
                 "family after = n_clusters_union")
    if not complete:
        parts.append("rows for the remaining identity thresholds are ABSENT (cd-hit "
                     "missing or failed), so this file is not the full census")
    if c.get("n_unclassified"):
        parts.append(f"INTERNAL CHECK FAILED: {c['n_unclassified']} cluster(s) held no "
                     "recognised member id")
    if c["n_union"] != c["n_shared"] + c["n_seed_only"] + c["n_new"] + c["n_unclassified"]:
        parts.append("INTERNAL CHECK FAILED: cluster classes do not sum to the union")
    return "; ".join(parts)


def _unlink(path: Path) -> bool:
    """Best-effort delete; True if the file is gone because we removed it."""
    try:
        path.unlink()
        return True
    except OSError:
        return False


def _purge_outputs(discovery: Path) -> list[Path]:
    """Remove every census file (and half-written .tmp) under `discovery`.

    Called on EVERY path that returns {}. A previous run's family_census.csv sitting in
    the discovery dir after a census that failed is indistinguishable — to a reader, to
    stage_summary, to the packager — from a census that succeeded, which is exactly the
    "confident result from a step that never ran" failure this file exists to stop.
    Returns what was actually removed so the caller can say so out loud.
    """
    removed: list[Path] = []
    try:
        candidates: list[Path] = []
        for pattern in _OUTPUT_GLOBS:
            candidates.extend(sorted(discovery.glob(pattern)))
        for rel in _OUTPUT_NAMES:
            candidates.append(discovery / rel)
            candidates.append(discovery / (rel + ".tmp"))
        seen: set[Path] = set()
        for path in candidates:
            if path in seen:
                continue
            seen.add(path)
            if path.is_file() and _unlink(path):
                removed.append(path)
    except Exception:
        pass
    return removed


def _abort(discovery: Path, log, reason: str) -> dict:
    """Log the one-line reason, leave no census behind, return {}."""
    try:
        log(f"  (family census skipped: {reason})")
    except Exception:
        pass
    removed = _purge_outputs(discovery)
    if removed:
        try:
            log("  (family census: removed stale " +
                ", ".join(p.name for p in removed) +
                " so nothing can quote a previous run's numbers as this run's)")
        except Exception:
            pass
    return {}


def _commit(discovery: Path, frames: list) -> dict:
    """Write `[(name, path, frame), ...]` all-or-nothing; return ``{name: path}``.

    Every frame goes to a sibling ``.tmp`` first and only becomes the real file once
    ALL of them are on disk, via os.replace (atomic on both POSIX and Windows). Writing
    the three CSVs straight to their names, as this used to, meant an exception in the
    third writer left the first two — or, after an earlier successful census, last run's
    complete and confidently wrong trio — in place while census() returned {}.
    """
    staged: list[tuple[str, Path, Path]] = []
    written: dict[str, str] = {}
    try:
        for name, out, frame in frames:
            out.parent.mkdir(parents=True, exist_ok=True)
            tmp = out.with_name(out.name + ".tmp")
            frame.to_csv(tmp, index=False)
            staged.append((name, out, tmp))
        for name, out, tmp in staged:
            os.replace(tmp, out)
            written[name] = str(out)
        return written
    except Exception:
        # Undo whatever landed: unreplaced .tmp files, and — if a replace succeeded
        # before a later one failed — the real files too, along with any stale
        # predecessor. Better no census than a partial one that reads as complete.
        for _name, _out, tmp in staged:
            _unlink(tmp)
        _purge_outputs(discovery)
        raise


def census(discovery, seeds_faa, cpu: int = 1, log=print) -> dict:
    """Run the family census; write the three CSVs; return a summary dict.

    Returns {} and logs a one-line reason on any failure — this is a reporting extra,
    and it must never take a finished discovery run down with it. Returning {} also
    means the discovery dir is left with NO census files at all (see `_purge_outputs`):
    the empty dict and the absent CSVs must always agree.
    """
    # Outside the try, because every failure path below needs a real Path to purge.
    try:
        discovery = Path(discovery)
    except Exception as e:
        try:
            log(f"  (family census skipped: unusable discovery dir {discovery!r}: {e})")
        except Exception:
            pass
        return {}

    try:
        seeds_faa = Path(seeds_faa)

        # --- the two inputs, each counted its own correct way --------------------
        seed_records = read_fasta(seeds_faa)
        if not seed_records:
            return _abort(discovery, log, f"no readable sequences in {seeds_faa}")
        reps, header_to_rep = dedup_by_sequence(seed_records)
        headers_per_rep: dict[str, list[str]] = {}
        for header, rep in header_to_rep.items():
            headers_per_rep.setdefault(rep, []).append(header)

        hits = homolog_table(discovery)
        if hits is None or hits.empty:
            return _abort(discovery, log,
                          "no homolog table — is this a finished run?")
        n_hit_rows = len(hits)

        # --- union members, on SYNTHETIC ids ------------------------------------
        # cd-hit splits names on whitespace and truncates them; seed headers carry
        # spaces, colons and pipes, so clustering on real labels silently collides
        # different seeds onto one name. S####/H#### are unambiguous, and
        # family_census_members.csv is the side map back to the real label.
        # The H ids are the dedup table's OWN homolog_id, so H0007 here and H0007 in
        # hits_deduplicated.csv / the tree / the paper table are the same protein.
        members: list[dict] = []
        for i, (header, seq) in enumerate(reps, 1):
            acc = nt_accession_in(header)
            members.append({"member_id": f"{SEED_PREFIX}{i:04d}", "member_type": "seed",
                            "label": header, "accession": acc, "seq": seq})
        for i, row in enumerate(hits.to_dict("records"), 1):
            hid = str(row.get("homolog_id") or "").strip() or f"{HIT_PREFIX}{i:04d}"
            label = (str(row.get("representative_organism") or "").strip()
                     or str(row.get("representative_genome") or "").strip() or hid)
            gid = str(row.get("representative_genome") or "").strip()
            members.append({"member_id": hid, "member_type": "homolog",
                            "label": label, "accession": nt_accession_in(gid) or gid,
                            "seq": normalize_aa(row.get("aa_sequence", ""))})
        members = [m for m in members if m["seq"]]
        seed_ids = [m["member_id"] for m in members if m["member_type"] == "seed"]
        hit_ids = [m["member_id"] for m in members if m["member_type"] == "homolog"]

        # A census of zero homologs is not a census. The drop above is silent by
        # design (one homolog row with a blank aa_sequence should not sink the table),
        # but if it takes ALL of them the run has nothing to cluster against the seeds,
        # and the clustering below would still succeed: every cluster would be
        # seed_only, n_new would be 0, and the file would announce in a confident
        # sentence that the family "grew" 79 -> 79. Stop here instead, and say how many
        # rows were read so the blank column is findable.
        if not hit_ids:
            return _abort(discovery, log,
                          f"{n_hit_rows} homolog row(s) read from the hit table but "
                          f"none carried a usable aa_sequence — nothing to compare "
                          f"the {len(seed_ids)} seed protein(s) against")
        if not seed_ids:
            return _abort(discovery, log,
                          f"{len(seed_records)} seed record(s) read from {seeds_faa} "
                          f"but none carried a usable sequence — a census with no "
                          f"'before' cannot say whether the family grew")

        # --- cluster the union at each threshold --------------------------------
        # `used_cdhit` is tracked PER ROW, not once for the run: if cd-hit clusters
        # 100 % and then fails at 95 %, the 100 % row must still advertise the real
        # cd-hit flags it was produced with rather than being relabelled a fallback.
        by_threshold: dict[float, dict] = {}
        used_cdhit: dict[float, bool] = {}
        with tempfile.TemporaryDirectory(prefix="family_census_") as td:
            work = Path(td)
            faa = work / "union.faa"
            faa.write_text("".join(f">{m['member_id']}\n{m['seq']}\n" for m in members))
            for t in THRESHOLDS:
                clusters = cluster_at(faa, t, work, cpu=cpu)
                if clusters:
                    by_threshold[t] = classify_clusters(clusters, seed_ids, hit_ids)
                    used_cdhit[t] = True
                    continue
                if t == THRESHOLDS[0]:
                    # cd-hit missing: emit the 100 % row from exact grouping so the
                    # user still gets a union count, then stop. Never abort.
                    log("  (family census: cd-hit unavailable — 100% row only, "
                        "grouped by exact sequence)")
                    by_threshold[t] = classify_clusters(
                        exact_clusters([(m["member_id"], m["seq"]) for m in members]),
                        seed_ids, hit_ids)
                    used_cdhit[t] = False
                break   # one failed threshold means cd-hit is unusable; stop cleanly
        have_cdhit = len(by_threshold) == len(THRESHOLDS) and all(used_cdhit.values())

        if not by_threshold:
            return _abort(discovery, log, "clustering produced nothing")

        # Build ALL THREE frames before ANY of them is written, then swap them into
        # place together: a frame that fails to build must not leave the other two,
        # nor a previous census, on disk (see `_commit`).
        written = _commit(discovery, [
            _census_frame(discovery, by_threshold, used_cdhit,
                          len(seed_records), len(seed_ids), len(hit_ids)),
            _members_frame(discovery, members, by_threshold),
            _seed_status_frame(discovery, members, reps, headers_per_rep, by_threshold),
        ])

        head_t = HEADLINE_THRESHOLD if HEADLINE_THRESHOLD in by_threshold \
            else max(by_threshold)
        head_approx = not used_cdhit.get(head_t, False)
        n_cdhit_rows = sum(1 for t in by_threshold if used_cdhit.get(t))
        method = ("cd-hit" if n_cdhit_rows == len(by_threshold)
                  else "exact_sequence" if n_cdhit_rows == 0 else "mixed")
        summary = {
            "n_seed_headers": len(seed_records),
            "n_seed_proteins": len(seed_ids),
            "n_homolog_proteins": len(hit_ids),
            "cdhit_flags": CDHIT_FLAGS,
            "used_cdhit": have_cdhit,
            # A downstream report has to be able to refuse to quote an approximate
            # census WITHOUT string-matching the headline, so the two facts that make
            # it unquotable are keys of their own: was every threshold clustered, and
            # by what. census_complete is False whenever any row came from the
            # exact-sequence fallback or any threshold is missing.
            "census_complete": have_cdhit,
            "clustering_method": method,
            "headline_approximate": head_approx,
            "headline_threshold": head_t,
            "headline": _headline(head_t, by_threshold[head_t],
                                  approximate=head_approx),
            "thresholds": {f"{t:.2f}": {
                "n_clusters_union": c["n_union"],
                "n_clusters_shared": c["n_shared"],
                "n_clusters_seed_only": c["n_seed_only"],
                "n_clusters_new": c["n_new"],
                "n_clusters_known": c["n_shared"] + c["n_seed_only"],
                "n_seeds_refound": len(c["seeds_refound"]),
            } for t, c in by_threshold.items()},
            "files": written,
        }
        log(f"Family census: {len(seed_records)} seed records = {len(seed_ids)} distinct "
            f"proteins; {len(hit_ids)} homolog proteins. " + summary["headline"] + ".")
        return summary
    except Exception as e:
        return _abort(discovery, log, f"{type(e).__name__}: {e}")


# ---------------------------------------------------------------------------
# Frame builders (each returns ONE (name, path, frame) triple for `_commit`;
# each already inside census's guard). They deliberately do not touch the disk:
# nothing may be written until all three have been built successfully.
# ---------------------------------------------------------------------------
def _census_frame(discovery: Path, by_threshold: dict, used_cdhit: dict,
                  n_headers: int, n_seeds: int, n_hits: int) -> tuple:
    """family_census.csv — the headline table, one row per identity threshold."""
    rows = []
    complete = len(by_threshold) == len(THRESHOLDS)
    for t in THRESHOLDS:
        c = by_threshold.get(t)
        if not c:
            continue
        cdhit_here = bool(used_cdhit.get(t))
        # Per-SEED, not per-cluster: the same definition seed_status.csv uses for
        # refound_*, so the percentage here always equals that column's mean.
        pct = (100.0 * len(c["seeds_refound"]) / n_seeds) if n_seeds else 0.0
        rows.append({
            "identity_threshold": f"{t:.2f}",
            "n_seed_proteins": n_seeds,
            "n_seed_headers": n_headers,
            "n_homolog_proteins": n_hits,
            "n_clusters_union": c["n_union"],
            "n_clusters_shared": c["n_shared"],
            "n_clusters_seed_only": c["n_seed_only"],
            "n_clusters_new": c["n_new"],
            "pct_seeds_refound": round(pct, 1),
            # Prefixed per ROW, from that row's own cdhit_here: the `note` column
            # already carries the caveat, but the headline is the string a report
            # quotes, usually without the note beside it.
            "headline": _headline(t, c, approximate=not cdhit_here),
            "cdhit_flags": f"-c {t:.2f} {CDHIT_FLAGS}" if cdhit_here
                           else "(none: cd-hit unavailable, exact-sequence grouping)",
            "note": _note(t, c, cdhit_here, complete),
        })
    return ("family_census", discovery / "family_census.csv", pd.DataFrame(rows))


def _members_frame(discovery: Path, members: list, by_threshold: dict) -> tuple:
    """family_census_members.csv — the S####/H#### side map plus per-threshold clusters.

    The `_95` columns stay BLANK when 95 % was not clustered (cd-hit absent) rather
    than being quietly filled from another identity — a column named for a threshold
    must never hold a different one's answer.
    """
    head = by_threshold.get(HEADLINE_THRESHOLD, {})
    rows = []
    for m in members:
        row = {
            "member_id": m["member_id"],
            "member_type": m["member_type"],
            "label": m["label"],
            "accession": m["accession"],
            "aa_len": len(m["seq"]),
        }
        for t in THRESHOLDS:
            c = by_threshold.get(t) or {}
            cid = (c.get("member_to_cluster") or {}).get(m["member_id"])
            row[f"cluster_{int(round(t * 100))}"] = "" if cid is None else cid
        cid95 = (head.get("member_to_cluster") or {}).get(m["member_id"])
        row["is_cluster_rep_95"] = m["member_id"] in (head.get("reps") or set())
        row["cluster_class_95"] = (head.get("cluster_class") or {}).get(cid95, "")
        rows.append(row)
    return ("family_census_members", discovery / "family_census_members.csv",
            pd.DataFrame(rows))


def _seed_status_frame(discovery: Path, members: list, reps: list,
                       headers_per_rep: dict, by_threshold: dict) -> tuple:
    """seed_qc/seed_status.csv — one row per DISTINCT seed protein.

    Written next to seed_recovery.csv because seed_qc/ is already packaged wholesale,
    and because the two files answer the halves of one question: seed_recovery says
    whether the MODEL scores a seed, this says whether the SEARCH found it in a genome.
    A seed can be recovered by the HMM and still be missing from the census (its genome
    was never in a searched database) — that gap is what the missed-seed feature hunts,
    and it can only see it if both columns sit in one row.
    """
    recovery = _read_seed_recovery(discovery)
    by_id = {m["member_id"]: m for m in members}
    head = by_threshold.get(HEADLINE_THRESHOLD, {})
    # `status` must never read "missed" for everything just because cd-hit was absent,
    # so it falls back to the tightest threshold that WAS computed. The refound_* columns
    # stay blank for thresholds that were not, so the fallback is visible in the file.
    status_t = HEADLINE_THRESHOLD if HEADLINE_THRESHOLD in by_threshold \
        else (max(by_threshold) if by_threshold else None)

    rows = []
    # reps and the seed members are built in the same order from the same list, so zip
    # rather than look up by header — two seed records CAN share a header string.
    seed_members = [m for m in members if m["member_type"] == "seed"]
    for (header, _seq), m in zip(reps, seed_members):
        acc = m["accession"]
        dup_headers = [h for h in headers_per_rep.get(header, []) if h != header]
        # seed_recovery keys on the first whitespace token; identical sequences score
        # identically, so any of this protein's headers answers for all of them.
        rec = {}
        for h in [header] + dup_headers:
            rec = recovery.get(h.split()[0] if h.split() else h, {})
            if rec:
                break

        # Homologs sharing this seed's headline cluster, best (= lowest homolog_id,
        # which the dedup table sorts by breadth then bit score) first.
        cid95 = (head.get("member_to_cluster") or {}).get(m["member_id"])
        mates = sorted(mid for mid, c in (head.get("member_to_cluster") or {}).items()
                       if c == cid95 and by_id.get(mid, {}).get("member_type") == "homolog")

        refound = {t: (m["member_id"] in (by_threshold[t]["seeds_refound"]))
                   for t in THRESHOLDS if t in by_threshold}
        rows.append({
            "seed_id": header.split()[0] if header.split() else header,
            "duplicate_of": ";".join(dup_headers),   # headers that ARE duplicates of this row
            "n_seed_headers": len(headers_per_rep.get(header, [header])),
            "organism": seed_organism(header, acc),
            "accession": acc,
            "accession_prefix": accession_prefix(acc),
            "accession_class": accession_class(header, acc),
            "aa_len": len(m["seq"]),
            "recovered_by_hmm": str(rec.get("after_recovered", "")),
            "after_bit": rec.get("after_bit", ""),
            "refound_100": refound.get(1.00, ""),
            "refound_95": refound.get(0.95, ""),
            "refound_90": refound.get(0.90, ""),
            "best_homolog_id_95": mates[0] if mates else "",
            # `status` is the HEADLINE verdict, i.e. at 95 % identity — a seed matched
            # only at 90 % is a different protein by any strain-level standard.
            "status": "refound" if refound.get(status_t) else "missed",
            "member_id": m["member_id"],   # join key into family_census_members.csv
        })
    # seed_qc/ is created by `_commit` alongside the .tmp, not here: this builder is
    # not allowed to touch the disk.
    return ("seed_status", discovery / "seed_qc" / "seed_status.csv",
            pd.DataFrame(rows))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--discovery-dir", type=Path, required=True)
    ap.add_argument("--seeds", type=Path, required=True)
    ap.add_argument("--cpu", type=int, default=1)
    args = ap.parse_args()
    summary = census(args.discovery_dir, args.seeds, cpu=args.cpu)
    for name, path in (summary.get("files") or {}).items():
        print(f"  wrote {path}")
    if not summary:
        print("  (no census written)")


if __name__ == "__main__":
    main()
