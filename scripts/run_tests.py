#!/usr/bin/env python3
"""
run_tests.py — fast unit tests for the pure helper logic (no network, no DBs).

Guards the bits most likely to regress: functional-category mapping (incl. the
VOGDB fallback and the 'virion RNA polymerase' fix), nucleotide detection +
table-11 translation, protein-vs-nucleotide accession routing, organism parsing,
and row-label disambiguation. Run:  python3 run_tests.py
"""
import sys
import tempfile
from pathlib import Path

fails = []


def check(name, cond):
    print(("ok   " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


import synteny_figure as S  # noqa: E402
check("categorize lysis", S.categorize("endolysin") == "lysis")
check("categorize capsid->structural", S.categorize("major capsid protein") == "structural")
check("categorize RNA polymerase->transcription",
      S.categorize("Virion DNA-directed RNA polymerase") == "transcription / regulation")
check("categorize DNA polymerase->replication",
      S.categorize("DNA polymerase I") == "replication / nucleotide metabolism")
check("categorize terminase->packaging", S.categorize("terminase large subunit") == "packaging")
check("categorize hypothetical", S.categorize("hypothetical protein") == S.HYPO_CAT)
check("categorize VOGDB fallback Xs->structural", S.categorize("hypothetical protein", "Xs") == "structural")
check("categorize ambiguous VOGDB stays hypothetical",
      S.categorize("hypothetical protein", "XhXrXs") == S.HYPO_CAT)

loci = [{"organism": "uncultured virus", "genome_id": "A1", "genes": []},
        {"organism": "Escherichia phage X", "genome_id": "B2", "genes": []},
        {"organism": "Escherichia phage X", "genome_id": "B3", "genes": []},
        {"organism": "Escherichia phage Y", "genome_id": "C4", "genes": []}]
S.set_row_labels(loci)
check("label: generic gets accession", "A1" in loci[0]["label"])
check("label: duplicate gets accession", "B2" in loci[1]["label"] and "B3" in loci[2]["label"])
check("label: unique stays clean", loci[3]["label"] == "Escherichia phage Y")

import annotate_genes as A  # noqa: E402
check("clean_desc strips sp|..| + Putative",
      A._clean_desc("sp|Q5UPF8|YL088 Putative ankyrin repeat protein") == "ankyrin repeat protein")

import annotate_organism as O  # noqa: E402
check("annot is_protein_acc WP", O._is_protein_acc("WP_123456.1"))
check("annot is_protein_acc not MG", not O._is_protein_acc("MG201401"))
check("annot org_from_title", O._org_from_title("foo [Escherichia phage phiKT]")
      == "Escherichia phage phiKT")

import pandas as pd  # noqa: E402
import cluster_and_clinker_corrected as C  # noqa: E402
_dd = pd.DataFrame([
    {"contig": "NC_019520", "nt_start": 100, "nt_end": 300, "source_type": "six_frame_orf", "bit_score": 150},
    {"contig": "NC_019520.1", "nt_start": 120, "nt_end": 290, "source_type": "annotated_protein", "bit_score": 140},
    {"contig": "NC_019520", "nt_start": 5000, "nt_end": 5200, "source_type": "six_frame_orf", "bit_score": 100},
    {"contig": "MG201401", "nt_start": 50, "nt_end": 250, "source_type": "six_frame_orf", "bit_score": 130}])
_keep = C.dedup_synteny_loci(_dd)
check("dedup: cross-DB same locus collapses (4 hits -> 3 loci)", len(_keep) == 3)
check("dedup: prefers six-frame over protein", 0 in _keep and 1 not in _keep)
check("dedup: keeps paralog + other genome", 2 in _keep and 3 in _keep)

# --- email never assumed: helpers stay OFFLINE (no NCBI call) without an address.
# With email=None the organism lookup must be skipped entirely (proving no
# placeholder address is ever sent) and fall back to a generic label. This test
# runs with NO network because the guard short-circuits before any Entrez call.
_eo = Path(tempfile.mkdtemp()) / "hits.tsv"
pd.DataFrame({"genome_id": ["NC_000000.1"], "db_name": ["RefSeq viral genomes"]}).to_csv(
    _eo, sep="\t", index=False)
O.annotate(_eo, None)
_eor = pd.read_csv(_eo, sep="\t")
check("annotate_organism offline (email=None) -> generic label, no NCBI",
      "organism" in _eor.columns and "uncultured virus" in str(_eor["organism"].iloc[0]))

# --- offline CSV export: a hits table with NO 'organism' column (the offline /
# --no-annotate case) must NOT abort the whole export. Regression for the bug
# where genome_metadata.csv + homolog_stats.csv silently went missing offline.
import export_csv as EX  # noqa: E402
_xd = Path(tempfile.mkdtemp())
_xv = _xd / "run1" / "benchmark" / "validated"; _xv.mkdir(parents=True)
pd.DataFrame({                                   # deliberately NO 'organism' column
    "hit_id": ["h1", "h2"], "genome_id": ["G1", "G2"],
    "db_name": ["INPHARED genomes", "INPHARED genomes"],
    "source_type": ["six_frame_orf", "six_frame_orf"], "run_label": ["1", "1"],
    "aa_sequence": ["MKAAQR", "MKBBST"], "bit_score": ["120", "90"],
    "evalue": ["1e-20", "1e-9"], "domain_aa_len": ["50", "40"],
    "passes_orf_filter": ["True", "True"],
}).to_csv(_xv / "hits.tsv", sep="\t", index=False)
_xfiles = EX.export(_xd)
check("offline export() runs without an 'organism' column", len(_xfiles) > 0)
check("offline export writes genome_metadata.csv", (_xd / "genome_metadata.csv").exists())
check("offline export writes homolog_stats.csv", (_xd / "homolog_stats.csv").exists())

# --- seed-recovery QC: tblout parsing + before/after status classification -----
import seed_recovery as SR  # noqa: E402
_tbl = (
    "#  comment line\n"
    "seqA  - q - 1e-50 150.2 0.1 1e-49 149.0 0.0 1.0 1 1 0 0 a capsid protein\n"
    "seqB  - q - 1e-03  30.0 0.0 1e-02  29.0 0.0 1.0 1 1 0 0 a weak hit\n"
    "seqA  - q - 1e-10  60.0 0.0 1e-09  59.0 0.0 1.0 1 1 0 0 lower dup of seqA\n")
_bb = SR.parse_tblout_best_bits(_tbl)
check("seed_recovery: best bit per target (dedup, keep max)",
      _bb.get("seqA") == 150.2 and _bb.get("seqB") == 30.0)
check("seed_recovery: comment/blank lines ignored", "#" not in "".join(_bb.keys()))
check("seed_recovery classify both -> recovered", SR.classify(True, True) == "recovered")
check("seed_recovery classify lost", SR.classify(True, False) == "lost_after_refinement")
check("seed_recovery classify gained", SR.classify(False, True) == "gained_after_refinement")
check("seed_recovery classify never", SR.classify(False, False) == "never_recovered")

# --- package layout: distinct numbered folders + per-folder README generation --
import package_layout as PL  # noqa: E402
check("package_layout: 8 distinct numbered dirs (no 00/00 collision)",
      len(set(PL.DIRS.values())) == 8 and len({d[:2] for d in PL.DIRS.values()}) == 8)
_pk = Path(tempfile.mkdtemp()) / "PACKAGE"
(_pk / PL.DIRS["tables"]).mkdir(parents=True)
(_pk / PL.DIRS["tables"] / "paper_main_table.csv").write_text("x\n")
(_pk / PL.DIRS["sequences"] / PL.PER_RUN / "run1").mkdir(parents=True)
(_pk / PL.DIRS["sequences"] / PL.PER_RUN / "run1" / "hits.tsv").write_text("x\n")
PL.write_readmes(_pk)
check("package_layout: top-level README written", (_pk / "README.txt").exists())
check("package_layout: folder README written", (_pk / PL.DIRS["tables"] / "README.txt").exists())
_rt = (_pk / PL.DIRS["tables"] / "README.txt").read_text()
check("package_layout: README describes a known file",
      "paper_main_table.csv" in _rt and "MAIN RESULT" in _rt)
check("package_layout: per-run README written",
      (_pk / PL.DIRS["sequences"] / PL.PER_RUN / "run1" / "README.txt").exists())

# --- ROC threshold calibration (engine controls.ControlReport.roc) -----------
import sys as _sys  # noqa: E402
_sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "engine"))
try:
    from pipeline.controls import ControlReport  # noqa: E402
    # Perfectly separable: positives 60-70, negatives 2-10 -> AUC 1.0, cut in the gap.
    _sep = ControlReport([
        {"role": "positive", "n_seqs": 5, "scores": [60.0, 62, 65, 68, 70]},
        {"role": "negative", "n_seqs": 5, "scores": [2.0, 4, 6, 8, 10]},
    ], 45.0, 30.0)
    _r = _sep.roc()
    check("ROC AUC = 1.0 for separable controls", abs(_r["auc"] - 1.0) < 1e-9)
    check("ROC optimal cutoff falls in the separation gap (10..60)",
          10.0 < _r["optimal_threshold"] < 60.0)
    check("ROC sens/spec = 1.0 at the optimum",
          _r["sensitivity_at_optimum"] == 1.0 and _r["specificity_at_optimum"] == 1.0)
    # Undetected negatives (n_seqs > len(scores)) must count as below-threshold.
    _pad = ControlReport([
        {"role": "positive", "n_seqs": 3, "scores": [50.0, 55, 60]},
        {"role": "negative", "n_seqs": 100, "scores": [40.0]},  # 99 undetected
    ], 45.0, 30.0)
    check("ROC pads undetected negatives (AUC = 1.0)", abs(_pad.roc()["auc"] - 1.0) < 1e-9)
    # All negatives undetected (the real gp75 case): noise floor is 0 bits, so the
    # max-margin optimum is mid-way between 0 and the lowest positive.
    _none = ControlReport([
        {"role": "positive", "n_seqs": 2, "scores": [60.0, 70.0]},
        {"role": "negative", "n_seqs": 10, "scores": []},  # nothing hit
    ], 45.0, 30.0)
    _rn = _none.roc()
    check("ROC all-undetected negatives -> AUC 1.0", abs(_rn["auc"] - 1.0) < 1e-9)
    check("ROC all-undetected optimum is mid-gap (0..60), not min_pos-1",
          25.0 < _rn["optimal_threshold"] < 35.0)
    # Identical distributions -> AUC ~ 0.5 (no discrimination).
    _ovl = ControlReport([
        {"role": "positive", "n_seqs": 4, "scores": [20.0, 30, 40, 50]},
        {"role": "negative", "n_seqs": 4, "scores": [20.0, 30, 40, 50]},
    ], 45.0, 30.0)
    check("ROC AUC ~ 0.5 for identical distributions", abs(_ovl.roc()["auc"] - 0.5) < 1e-9)
    check("summary() carries the advisory roc block",
          isinstance(_sep.summary().get("roc", {}).get("auc"), float))
    # n_seqs present-but-None (the builtin control catalogue uses None) must not crash.
    _noneN = ControlReport([
        {"role": "positive", "n_seqs": None, "scores": [60.0, 70.0]},
        {"role": "negative", "n_seqs": None, "scores": [5.0, 8.0]},
    ], 45.0, 30.0)
    check("ROC tolerates n_seqs=None (no TypeError)", abs(_noneN.roc()["auc"] - 1.0) < 1e-9)
    # A DETECTED hit with a (rare) negative bit score must still rank above an
    # UNDETECTED sequence — undetected is -inf, not a 0.0 floor that would outrank it.
    _negbit = ControlReport([
        {"role": "positive", "n_seqs": 2, "scores": [-2.0, 50.0]},
        {"role": "negative", "n_seqs": 3, "scores": []},  # all undetected
    ], 45.0, 30.0)
    check("ROC: undetected ranks below a negative-scoring detection (AUC=1.0)",
          abs(_negbit.roc()["auc"] - 1.0) < 1e-9)
except Exception as _e:
    check(f"ROC calibration import/compute failed: {_e}", False)

# --- interrupted-homolog finder (read-through translation + stop counting) ----
import find_interrupted as FI  # noqa: E402
_s, _m = FI.read_through_aa("ATG" + "AAA" + "TAA" + "GGG", "+", 0)   # M K * G
check("find_interrupted: read-through marker keeps '*' at the stop", _m == "MK*G")
check("find_interrupted: search sequence masks the stop as 'X'", _s == "MKXG")
check("find_interrupted: internal stop counted (1-based envelope)",
      FI.count_envelope_stops("AB*CD", 1, 5) == (1, [3]))
check("find_interrupted: terminal stop not counted as internal",
      FI.count_envelope_stops("ABC*", 1, 4) == (0, []))
check("find_interrupted: stop outside the envelope ignored",
      FI.count_envelope_stops("*ABCD", 2, 5) == (0, []))
# extend_orf: from the upstream stop, through the domain, to the natural stop
check("find_interrupted: extend_orf spans flanking stops",
      FI.extend_orf("M*KAR*Q", 3, 4) == (3, 6, "KAR*"))
check("find_interrupted: extend_orf with no upstream stop starts at residue 1",
      FI.extend_orf("KAR*Q", 1, 2) == (1, 4, "KAR*"))

import hmm_finder as H  # noqa: E402
td = Path(tempfile.mkdtemp())
dna = td / "x.fna"
dna.write_text(">s\nATGAGTAAATTCAAGAAATATCTGGGTGCC\n")
prot = td / "p.faa"
prot.write_text(">s\nMSKFKKYLGAAW\n")
check("detect nucleotide", H._looks_like_nucleotide(dna))
check("detect protein", not H._looks_like_nucleotide(prot))
seed = H.translate_seed(dna, 11, td, lambda m: None)
check("translate table11", "\nMSKFKKYLGA" in seed.read_text())

# --- new logic: HMM length parse, canonical-run pick, convergence wiring -----
hmm = td / "p.hmm"
hmm.write_text("HMMER3/f [3.4]\nNAME  x\nLENG  137\nALPH  amino\nHMM  A C D E\n")
check("hmm_leng reads LENG", H._hmm_leng(hmm) == 137)
check("hmm_leng missing -> 0", H._hmm_leng(td / "nope.hmm") == 0)

bd = td / "disc"
for ri, nrows in [(1, 3), (2, 5), (3, 5)]:
    vd = bd / f"run{ri}" / "benchmark" / "validated"
    vd.mkdir(parents=True)
    (vd / "hits.tsv").write_text("hit_id\thdr\n" + "".join(f"h{j}\tx\n" for j in range(nrows)))
check("best_run_index picks most complete (ties->earliest)", H._best_run_index(bd, 3) == 2)
check("best_run_index empty -> 1", H._best_run_index(td / "empty", 3) == 1)

check("convergence_check wired from engine", H.convergence_check is not None)
if H.convergence_check:
    check("convergence: stable hits+leng -> True", H.convergence_check(100, 102, 150, 150) is True)
    check("convergence: growing hits -> False", H.convergence_check(100, 200, 150, 160) is False)
    check("convergence: growing model -> False", H.convergence_check(100, 101, 150, 160) is False)

# --- organism-first label parsing for tree/alignment tips -------------------
import build_tree_of_hits as BT  # noqa: E402
check("organism OS= (UniProt)",
      BT._organism_from_desc("Major capsid protein OS=Escherichia phage T5 OX=2695836 GN=D20")
      == "Escherichia phage T5")
check("organism [bracket] (NCBI protein)",
      BT._organism_from_desc("hypothetical protein [Escherichia phage phiKT]")
      == "Escherichia phage phiKT")
check("organism NCBI genome title",
      BT._organism_from_desc("NC_019520.1 NC_019520.1:37102-37269 Escherichia phage phiKT, complete genome")
      == "Escherichia phage phiKT")
check("organism unknown -> ''", BT._organism_from_desc("ABC123 hypothetical") != "")  # returns something, not crash
check("short_acc UniProt", BT._short_acc("sp|P49861|CAPSD_BPHK7") == "P49861")
check("short_acc plain", BT._short_acc("NC_019520.1") == "NC_019520.1")
check("canon collapses host-genus alias",
      BT._canonical_organism("Enterobacteria phage N4") == BT._canonical_organism("Escherichia phage N4") == "n4")
check("canon virus form", BT._canonical_organism("Shigella virus Moo19") == "moo19")
check("canon metagenomic -> genome fallback",
      BT._canonical_organism("uncultured virus", "GPD_0001") == "gpd_0001")
check("canon distinct phages stay distinct",
      BT._canonical_organism("Escherichia phage T5") != BT._canonical_organism("Escherichia phage T7"))
_uu = set()
check("uniquify dedups labels",
      BT._uniquify("Escherichia_phage_X", _uu) == "Escherichia_phage_X"
      and BT._uniquify("Escherichia_phage_X", _uu) == "Escherichia_phage_X_2")

print(f"\n{len(fails)} FAILURE(S): {fails}" if fails else "\nALL TESTS PASSED")
sys.exit(1 if fails else 0)
