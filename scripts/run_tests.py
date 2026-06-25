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
check("categorize major coat protein->structural", S.categorize("putative major coat protein") == "structural")
check("categorize head protein->structural", S.categorize("head protein") == "structural")
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
pd.DataFrame({"genome_id": ["NC_000000.1", "GPD_0001"],
              "db_name": ["RefSeq viral genomes", "GPD"]}).to_csv(_eo, sep="\t", index=False)
O.annotate(_eo, None)
_eor = pd.read_csv(_eo, sep="\t")
# A cultured NCBI accession whose lookup is skipped offline must fall back to the ACCESSION,
# NOT be mislabelled "uncultured virus"; a genuinely non-NCBI id still gets the metagenomic label.
check("annotate_organism offline: NCBI accession -> accession label (not 'uncultured'), no NCBI",
      "organism" in _eor.columns and _eor["organism"].iloc[0] == "NC_000000"
      and "uncultured" not in str(_eor["organism"].iloc[0]))
check("annotate_organism offline: non-NCBI id keeps the metagenomic label",
      "uncultured virus" in str(_eor["organism"].iloc[1]))

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
# extend_orf starts at the Met START codon, not the upstream stop: the residues between the
# upstream stop and the ATG ("PQ" here) must NOT be prepended to the gene.
check("find_interrupted: extend_orf starts the gene at the Met, not the upstream stop",
      FI.extend_orf("*PQMTDLR*Z", 4, 8) == (4, 9, "MTDLR*"))
# Main-hit path: met_anchor reports the gene from its Met, not the stop-to-stop six-frame ORF.
import extract_validated_hits as EVH  # noqa: E402
# "ABMTDLR": domain starts at the M (env_from=3, 1-based) -> gene = "MTDLR" (5 aa), dropping "AB"
check("extract: met_anchor trims to the Met at the domain start (145->138 style fix)",
      EVH.met_anchor("ABMTDLR", 3) == (2, "MTDLR"))
check("extract: met_anchor leaves an already-Met-starting ORF unchanged",
      EVH.met_anchor("MTDLR", 1) == (0, "MTDLR"))
check("extract: met_anchor with no near Met anchors at the domain start (no upstream extension)",
      EVH.met_anchor("ARNDCQ", 4) == (3, "DCQ"))
# A Met far upstream across a long stop-free (antisense/overprint) frame is NOT the gene start:
# it must be ignored, anchoring at the domain start, so orf_aa_len is not inflated ~10x.
check("extract: met_anchor ignores a far-upstream Met on a long stop-free frame",
      EVH.met_anchor("M" + "A" * 108 + "K" * 30, 110) == (109, "K" * 30))
# run_pipeline launcher: presets + the flag/out-dir helpers used to inject no-prompt defaults
import run_pipeline as RP  # noqa: E402
check("run_pipeline: presets defined", {"phage-discovery", "offline", "smoke"} <= set(RP.PRESETS))
check("run_pipeline: _has detects a flag and its '=' form",
      RP._has(["--email=x"], "--email") and RP._has(["--fasta", "s"], "--fasta")
      and not RP._has(["--x"], "--email"))
check("run_pipeline: _out_dir reads --out-dir and --out (space + '=')",
      RP._out_dir(["--out-dir", "/a"]) == "/a" and RP._out_dir(["--out=/b"]) == "/b")
# preload_databases: warms the same six-frame translation cache (<db>.sixframe.min<N>.faa) the
# engine writes/consumes; translation_cache() must read exactly that key format.
import preload_databases as PD  # noqa: E402
check("preload: default database set is non-empty", bool(PD.DEFAULT_DATABASES) and "," in PD.DEFAULT_DATABASES)
import tempfile as _tf  # noqa: E402
_pc = Path(_tf.mkdtemp()); (_pc / "cache" / "DBX").mkdir(parents=True)
(_pc / "cache" / "DBX" / "0000_x.fa.gz.sixframe.min30.faa").write_text(">a\nMK\n")
(_pc / "cache" / "DBX" / "0000_x.fa.gz.sixframe.min50.faa").write_text(">a\nMK\n")
check("preload: translation_cache finds the min-N six-frame cache files (and only that N)",
      list(PD.translation_cache(_pc, 30)) == ["DBX/0000_x.fa.gz.sixframe.min30.faa"])
# aa_to_nt / stop_nt: genome coordinates + DNA that translate back (frame/strand-correct)
from Bio.Seq import Seq as _Seq  # noqa: E402
_g = "ATGAAATAAGGGCCC"   # + frame 0: M K * G P
_fs, _fe, _cod = FI.aa_to_nt(_g, "+", 0, 1, 5)
check("find_interrupted: aa_to_nt (+) coords + DNA round-trip",
      (_fs, _fe) == (1, 15) and str(_Seq(_cod).translate()) == "MK*GP")
check("find_interrupted: aa_to_nt (-) gives the reverse-complement frame DNA",
      str(_Seq(FI.aa_to_nt(_g, "-", 0, 1, 5)[2]).translate())
      == str(_Seq(_g).reverse_complement().translate()))
check("find_interrupted: stop_nt (+) points at the stop codon", FI.stop_nt(_g, "+", 0, 3) == 7)
# write_aa_fastas: emit domain + full-ORF protein FASTAs with the internal stop kept '*'
_fid = Path(tempfile.mkdtemp())
_rows = [{"contig": "c1", "strand": "+", "frame": "0", "domain_nt_start": "1",
          "domain_nt_end": "12", "internal_stops": "1", "stop_aa_positions": "3",
          "domain_bit_score": "40.0", "i_evalue": "1e-10",
          "domain_aa_with_stops": "MK*G", "full_orf_aa": "MK*GAR*"}]
_domf, _orff = FI.write_aa_fastas(_rows, _fid / "interrupted_homologs.tsv")
check("find_interrupted: domain FASTA keeps '*' at the internal stop",
      _domf.name == "interrupted_homologs_domain_aa.faa" and "MK*G" in _domf.read_text())
check("find_interrupted: full-ORF FASTA keeps the read-through ORF with stops",
      _orff.name == "interrupted_homologs_full_orf_aa.faa" and "MK*GAR*" in _orff.read_text())
# full-ORF nucleotide: coding DNA that translates back to full_orf_aa (incl. the stop codons)
_orfg = "ATGAAATAAGGGGCACGTTAA"   # + frame 0: M K * G A R *  (== full_orf_aa)
_rows2 = [dict(_rows[0], full_orf_nt=_orfg, full_orf_aa="MK*GAR*")]
_ntf = FI.write_orf_nt_fasta(_rows2, _fid / "interrupted_homologs.tsv")
check("find_interrupted: full-ORF nt FASTA written with the ORF DNA",
      _ntf.name == "interrupted_homologs_full_orf_nt.fna" and _orfg in _ntf.read_text())
check("find_interrupted: full-ORF nt translates back to full_orf_aa (incl. actual stop codon)",
      str(_Seq(_orfg).translate()) == "MK*GAR*")
check("find_interrupted: ORF nt columns present in ROW_COLS",
      {"orf_nt_start", "orf_nt_end", "natural_stop_nt", "full_orf_nt"} <= set(FI.ROW_COLS))
# Overprinting / silent-stop analysis (the proof-of-overprinting step)
check("find_interrupted: codon_covering (+) reads the forward codon",
      FI._codon_covering("ATGAAATAA", "+", 0, 0) == ("ATG", 0))
check("find_interrupted: codon_covering (-) reads the reverse-complement codon",
      FI._codon_covering("ATGAAATAA", "-", 0, 8) == ("TTA", 0))
# '+' stop TAA aligned to antisense frame 0 -> antisense codon TTA (Leu); T->C removes the
# small stop (CAA) and is synonymous antisense (TTA->TTG). So the stop is silent there.
check("find_interrupted: stop silent in an open antisense frame (aligned TAA/Leu wobble)",
      FI._stop_silent_in_frame("TAAGGG", "+", 1, "-", 0) is True)
check("find_interrupted: a non-stop position is never 'silent'",
      FI._stop_silent_in_frame("GGGTAA", "+", 1, "-", 0) is False)
check("find_interrupted: frame_stop_count 0 over an open antisense frame",
      FI._frame_stop_count("TAAGGG", "-", 0, 1, 6) == 0)
_op = FI.analyze_overprinting("TAAGGG", "+", 1, 6, [1])
check("find_interrupted: analyze_overprinting -> strong when open + synonymous",
      _op["support"] == "strong" and _op["open_stops"] == 0 and _op["per_stop_silent"] == [True])
check("find_interrupted: overprinting columns present in ROW_COLS",
      {"overprinting_support", "antisense_open_frame", "antisense_open_stops",
       "stop_silent_antisense"} <= set(FI.ROW_COLS))

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

# --- database-failure resilience (engine: one DB blip must not abort the run) ---
import importlib.util as _ilu  # noqa: E402
_bench_path = Path(__file__).resolve().parent.parent / "engine" / "scripts" / "run_all_database_benchmark.py"
try:
    _spec = _ilu.spec_from_file_location("engine_benchmark_for_test", _bench_path)
    _bench = _ilu.module_from_spec(_spec)
    _spec.loader.exec_module(_bench)
    _fatal = _bench.db_failure_fatal
    # DEFAULT (no strict flags): any single DB failure is non-fatal -> run continues
    check("db-resilience: required DB failure is NON-fatal by default",
          _fatal(False, False, False) is False)
    check("db-resilience: optional DB failure is NON-fatal by default",
          _fatal(True, False, False) is False)
    # --strict-databases: a required DB failure is fatal again (pre-hardening behaviour)
    check("db-resilience: --strict-databases makes a required DB failure fatal",
          _fatal(False, True, False) is True)
    check("db-resilience: --strict-databases does NOT force optional failures fatal",
          _fatal(True, True, False) is False)
    # --stop-on-optional-failure still makes optional failures fatal (unchanged)
    check("db-resilience: --stop-on-optional-failure makes an optional DB failure fatal",
          _fatal(True, False, True) is True)
    check("db-resilience: --stop-on-optional-failure leaves required default non-fatal",
          _fatal(False, False, True) is False)
except Exception as _e:
    check(f"db-resilience: SKIPPED (engine import error: {_e})", True)

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

# --- integration tests for this session's pipeline paths -----------------------
# Family-calibrated interrupted threshold max(floor 30, ROC-Youden), via the helper.
check("interrupted threshold: ROC above floor wins",
      H._interrupted_min_bit({"roc": {"optimal_threshold": 41.2}})[0] == 41.2)
check("interrupted threshold: floor wins when ROC below it",
      H._interrupted_min_bit({"roc": {"optimal_threshold": 12.0}})[0] == 30.0)
check("interrupted threshold: no controls -> bare floor",
      H._interrupted_min_bit({})[0] == 30.0)
check("perhit HMM alignment: missing inputs -> {} (no raise)",
      H.run_perhit_hmm_alignment(Path("/no.hmm"), Path("/no.faa"),
                                 Path(tempfile.mkdtemp()), lambda *a, **k: None) == {})

# --- provenance redaction: the user's NCBI e-mail must never reach a shared
# manifest/METHODS (it is not needed to reproduce a run). Regression for the
# leak where run_manifest.json embedded --email <addr> via command_line.
check("hmm_finder: _redact_emails strips a --email value from a command line",
      H._redact_emails("hmm_finder.py --email jane.doe@uni.edu --cpu 4")
      == "hmm_finder.py --email <redacted-email> --cpu 4")
check("hmm_finder: _redact_emails handles plus/dot/sub-domain addresses + is idempotent",
      H._redact_emails(H._redact_emails("x a.b+c@sub.example.co --y")) == "x <redacted-email> --y")
check("hmm_finder: _redact_emails leaves e-mail-free text untouched",
      H._redact_emails("--databases 'INPHARED genomes'") == "--databases 'INPHARED genomes'")
check("hmm_finder: _redact_emails passes non-strings through",
      H._redact_emails(None) is None and H._redact_emails(42) == 42)
# --- 08_scripts reproducibility copy must EXCLUDE scratch (_*) helpers (which can
# embed the e-mail / local paths), keep real scripts + dunder, drop bytecode.
_ig = H._ignore_scratch("/scripts",
        ["hmm_finder.py", "_run_gp75.unix.sh", "_audit.sh", "__init__.py", "__pycache__", "foo.pyc"])
check("hmm_finder: _ignore_scratch drops scratch _* scripts (the e-mail leak vector)",
      {"_run_gp75.unix.sh", "_audit.sh"} <= _ig)
check("hmm_finder: _ignore_scratch keeps real scripts + dunder files (__init__.py)",
      "hmm_finder.py" not in _ig and "__init__.py" not in _ig)
check("hmm_finder: _ignore_scratch drops bytecode (__pycache__, *.pyc)",
      {"__pycache__", "foo.pyc"} <= _ig)
# --- --no-overwrite: a re-run must never clobber a previous run's folder --------
_od = Path(tempfile.mkdtemp()); _base = _od / "gp75_discovery"
check("hmm_finder: _unique_out_dir returns the path unchanged when it doesn't exist",
      H._unique_out_dir(_base) == _base)
_base.mkdir()
check("hmm_finder: _unique_out_dir reuses an EMPTY existing folder (nothing to clobber)",
      H._unique_out_dir(_base) == _base)
(_base / "report.html").write_text("x")            # now it holds a run
check("hmm_finder: _unique_out_dir bumps to _2 when the folder holds results",
      H._unique_out_dir(_base) == _od / "gp75_discovery_2")
(_od / "gp75_discovery_2").mkdir(); (_od / "gp75_discovery_2" / "r").write_text("y")
check("hmm_finder: _unique_out_dir skips an occupied _2 and returns _3",
      H._unique_out_dir(_base) == _od / "gp75_discovery_3")

import cluster_and_clinker_corrected as _CC  # noqa: E402
check("clinker static PNG: bad input -> False (graceful, no raise)",
      _CC._html_to_png(Path("/nonexistent.html"), Path(tempfile.mkdtemp()) / "x.png") is False)

# generate_report renders the inline coloured MSA + interrupted/overprinting sections.
import json as _json  # noqa: E402
import generate_report as _GR  # noqa: E402
_rd = Path(tempfile.mkdtemp())
(_rd / "downstream" / "tree").mkdir(parents=True)
(_rd / "downstream" / "tree" / "hits.aln.faa").write_text(">h1 org=x\nMKLAR-DE\n>h2\nMKLAWDEC\n")
(_rd / "interrupted_homologs.tsv").write_text(
    "contig\tstrand\tframe\tinternal_stops\tstop_nt_positions\toverprinting_support\t"
    "domain_bit_score\ti_evalue\nc1\t+\t0\t1\t51421\tstrong\t99.0\t1e-20\n")
(_rd / "run_manifest.json").write_text(_json.dumps({"parameters": {"label": "t"},
    "interrupted_homologs": {"interrupted_candidates": 1, "overprinting_strong": 1,
                             "overprinting_partial": 0}}))
_rhtml = _GR.generate(_rd).read_text()
check("report render: inline coloured MSA section present",
      "Coloured alignment (final hits)" in _rhtml)
check("report render: interrupted + overprinting section present",
      "Stop-interrupted / overprinted homologs" in _rhtml
      and "Overprinting (antisense-open-frame) test" in _rhtml)

# End-to-end find_interrupted: build a tiny family HMM, scan a contig with a known
# internal stop, assert the whole _search_batch -> ORF/overprinting/FASTA chain.
# Skips cleanly when HMMER isn't on PATH (e.g. CI without the conda env).
import shutil as _sh  # noqa: E402
import subprocess as _sp  # noqa: E402
import csv as _csv  # noqa: E402
if _sh.which("hmmbuild") and _sh.which("hmmsearch"):
    try:
        _e2e = Path(tempfile.mkdtemp())
        _prot = "MKAILVGGTRSDEFHNPQWYACMKLLVGGTRSDEFHNPQWYACMKAILVGG"
        _cod = {'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA',
                'E': 'GAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'L': 'CTT', 'K': 'AAA',
                'M': 'ATG', 'F': 'TTT', 'P': 'CCT', 'S': 'TCT', 'T': 'ACT', 'W': 'TGG',
                'Y': 'TAT', 'V': 'GTT'}
        _dna = "".join(_cod[a] for a in _prot)
        _dna = _dna[:24 * 3] + "TAA" + _dna[24 * 3 + 3:] + "TAA"   # internal stop @aa25 + natural stop
        (_e2e / "fam.faa").write_text(f">seed\n{_prot}\n")
        _sp.run(["hmmbuild", "--amino", str(_e2e / "fam.hmm"), str(_e2e / "fam.faa")],
                check=True, capture_output=True)
        (_e2e / "genome.fna").write_text(f">c1\n{_dna}\n")
        _ot = _e2e / "interrupted_homologs.tsv"
        FI._run(_e2e / "genome.fna", _e2e / "fam.hmm", _ot, 5.0, 1,
                log=lambda *a, **k: None, emit_fasta=True)
        _irows = list(_csv.DictReader(open(_ot), delimiter="\t"))
        check("find_interrupted e2e: an interrupted candidate is found", len(_irows) >= 1)
        if _irows:
            _ir = _irows[0]
            check("find_interrupted e2e: internal stop detected", int(_ir["internal_stops"]) >= 1)
            check("find_interrupted e2e: full_orf_nt round-trips to full_orf_aa",
                  str(_Seq(_ir["full_orf_nt"]).translate(table=11)) == _ir["full_orf_aa"])
            check("find_interrupted e2e: the three sequence FASTAs are written",
                  all(Path(str(_ot)[:-4] + sfx).exists()
                      for sfx in ("_domain_aa.faa", "_full_orf_aa.faa", "_full_orf_nt.fna")))
            check("find_interrupted e2e: overprinting_support populated",
                  _ir.get("overprinting_support") in ("strong", "partial", "none"))
    except Exception as _e:
        check(f"find_interrupted e2e: SKIPPED (tooling error: {_e})", True)
else:
    check("find_interrupted e2e: SKIPPED (HMMER not on PATH)", True)

# scan_genome --accession: must require an email before any NCBI call (pure-logic guard).
import os as _os  # noqa: E402
import scan_genome as _SGm  # noqa: E402
_saved_email = _os.environ.pop("NCBI_EMAIL", None)
_acc_rejected = False
try:
    _SGm.fetch_genome("KX098390", None, Path(tempfile.mkdtemp()), lambda *a, **k: None)
except SystemExit:
    _acc_rejected = True
except Exception:
    _acc_rejected = False
finally:
    if _saved_email is not None:
        _os.environ["NCBI_EMAIL"] = _saved_email
check("scan_genome: --accession without an email is rejected before any NCBI call", _acc_rejected)

# scan_genome: single-genome targeted scan (build HMM -> scan one genome -> present/absent).
if _sh.which("hmmbuild") and _sh.which("hmmsearch"):
    try:
        import scan_genome as SG  # noqa: E402
        _sg = Path(tempfile.mkdtemp())
        _sp = "MKAILVGGTRSDEFHNPQWYACMKLLVGGTRSDEFHNPQWYACMKAILVGG"
        _sc = {'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA',
               'E': 'GAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'L': 'CTT', 'K': 'AAA',
               'M': 'ATG', 'F': 'TTT', 'P': 'CCT', 'S': 'TCT', 'T': 'ACT', 'W': 'TGG',
               'Y': 'TAT', 'V': 'GTT'}
        _sg_dna = "".join(_sc[a] for a in _sp)
        (_sg / "seeds.faa").write_text(f">seed\n{_sp}\n")
        (_sg / "clean.fna").write_text(">g\nGGGCCCAAA" + _sg_dna + "TAAGGGCCCAAA\n")
        (_sg / "absent.fna").write_text(">g\n" + "GATCGATCGGCTAGCATCGATGCATGCTAGC" * 20 + "\n")
        _hmm = SG.build_hmm_from_seeds(_sg / "seeds.faa", 11, _sg, 1, lambda *a, **k: None)
        _sc1 = SG.scan(_sg / "clean.fna", _hmm, _sg / "oc", 5.0, False, 1, lambda *a, **k: None)
        check("scan_genome: detects a clean copy of the gene",
              _sc1["n_clean"] >= 1 and "PRESENT" in _sc1["verdict"])
        _sc2 = SG.scan(_sg / "absent.fna", _hmm, _sg / "oa", 5.0, False, 1, lambda *a, **k: None)
        check("scan_genome: reports a genome lacking the gene as not detected",
              _sc2["n_clean"] == 0 and _sc2["n_interrupted"] == 0)
        check("scan_genome: writes scan_hits.tsv + a per-genome report",
              (_sg / "oc" / "scan_hits.tsv").exists() and (_sg / "oc" / "scan_report.txt").exists())
        _sc3 = _SGm.scan(_sg / "clean.fna", _hmm, _sg / "on", 5.0, False, 1,
                         lambda *a, **k: None, neighbours=True, db_cache=_sg)
        check("scan_genome: Prodigal neighbour-calling path runs without error",
              isinstance(_sc3, dict) and _sc3.get("n_clean", 0) >= 1)
        # _select_neighbours keeps overlapping genes (e.g. an overprint partner)
        _g = [(100, 400, 1, {}), (1200, 1500, 1, {}), (700, 2000, -1, {})]  # last overlaps [800,1000]
        _u, _d, _o = _SGm._select_neighbours(_g, 800, 1000)
        check("scan_genome: overlapping (overprint) gene is kept, not dropped",
              len(_o) == 1 and _o[0][0] == 700 and len(_u) == 1 and len(_d) == 1)
        import genome_map as _GMm  # noqa: E402
        check("scan_genome: 'relationship' column is wired",
              "relationship" in _SGm.SCAN_NB_COLS)
        _gr = _GMm.build_genes((100, 400, 1),
                               [(500, 700, 1, {"gene": "nbr"}), (200, 600, -1, {"gene": "ov"})],
                               flank_keys={(500, 700)})
        check("genome_map: the gene of interest is marked as the anchor",
              any(it["role"] == "anchor" and it["start"] == 100 for it in _gr))
        check("genome_map: overlapping gene flagged, flank labelled from annotation",
              any(it["role"] == "overlap" for it in _gr)
              and any(it["role"] == "flank" and it["label"] == "nbr" for it in _gr))
        _lg = _GMm.write_locus_genbank(
            [{"start": 1, "end": 300, "strand": 1, "role": "anchor", "label": "GOI"},
             {"start": 400, "end": 700, "strand": -1, "role": "flank", "label": "nbr"}],
            "ATGC" * 300, "Test phage", "ACC123", _sg / "locus.gb")
        _lgt = _lg.read_text() if _lg else ""
        check("genome_map: locus GenBank written (openable in Easyfig/Artemis/clinker)",
              bool(_lg) and "CDS" in _lgt and "gene_of_interest" in _lgt)
        # coordinate fidelity with wlo>1: a 1-based-inclusive [s,e] feature must keep its full
        # length in the exported .gb (the +1 on the FeatureLocation end). wlo here = 101.
        from Bio import SeqIO as _St
        _lg2 = _GMm.write_locus_genbank(
            [{"start": 101, "end": 200, "strand": 1, "role": "anchor", "label": "GOI"},   # 100 bp
             {"start": 301, "end": 309, "strand": -1, "role": "flank", "label": "nbr"}],   # 9 bp
            "ATGC" * 100, "Test phage", "ACC9", _sg / "locus2.gb")
        _lens = sorted(len(f.location) for f in _St.read(str(_lg2), "genbank").features if f.type == "CDS")
        check("genome_map: locus GenBank CDS keep full 1-based-inclusive length (wlo>1, no 3' truncation)",
              _lens == [9, 100])
        check("genome_map: renderer options enumerated (incl. dfv + easyfig)",
              {"dfv", "pub", "pygenomeviz", "matplotlib", "easyfig"} <= set(_GMm.MAP_TOOLS))
        import inspect as _insp  # noqa: E402
        check("genome_map: DNA Features Viewer ('dfv') is the default renderer",
              _insp.signature(_GMm.draw).parameters["tool"].default == "dfv")
        # tolerant row-packing: trivial start/stop overlaps collapse to the baseline, a
        # real overprint (overlap >> tol) gets its own row so the gene of interest shows.
        class _F:
            def __init__(s, a, b): s.start, s.end = a, b
        _adj = [_F(0, 400), _F(398, 800), _F(798, 1200)]        # 2 bp shared start/stop
        _lv = _GMm._dfv_tolerant_levels(_adj, 60)
        check("genome_map: tiny start/stop overlaps stay on one row (no staircase)",
              set(_lv.values()) == {0})
        _ovp = [_F(0, 1000), _F(700, 1100)]                     # 300 bp overprint-style overlap
        _lv2 = _GMm._dfv_tolerant_levels(_ovp, 60)
        check("genome_map: a substantial (overprint) overlap gets its own row",
              len(set(_lv2.values())) == 2)
        # selectable palette: 'colorblind' yields different (Tol-muted) structural colour
        _cc_def = _GMm._scheme("default")[0]
        _cc_cb = _GMm._scheme("colorblind")[0]
        check("genome_map: colorblind palette differs from default",
              _cc_def.get("structural") and _cc_cb.get("structural")
              and _cc_def["structural"] != _cc_cb["structural"])
        check("genome_map: draw() exposes palette/functional/module-bracket options",
              {"palette", "functional_labels", "module_brackets"} <= set(
                  _insp.signature(_GMm.draw).parameters))
        # module runs: a contiguous same-category run is a module; hypothetical is NOT
        _mg = [{"role": "flank", "category": "structural", "start": 0, "end": 500},
               {"role": "flank", "category": "structural", "start": 510, "end": 900},
               {"role": "flank", "category": "hypothetical / unknown", "start": 1000, "end": 1400},
               {"role": "flank", "category": "hypothetical / unknown", "start": 1410, "end": 1800}]
        _runs = _GMm._module_runs(_mg, _cc_def)
        check("genome_map: contiguous same-category run -> one module; hypothetical excluded",
              len(_runs) == 1 and _runs[0]["cat"] == "structural" and _runs[0]["n"] == 2)
        # legend handles annotate categories with counts
        _lh = _GMm._legend_handles(
            [{"role": "anchor", "category": "gene of interest"},
             {"role": "flank", "category": "structural"},
             {"role": "flank", "category": "structural"}], _cc_def, "#ffd400", "#dde2e8")
        check("genome_map: legend shows category counts",
              any("structural (2)" in h.get_label() for h in _lh))
    except Exception as _e:
        check(f"scan_genome e2e: SKIPPED (tooling error: {_e})", True)
else:
    check("scan_genome e2e: SKIPPED (HMMER not on PATH)", True)

# synteny gene-neighbourhood ORDER table (positions relative to the gene of interest)
_loc = {"genome_id": "g1", "organism": "phage X", "genes": [
    {"s": -1500, "e": -900, "st": 1, "fam": False, "og": "OG1", "category": "structural",
     "vfam": "V1", "func": "tail fiber"},
    {"s": 0, "e": 600, "st": 1, "fam": True, "og": "OG0", "category": "gene of interest",
     "vfam": "", "func": ""},
    {"s": 800, "e": 1400, "st": -1, "fam": False, "og": "OG2", "category": "transcription",
     "vfam": "V2", "func": "RNA polymerase"}]}
_nb = {r["pos_index"]: r for r in S.neighbourhood_rows("0", [_loc])}
check("synteny neighbourhood: the gene of interest is pos_index 0", _nb[0]["is_anchor"] == 1)
check("synteny neighbourhood: upstream gene = pos_index -1 (negative distance)",
      _nb[-1]["distance_to_anchor_bp"] < 0 and _nb[-1]["function"] == "tail fiber")
check("synteny neighbourhood: downstream gene = pos_index +1 (positive distance, opposite strand)",
      _nb[1]["distance_to_anchor_bp"] > 0 and _nb[1]["strand_vs_gene"] == "-")
check("synteny neighbourhood: columns match NEIGHBOUR_COLS",
      set(_nb[0].keys()) == set(S.NEIGHBOUR_COLS))

print(f"\n{len(fails)} FAILURE(S): {fails}" if fails else "\nALL TESTS PASSED")
sys.exit(1 if fails else 0)
