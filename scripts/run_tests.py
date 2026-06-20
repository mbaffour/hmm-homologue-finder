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

print(f"\n{len(fails)} FAILURE(S): {fails}" if fails else "\nALL TESTS PASSED")
sys.exit(1 if fails else 0)
