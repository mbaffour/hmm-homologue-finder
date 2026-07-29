#!/usr/bin/env python3
"""package_layout.py — the single source of truth for the PACKAGE/ folder layout.

Both the assembler (hmm_finder.assemble_package) and the table exporter
(export_csv) import the folder names from here so the structure stays consistent,
and `write_readmes()` drops a plain-text README.txt into every package folder
describing each file and its purpose (and one at the package root).

A registry entry is either a plain description string or, for a table, a dict:

    "hits_deduplicated.csv": {"desc": "...", "columns": {"homolog_id": "...", ...}}

`write_readmes()` renders BOTH the per-folder README.txt (with a COLUMN DICTIONARY
block under each such file) and PACKAGE/DATA_DICTIONARY.md from that one structure,
so the folder README and the data dictionary cannot drift apart.
"""
from __future__ import annotations

import csv
import fnmatch
import textwrap
from pathlib import Path

# Canonical PACKAGE/ subfolders, in reading order. Names live ONLY here.
# NOTE: the numbering is intentionally non-contiguous — "08_scripts" predates the two
# folders added after it. Renumbering would rewrite every path in the docs and in
# published links for no gain, so the new folders take 09/10 and the order below (not
# the number) is what determines reading order.
DIRS = {
    "tables":    "01_summary_tables",
    "sequences": "02_sequences",
    "hmm":       "03_hmm_profile",
    "phylo":     "04_alignment_phylogeny",
    "synteny":   "05_synteny",
    "dbsum":     "06_database_summaries",
    "seedqc":    "07_seed_qc",
    "scripts":   "08_scripts",
    "controls":  "09_controls",
    "overprint": "10_overprinting",
}
PER_RUN = "per_run"  # iteration subfolder under the sequences folder

# Columns shared by every per-stage summary table AND by pipeline_stage_summary.csv —
# defined once because the tables are deliberately one schema (see stage_summary.py).
_STAGE_COLUMNS = {
    "stage": "Two-digit stage number (00 input … 08 downstream), in pipeline order.",
    "stage_name": "Stage slug (input, model, search, validation, homologs, controls, seeds, overprint, downstream).",
    "metric": "Name of the quantity. Unique within a stage.",
    "value": "The measured value. Blank means the run did not produce it.",
    "units": "What `value` counts (hits, loci, proteins, organisms, genomes, bits, aa, %, seconds, flag, text).",
    "note": "Why the metric means what it means, and any caveat about how it must NOT be read.",
    "source_file": "The file in this run the value was read from, relative to the run directory.",
}

_STAGE_TABLE = {
    "desc": "Every per-stage table concatenated — the single table to read for a run at a "
            "glance. One row per metric, blocked by pipeline stage (00 input … 08 downstream).",
    "columns": _STAGE_COLUMNS,
}

# Same `columns` OBJECT as above on purpose: the data dictionary groups by the column list,
# so all ten stage tables render as one section instead of ten identical ones.
_STAGE_TABLE_ONE = {
    "desc": "One pipeline stage's metrics, in the shared stage schema; the stage tables "
            "concatenated are pipeline_stage_summary.csv.",
    "columns": _STAGE_COLUMNS,
}

# ---- the result tables -----------------------------------------------------------------
# Column lists below are the REAL headers written by the code that produces each file, not
# a paraphrase: paper_main_table/hits_deduplicated from export_csv._paper_table/_dedup_hits,
# seed_recovery from seed_recovery.py, interrupted_homologs from find_interrupted.py.
_HITS_DEDUP = {
    "desc": "One row per distinct homolog PROTEIN in the canonical converged run, collapsed "
            "across databases and iterations, with full provenance. Identity is the genomic "
            "LOCUS first (organism + strand + overlapping ORF interval), then the protein — "
            "never the raw aa_sequence, which is the HMM envelope slice and is re-trimmed "
            "between iterations.",
    "columns": {
        "homolog_id": "Stable id within this run (H0001…), assigned in output order (most widespread protein first).",
        "representative_organism": "Organism of the best-scoring copy of this protein.",
        "representative_genome": "Genome accession of that best-scoring copy.",
        "representative_db": "Database the representative copy came from.",
        "source_type": "How it was found: six_frame_orf (translated genome) or annotated_protein (protein database).",
        "n_organisms": "Distinct PHAGES carrying this protein, after collapsing accession aliases of the same phage. The honest breadth metric.",
        "organisms": "Semicolon-separated raw organism names (display; not deduplicated).",
        "n_loci": "Genomic gene copies of this protein — one per physical locus. Several phages can carry an identical protein, so loci >= proteins.",
        "n_genomes": "Distinct genome accessions (version suffix stripped, so a RefSeq NC_ mirror and its GenBank original count once).",
        "n_databases": "Number of DATABASE RECORDS, not independent corroboration — INPHARED redistributes RefSeq NC_ records, so one physical genome routinely appears in both.",
        "databases": "Semicolon-separated database names the protein was recovered from.",
        "n_runs": "How many search iterations recovered it.",
        "runs": "Semicolon-separated iteration labels.",
        "n_copies": "Raw hit rows collapsed into this row (records, not gene copies).",
        "domain_aa_len": "Length of the HMM-matched domain in the canonical run, in amino acids.",
        "max_domain_aa_len_any_round": "Longest domain ANY iteration called for this gene. Larger than domain_aa_len means a later, refined model trimmed it shorter; the gene is not actually short.",
        "full_length_aa": "Length of the surrounding ORF (orf_aa_len of the representative hit).",
        "domain_coverage": "Fraction of the ORF covered by the matched domain (0-1).",
        "best_evalue": "Best (lowest) HMMER E-value across all copies.",
        "best_bit_score": "Best (highest) HMMER bit score across all copies.",
        "confidence_tier": "Confidence tier of the representative hit (high_confidence / moderate / low).",
        "example_hit_id": "hit_id of the representative copy — join key into all_runs_hits.csv and homolog_stats.csv.",
        "aa_sequence": "Amino-acid sequence of the matched domain (the HMM envelope slice).",
    },
}

_PAPER_MAIN = {
    "desc": "MAIN RESULT — a projection of hits_deduplicated.csv (same rows, fewer columns, "
            "ranked by bit score) so the paper table and the homolog table can never "
            "disagree about how many homologs were found.",
    "columns": {
        "rank": "Row order, 1 = highest best_bit_score.",
        "representative_organism": "Organism of the best-scoring copy.",
        "accession": "Genome accession of that copy (hits_deduplicated.representative_genome).",
        "database": "Database it came from (hits_deduplicated.representative_db).",
        "n_loci": "Genomic gene copies of this protein.",
        "database_records": "Database RECORDS for this exact sequence (hits_deduplicated.n_copies) — NOT gene copies and NOT independent support.",
        "n_genomes": "Distinct genome accessions carrying it.",
        "n_organisms": "Distinct phages carrying it — the breadth metric.",
        "domain_aa_len": "Matched domain length in the canonical run (aa).",
        "max_domain_aa_len_any_round": "Longest domain any iteration called for this gene (aa).",
        "full_length_aa": "Surrounding ORF length (aa).",
        "domain_coverage": "Fraction of the ORF covered by the domain (0-1).",
        "best_evalue": "Best HMMER E-value.",
        "best_bit_score": "Best HMMER bit score.",
        "confidence_tier": "Confidence tier of the representative hit.",
        "example_hit_id": "hit_id of the representative copy (join key).",
    },
}

_HIT_SUMMARY = {
    "desc": "Per-iteration totals — how the search grew from one round to the next.",
    "columns": {
        "run": "Iteration number (1 = search with the model built from the input seeds).",
        "total_hits": "Validated hit rows in that iteration (one gene copy can appear in several databases).",
        "passed_filter": "Rows passing the ORF filter (passes_orf_filter = True).",
        "six_frame_hits": "Rows found by six-frame translation of a nucleotide database.",
        "protein_db_hits": "Rows found in an annotated protein database.",
        "unique_sequences": "Distinct homolog PROTEINS after locus-then-protein deduplication — the same definition used by every other table.",
        "unique_organisms": "Distinct phages (canonical organism) hit in that iteration.",
        "databases": "Semicolon-separated databases contributing hits.",
    },
}

_DB_HIT_SUMMARY = {
    "desc": "Every database searched in the canonical iteration, INCLUDING those returning "
            "zero hits — a zero is a result (an unannotated, antisense-overprinted gene is "
            "absent from the protein databases by construction). The trailing ALL row is "
            "deduplicated across databases, so it is NOT the column sum.",
    "columns": {
        "database": "Database name; the final row 'ALL (deduplicated across databases)' is the cross-database total.",
        "type": "'nucleotide (six-frame)' or 'protein'.",
        "status": "Engine status for that database in this iteration (complete / failed / skipped).",
        "hits": "Hit rows returned.",
        "strict_hits": "Hits at or above the strict bit-score threshold.",
        "unique_sequences": "Distinct homolog proteins from that database, deduplicated the same way as everywhere else.",
        "unique_organisms": "Distinct phages from that database.",
        "runtime_seconds": "Wall-clock seconds the search against that database took.",
    },
}

_GENOME_METADATA = {
    "desc": "Supplementary S1 — one row per PHYSICAL genome. Accessions are collapsed by "
            "base accession, so the same genome catalogued as NC_023589.1 and NC_023589 "
            "(RefSeq + INPHARED) is one row, not two.",
    "columns": {
        "genome_id": "Base accession (version suffix stripped) — the physical genome.",
        "accessions": "Semicolon-separated accessions that resolved to this genome.",
        "n_accessions": "How many distinct accessions that was (>1 means cross-database aliasing).",
        "organism": "Organism name as reported by the source database.",
        "host": "Host genus parsed from the organism name ('Escherichia phage X' -> Escherichia); blank when the name does not follow that form.",
        "databases": "Semicolon-separated databases the genome was found in.",
        "source_type": "Semicolon-separated hit source types for this genome (six_frame_orf / annotated_protein).",
        "n_hits": "Hit rows on this genome, across all iterations.",
    },
}

_HOMOLOG_STATS = {
    "desc": "Supplementary S3 — per-hit homology statistics, one row per hit row across ALL "
            "iterations (un-collapsed; join to hits_deduplicated.csv on example_hit_id).",
    "columns": {
        "hit_id": "Unique hit identifier (genome accession + source + ORF offset).",
        "organism": "Organism of the source genome.",
        "genome_id": "Source genome accession.",
        "db_name": "Database the hit came from.",
        "source_type": "six_frame_orf or annotated_protein.",
        "run_label": "Iteration that produced the hit.",
        "evalue": "HMMER E-value for this hit.",
        "bit_score": "HMMER bit score for this hit.",
        "orf_aa_len": "Length of the surrounding ORF (aa).",
        "domain_aa_len": "Length of the HMM-matched domain (aa).",
        "domain_coverage": "domain_aa_len / orf_aa_len (0-1).",
        "confidence_tier": "Confidence tier assigned to the hit.",
    },
}

_SEED_RECOVERY = {
    "desc": "Per-seed QC — does the model still recognise the sequences it was built from? "
            "This measures MODEL DRIFT across refinement. It does NOT say whether the "
            "seed's genomic locus was re-found by the database search (that is the family "
            "census / seed_status).",
    "columns": {
        "seed_id": "Identifier of the input seed (FASTA header, truncated at the first space).",
        "before_bit": "Best bit score of this seed against the INITIAL model.",
        "before_recovered": "True if before_bit is at or above the strict threshold.",
        "after_bit": "Best bit score of this seed against the FINAL (refined) model.",
        "after_recovered": "True if after_bit is at or above the strict threshold.",
        "status": "recovered / lost_after_refinement / gained_after_refinement / never_recovered.",
    },
}

_INTERRUPTED = {
    "desc": "(--find-interrupted) Homolog copies carrying a premature stop codon, found by "
            "read-through translation — the stop-to-stop six-frame search cannot see these. "
            "Candidate overprinted / pseudogenised genes, with the silent-stop test that "
            "distinguishes the two.",
    "columns": {
        "contig": "Source contig / genome accession the copy sits on.",
        "organism": "Organism of that genome.",
        "strand": "Strand of the homolog (+ / -).",
        "frame": "Reading frame of the homolog on that strand.",
        "domain_nt_start": "Genome coordinate where the matched domain starts.",
        "domain_nt_end": "Genome coordinate where the matched domain ends.",
        "domain_aa_len": "Length of the matched domain (aa, stops included).",
        "internal_stops": "Number of premature stop codons inside the domain.",
        "stop_nt_positions": "Genome coordinates of those stops (semicolon-separated).",
        "stop_aa_positions": "Positions of those stops within the domain protein.",
        "overprinting_support": "strong / partial / none — how consistently the premature stops are synonymous in the overlapping antisense frame.",
        "antisense_open_frame": "Which antisense frame (0/1/2) has fewest stops over the domain — the candidate overprinted gene.",
        "antisense_open_stops": "Stop codons in that antisense frame across the domain; 0 = a fully open overlapping antisense ORF.",
        "stop_silent_antisense": "Per premature stop, 1/0 = that stop is synonymous in the antisense frame (semicolon-separated, same order as stop_nt_positions).",
        "domain_bit_score": "HMMER bit score of the read-through domain match.",
        "i_evalue": "Independent E-value of the domain match.",
        "orf_aa_len": "Length of the full read-through ORF (aa).",
        "aa_before_first_stop": "Residues between the ORF start and the first premature stop.",
        "aa_after_last_stop": "Residues between the last premature stop and the natural stop.",
        "orf_nt_start": "Genome coordinate where the read-through ORF starts.",
        "orf_nt_end": "Genome coordinate where the read-through ORF ends.",
        "natural_stop_nt": "Genome coordinate of the ACTUAL terminating stop codon.",
        "domain_nt": "Nucleotide sequence of the matched domain.",
        "domain_aa_with_stops": "Domain protein with every internal stop shown as '*'.",
        "full_orf_aa": "Full read-through ORF protein; internal '*' = premature stops, trailing '*' = the natural gene end.",
        "full_orf_nt": "Full read-through ORF nucleotide (coding 5'->3', ending in the actual stop codon).",
    },
}

# Written by family_census.py and overprint_report.py — registered here so the mirror and
# the data dictionary describe them the moment they appear; absent files are skipped.
_FAMILY_CENSUS = {
    "desc": "Family size as a function of sequence-identity threshold: the union of the "
            "input seeds and the discovered homologs, clustered at 100 / 95 / 90 % identity. "
            "This answers 'did the family grow?', which the search result alone cannot — the "
            "search reports only what it FOUND, never the seed-plus-hit union.",
    "columns": {
        "identity_threshold": "Clustering identity (1.00 / 0.95 / 0.90).",
        "n_seed_proteins": "Distinct seed PROTEINS (exact duplicate headers collapsed).",
        "n_seed_headers": "Raw seed FASTA records, duplicates included.",
        "n_homolog_proteins": "Distinct homolog proteins from the canonical run.",
        "n_clusters_union": "Clusters over seeds + homologs together = the family size at this threshold.",
        "n_clusters_shared": "Clusters containing both a seed and a homolog (re-found).",
        "n_clusters_seed_only": "Clusters containing only seeds (not re-found by the search).",
        "n_clusters_new": "Clusters containing only homologs (new to the family).",
        "pct_seeds_refound": "Percentage of seed clusters that a homolog also landed in.",
        "headline": "True on the threshold quoted in the paper/report.",
        "cdhit_flags": "Exact cd-hit flags used, so the clustering is reproducible.",
        "note": "Free-text caveat for that row.",
    },
}

# Column names and order are family_census._write_members' row dict verbatim. It emits
# member_type and cluster_class_95; this entry used to say `source` and `class`, so two
# columns were documented that no file has ever contained while `member_type`, `accession`
# and `cluster_class_95` went undescribed.
_FAMILY_CENSUS_MEMBERS = {
    "desc": "Audit trail behind family_census.csv — one row per member sequence, showing "
            "which cluster it fell in at each threshold.",
    "columns": {
        "member_id": "Synthetic id (S#### seed, H#### homolog) — synthetic because pipe/space-bearing FASTA headers break cd-hit name parsing. For a homolog this is its homolog_id, so the row joins straight to hits_deduplicated.csv.",
        "member_type": "seed or homolog.",
        "label": "Original FASTA header (seed) or representative organism/genome (homolog).",
        "accession": "Nucleotide accession parsed out of the seed header / the homolog's representative genome; blank when the label carries none.",
        "aa_len": "Sequence length (aa).",
        "cluster_100": "Cluster id at 100 % identity.",
        "cluster_95": "Cluster id at 95 % identity. Blank when 95 % was not clustered (cd-hit absent) — a column named for a threshold never holds another one's answer.",
        "cluster_90": "Cluster id at 90 % identity.",
        "is_cluster_rep_95": "True if cd-hit chose this member as the representative of its 95 % cluster.",
        "cluster_class_95": "shared / seed_only / new — the class of this member's 95 % cluster.",
    },
}

# overprinted_loci.csv is NOT interrupted_homologs.tsv plus host columns. overprint_report
# carries over only 13 of that file's 26 columns; the read-through scores and the bulk
# sequence columns (stop_aa_positions, domain_bit_score, i_evalue, orf_aa_len,
# aa_before_first_stop, aa_after_last_stop, orf_nt_start, orf_nt_end, natural_stop_nt,
# domain_nt, domain_aa_with_stops, full_orf_aa, full_orf_nt) stay in interrupted_homologs.tsv,
# which ships beside it. This entry used to splice in the whole _INTERRUPTED column list and
# so published those 13 as if overprinted_loci.csv had them.
#
# The 13 inherited names are listed here and their text is pulled from _INTERRUPTED by name,
# so the two tables cannot define the same column two different ways; everything else gets
# its own description below. Order follows overprint_report.LOCI_COLS, i.e. the real file.
_OVERPRINTED_LOCI_OWN = {
    "locus_id": "Stable id of the locus within this run (OP001…), in descending domain bit score, assigned AFTER deduplication — one id per genomic locus, not per database record.",
    "accession": "Nucleotide accession of the representative record for this locus.",
    "accession_aliases": "Semicolon-separated OTHER accessions the same locus was found under (RefSeq vs GenBank, versioned vs unversioned); blank for a locus seen under a single accession. accession + accession_aliases is the complete record set.",
    "host_gene": "Gene name of the antisense CDS the homolog is overprinted inside.",
    "host_product": "Product/function annotation of that host gene.",
    "host_locus_tag": "Locus tag of the host gene in the source record.",
    "host_protein_id": "Protein accession of the host gene.",
    "host_gene_start": "Genome coordinate where the host gene starts.",
    "host_gene_end": "Genome coordinate where the host gene ends.",
    "host_gene_strand": "Strand of the host gene (opposite the homolog, by definition of antisense overprinting).",
    "host_gene_aa_len": "Length of the host protein (aa).",
    "host_annotation_source": "Where the host annotation came from (NCBI record, catalogue, or Prodigal call); 'none' when no annotation was available — the row is kept either way.",
    "host_category": "Functional category assigned to the host gene.",
    "nested_fully": "True when the homolog domain lies entirely inside the host gene.",
    "overlap_bp": "Nucleotides of overlap between the domain and the host gene.",
    "overlap_pct_of_domain": "overlap_bp as a percentage of the domain length.",
    "antisense_orf_is_open": "1 = the antisense frame is open across the whole domain and the ORF extent below was measured; 0 = it is not open, so there is no ORF and the antisense_orf_* coordinates are blank; blank = not determined (no contig sequence, or no antisense frame recorded).",
    "antisense_orf_nt_start": "Start of the computed open antisense ORF over the domain. Blank unless antisense_orf_is_open = 1.",
    "antisense_orf_nt_end": "End of that computed antisense ORF. Blank unless antisense_orf_is_open = 1.",
    "antisense_orf_aa_len": "Length of that computed antisense ORF (aa). Blank unless antisense_orf_is_open = 1.",
    "antisense_orf_matches_host_gene": "True when the computed antisense ORF coincides with the annotated host CDS — the single strongest line of evidence that the overlap is a real, expressed gene.",
    "figure": "Filename of this locus's diagram in this folder (blank if the figure cap was hit).",
}

_OVERPRINTED_LOCI_ORDER = (
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
)

_OVERPRINTED_LOCI = {
    "desc": "One row per interrupted/overprinted LOCUS (deduplicated across database "
            "records) WITH the antisense host gene it is printed inside — the evidence "
            "behind the overprinting claim. Carries the coordinate and silent-stop columns "
            "of interrupted_homologs.tsv; that file keeps the read-through scores and the "
            "domain/ORF sequences, which are NOT repeated here.",
    "columns": {c: (_OVERPRINTED_LOCI_OWN.get(c) or _INTERRUPTED["columns"].get(c)
                    or "See interrupted_homologs.tsv.")
                for c in _OVERPRINTED_LOCI_ORDER},
}

# NOT the seven-column stage schema. overprint_report.SUMMARY_COLS is (metric, value, note):
# stage_summary._stage07 READS this file and re-emits its rows into the stage schema, which
# is a different thing from this file being in it. Documenting the stage schema here
# published four columns (stage, stage_name, units, source_file) that overprinting_summary.csv
# has never had.
_OVERPRINTING_SUMMARY = {
    "desc": "Family-level overprinting rollup: how many homolog copies are stop-interrupted, "
            "how well the stops are explained by an overlapping antisense gene, and what "
            "those host genes are. Three columns, NOT the seven-column stage schema — "
            "stage 07 of pipeline_stage_summary.csv is built by re-reading these rows.",
    "columns": {
        "metric": "Name of the quantity.",
        "value": "The measured value.",
        "note": "Definition of the metric and the caveat on how it must NOT be read; every row carries one, so this file is readable without the methods section.",
    },
}

# Per-folder purpose + per-file descriptions. Files are matched by exact name,
# then by glob pattern; anything unmatched gets a generic by-extension label.
_REGISTRY = {
    DIRS["tables"]: {
        "purpose": "The headline result tables and the per-database hit chart. Start here.",
        "files": {
            "paper_main_table.csv": _PAPER_MAIN,
            "hits_deduplicated.csv": _HITS_DEDUP,
            "hit_summary.csv": _HIT_SUMMARY,
            "database_hit_summary.csv": _DB_HIT_SUMMARY,
            "database_hits.png": "Bar chart of hits per database (grey = searched, 0 hits).",
            "database_hits.svg": "Bar chart, editable vector (Inkscape/Illustrator).",
            "database_hits.pdf": "Bar chart, print-ready vector.",
            "genome_metadata.csv": _GENOME_METADATA,
            "homolog_stats.csv": _HOMOLOG_STATS,
            "all_runs_hits.csv": "The complete, un-collapsed hit table across all iterations (every hit, every column; schema in docs/OUTPUTS.md). Counting distinct aa_sequence strings HERE is wrong — the sequence is the HMM envelope slice and is re-trimmed between iterations; use hits_deduplicated.csv.",
            "database_summary.csv": "Raw engine per-database summary across all iterations (one row per database per iteration, with status, hit counts, runtime and source provenance).",
            "pipeline_stage_summary.csv": _STAGE_TABLE,
            "stage*_summary.csv": _STAGE_TABLE_ONE,
            "family_census.csv": _FAMILY_CENSUS,
            "family_census_members.csv": _FAMILY_CENSUS_MEMBERS,
            "overprinted_loci.csv": _OVERPRINTED_LOCI,
            "overprinting_summary.csv": _OVERPRINTING_SUMMARY,
            "interrupted_homologs.tsv": _INTERRUPTED,
        },
    },
    DIRS["sequences"]: {
        "purpose": "All discovered sequences — combined multi-FASTAs plus per-iteration files.",
        "files": {
            "all_hits_aa.faa": "Every validated hit, amino-acid (all iterations).",
            "all_hits_nt.fna": "Every validated hit, nucleotide (all iterations).",
            "unique_homologs_aa.faa": "One sequence per unique homolog, rich headers (organism, #organisms, #databases).",
            "interrupted_homologs_domain_aa.faa": "(--find-interrupted) Protein of each stop-interrupted/overprinted homolog's matched domain, with every internal stop shown as '*'.",
            "interrupted_homologs_full_orf_aa.faa": "(--find-interrupted) Full read-through ORF protein for each interrupted homolog — premature stops kept as '*', terminal '*' = the natural gene end.",
            "interrupted_homologs_full_orf_nt.fna": "(--find-interrupted) Full read-through ORF nucleotide (coding 5'->3', ending in the actual stop codon triplet; translates back to the full ORF protein).",
            PER_RUN: "Per-iteration sequences and the full per-hit evidence table (see per_run/README.txt).",
        },
    },
    PER_RUN: {
        "purpose": "Outputs for a single search iteration.",
        "files": {
            "hits.tsv": "Full per-hit evidence table (column schema in docs/OUTPUTS.md).",
            "hits.csv": "Same table, comma-separated for Excel.",
            "hits.gff3": "Genome-browser track of every hit (IGV / JBrowse / Artemis).",
            "hits_aa.faa": "Homolog domain, amino-acid.",
            "hits_nt.fna": "Homolog domain, nucleotide.",
            "orfs_aa.faa": "Full surrounding ORF, amino-acid.",
            "orfs_nt.fna": "Full surrounding ORF, nucleotide.",
            "hits_unique_aa.faa": "Deduplicated domains used to seed the next iteration.",
        },
    },
    DIRS["hmm"]: {
        "purpose": "The calibrated profile HMM — the model the whole run is built on.",
        "files": {
            "profile.hmm": "Profile HMM (most-complete run). Search with HMMER; submit to Pfam / NCBI CDD / VOGDB.",
        },
    },
    DIRS["phylo"]: {
        "purpose": "Multiple-sequence alignment of the homologs and the maximum-likelihood tree.",
        "files": {
            "hits.treefile": "ML tree (Newick) of the homologs — the input seeds are included for context and labelled with a '_seed' suffix. Open in iTOL/FigTree.",
            "hits.contree": "Bootstrap consensus tree (Newick).",
            "hits_tree.png": "Rendered tree figure (300 dpi).",
            "hits_tree.svg": "Rendered tree, editable vector.",
            "hits_tree.pdf": "Rendered tree, print-ready vector.",
            "hits_tree_homologs_only.*": "Tree with the seed tips pruned (homologs only).",
            "hits.aln.faa": "Multiple sequence alignment of the homologs (MAFFT).",
            "hits.aln.trim.faa": "Trimmed alignment (trimAl) — the tree input.",
            "hits_hmmalign.sto": "Per-hit alignment of every homolog to the family HMM (hmmalign, Stockholm) — shows each hit's match states vs insertions relative to the model.",
            "hits_hmmalign.a2m": "Same per-hit HMM alignment in A2M (aligned FASTA; UPPERCASE = match columns, lowercase/'.' = insertions).",
            "hits.aln.stats.json": "Alignment quality: length, gap %, conserved columns, mean pairwise identity.",
            "alignment_figure.*": "ClustalX-coloured alignment figure (PNG/SVG/PDF).",
            "hits.iqtree": "IQ-TREE report (model selection, parameters, support).",
            "hits.log": "IQ-TREE run log.",
            "tree_input.faa": "Exact sequences fed to the aligner (organism-labelled).",
        },
    },
    DIRS["synteny"]: {
        "purpose": "Gene-neighbourhood (synteny) comparisons and real-sequence GenBank files.",
        "files": {
            "clinker": "clinker gene-neighbourhood comparison — interactive cluster_*.html (open in any browser) plus, when a headless browser is installed, a static cluster_*.png per cluster (the 'Static PNG' column in clinker/index.html).",
            "publication_figures": "Publication synteny panels (PNG/SVG/PDF) + neighbour_gene_annotations.csv — the ordered gene-neighbourhood table: each neighbour's pos_index (order; 0 = your gene, - upstream, + downstream), position relative to your gene (rel_start/rel_end, distance_to_anchor_bp), strand vs. your gene, and function, so you can describe/label the bordering genes.",
            "genbank_with_sequence": "Real-sequence GenBank neighbourhoods, named by phage (Artemis/Geneious/UGENE), each with a <name>_genome_map.png/.svg marking the gene of interest (the HMM hit, red) among its neighbours.",
        },
    },
    DIRS["dbsum"]: {
        "purpose": "The engine's per-database summary for each iteration (status, hit counts, provenance).",
        "files": {
            "run*_summary.tsv": "Per-database results for one iteration: database, status, hit/strict counts, runtime, provenance.",
        },
    },
    DIRS["seedqc"]: {
        "purpose": "Quality control of the INPUT seeds (before committing to the search).",
        "files": {
            "seed_recovery.csv": _SEED_RECOVERY,
            "seed_status.csv": {
                "desc": "One row per DISTINCT seed protein: whether the search re-found its "
                        "locus, and at which identity threshold. This is the 'was this seed "
                        "re-found?' table; seed_recovery.csv is the different question of "
                        "whether the model still scores the seed.",
                "columns": {
                    "seed_id": "Identifier of the distinct seed protein.",
                    "duplicate_of": "seed_id this record duplicates exactly (blank if it is the representative).",
                    "n_seed_headers": "How many input FASTA records carry this exact protein.",
                    "organism": "Organism parsed from the seed header.",
                    "accession": "Nucleotide/protein accession parsed from the seed header.",
                    "accession_prefix": "Accession prefix (NC_, MT, OU, …).",
                    "accession_class": "refseq / genbank / metagenome / unresolved.",
                    "aa_len": "Seed protein length (aa).",
                    "recovered_by_hmm": "True if the final model still scores this seed above threshold.",
                    "after_bit": "That seed's bit score against the final model.",
                    "refound_100": "True if a discovered homolog clusters with it at 100 % identity.",
                    "refound_95": "…at 95 % identity (the headline threshold).",
                    "refound_90": "…at 90 % identity.",
                    "best_homolog_id_95": "homolog_id of the homolog it clustered with at 95 %.",
                    "status": "refound / missed.",
                },
            },
            "tree_input.faa": "The seed sequences as aligned (organism-labelled).",
            "hits.aln.faa": "Alignment of just the input seeds (sanity-check the set).",
            "hits.treefile": "Seed-only QC tree (Newick).",
            "hits_tree.*": "Seed-only QC tree figure.",
        },
    },
    DIRS["scripts"]: {
        "purpose": "A copy of the exact scripts that produced this run (for reproducibility).",
        "files": {},
    },
    DIRS["controls"]: {
        "purpose": "The specificity evidence: positive/negative control sets, the ROC curve, "
                   "and the six-frame decoy false-discovery-rate test. The report quotes "
                   "these numbers — this folder is where they come from.",
        "files": {
            "controls_summary.csv": {
                "desc": "One row per control set — the table behind the sensitivity/specificity claim.",
                "columns": {
                    "name": "Control set name.",
                    "role": "positive (must be detected) or negative (must NOT be detected).",
                    "desc": "What the set contains.",
                    "n_sequences": "Sequences in the set.",
                    "n_hits": "Sequences the model detected.",
                    "n_hits_strict": "Detections at or above the strict bit-score threshold.",
                    "hit_rate_pct": "n_hits / n_sequences as a percentage.",
                    "min_score": "Lowest bit score in the set (0 when nothing was detected).",
                    "max_score": "Highest bit score in the set.",
                    "mean_score": "Mean bit score over detected sequences.",
                    "pass": "True when the set behaved as its role requires.",
                },
            },
            "control_report.json": "Machine-readable control results: sensitivity, specificity, false-positive rate, per-set detail, and the ROC block (auc, strict_threshold, Youden optimum and whether that optimum is even defined).",
            "sixframe_decoy_control.json": {
                "desc": "The six-frame decoy FDR test — reversed six-frame ORFs from the SAME "
                        "genomes (identical length and composition, no homology). It is the "
                        "only control that measures the false discovery rate of a six-frame "
                        "search, which protein-database negative controls cannot.",
                "columns": {
                    "n_decoy_sequences": "Decoy ORFs scored.",
                    "n_sixframe_orfs_total": "Size of the real six-frame search space.",
                    "sampled_fraction": "n_decoy_sequences / n_sixframe_orfs_total.",
                    "threshold": "Strict bit-score cutoff the FDR is quoted at.",
                    "best_decoy_bit_score": "Best-scoring decoy — the noise ceiling.",
                    "weakest_true_positive_bit_score": "Weakest accepted real hit — the signal floor.",
                    "gap_bits": "Signal floor minus noise ceiling; > 0 means the two do not overlap.",
                    "clean_separation": "True when no decoy scores as high as any accepted hit.",
                    "decoys_at_or_above_threshold": "Decoys reaching the strict threshold.",
                    "expected_decoys_in_full_search_space": "That count scaled from the sample to the whole search space.",
                    "true_positives_at_or_above_threshold": "Real hits at or above the threshold.",
                    "empirical_fdr": "expected decoys / (expected decoys + true positives).",
                    "hmmsearch_filters": "Reporting filters used — deliberately opened so weak decoys are not censored out of the comparison.",
                },
            },
            "roc_curve.png": "ROC curve over the control sets (raster, 300 dpi).",
            "roc_curve.svg": "ROC curve, editable vector.",
            "roc_curve.pdf": "ROC curve, print-ready vector.",
            "ctrl_*.tbl": "Raw HMMER per-sequence output for one control set.",
            "shuffled_ctrl.faa": "The amino-acid-shuffled seeds used as a composition-matched negative control.",
            "sixframe_decoy.faa": "The sampled reversed six-frame decoy ORFs.",
            "sixframe_decoy.tbl": "Raw HMMER output over the decoys.",
        },
    },
    DIRS["overprint"]: {
        "purpose": "Overprinting evidence: each interrupted homolog copy shown against the "
                   "antisense gene it is printed inside, plus the family-level rollup.",
        "files": {
            "overprinted_loci.csv": _OVERPRINTED_LOCI,
            "overprinting_summary.csv": _OVERPRINTING_SUMMARY,
            "overprinting_overview.png": "One line per locus (domain length, coloured by support, labelled with the host product) — the whole family's overprinting at a glance.",
            "overprinting_overview.svg": "Overview figure, editable vector.",
            "overprinting_overview.pdf": "Overview figure, print-ready vector.",
            "OP*_antisense.png": "Per-locus diagram: the homolog (sense) drawn over the antisense host gene it is nested in, with each premature stop marked.",
            "OP*_antisense.svg": "Per-locus diagram, editable vector.",
            "OP*_antisense.pdf": "Per-locus diagram, print-ready vector.",
        },
    },
}

_EXT_DESC = {
    ".csv": "table (comma-separated)", ".tsv": "table (tab-separated)",
    ".faa": "protein FASTA", ".fna": "nucleotide FASTA", ".fasta": "FASTA",
    ".png": "figure (raster, 300 dpi)", ".svg": "figure (editable vector)",
    ".pdf": "figure (print vector)", ".hmm": "HMM profile (HMMER)",
    ".gbk": "GenBank record", ".gff3": "genome-browser track",
    ".json": "machine-readable data", ".treefile": "phylogenetic tree (Newick)",
    ".contree": "bootstrap consensus tree (Newick)", ".nwk": "phylogenetic tree (Newick)",
    ".md": "documentation", ".html": "interactive web page", ".log": "log file",
    ".txt": "text", ".py": "Python script",
    ".sto": "Stockholm alignment (per-hit HMM alignment)",
    ".a2m": "A2M alignment (aligned FASTA; UPPERCASE = HMM match columns)",
}


# Folder reading order for the top-level README and the data dictionary. Kept separate from
# DIRS because the folder NUMBERS are not contiguous (see the note on DIRS) — this tuple, not
# the prefix, is what "in reading order" means.
_ORDER = ("tables", "sequences", "hmm", "phylo", "synteny", "dbsum", "seedqc", "scripts",
          "controls", "overprint")

DATA_DICTIONARY = "DATA_DICTIONARY.md"


def _entry(name: str, is_dir: bool, files_map: dict):
    """The registry entry for a file — a plain description string, or a
    {"desc":…, "columns":{…}} dict for a table. Exact name wins over a glob pattern so a
    specific file can override a family (pipeline_stage_summary.csv vs stage*_summary.csv)."""
    if name in files_map:
        return files_map[name]
    for pat, desc in files_map.items():
        if ("*" in pat or "?" in pat) and fnmatch.fnmatch(name, pat):
            return desc
    if is_dir:
        return "folder"
    return _EXT_DESC.get(Path(name).suffix.lower(), "output file")


def _describe(name: str, is_dir: bool, files_map: dict) -> str:
    """One-line description, whichever form the registry entry takes."""
    e = _entry(name, is_dir, files_map)
    return e.get("desc", "table") if isinstance(e, dict) else e


def _columns_of(name: str, is_dir: bool, files_map: dict) -> dict:
    """{column: meaning} for a registered table, {} for everything else."""
    e = _entry(name, is_dir, files_map)
    cols = e.get("columns") if isinstance(e, dict) else None
    return cols if isinstance(cols, dict) else {}


def _actual_columns(path: Path) -> list:
    """Header of a delimited file on disk, [] if it is not one / cannot be read.

    Used to cross-check the registry against reality: a documented column that vanished, or
    a new column nobody documented, is exactly the drift a data dictionary exists to catch,
    and silently publishing a dictionary that disagrees with the file is worse than none."""
    path = Path(path)                       # tolerate a str from a caller
    suf = path.suffix.lower()
    if suf not in (".csv", ".tsv"):
        return []
    try:
        with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
            row = next(csv.reader(fh, delimiter="\t" if suf == ".tsv" else ","), [])
        return [c.strip() for c in row if str(c).strip()]
    except Exception:
        return []


def _wrap(text: str, width: int = 96, indent: str = "") -> list:
    return textwrap.wrap(str(text), width=width,
                         initial_indent=indent, subsequent_indent=indent) or [indent.rstrip()]


def _folder_readme(folder: Path, key: str) -> str:
    info = _REGISTRY.get(key, {})
    purpose = info.get("purpose", "Run outputs.")
    files_map = info.get("files", {})
    lines = [f"{folder.name}/  —  {purpose}", "=" * 70, "",
             "Contents (file : what it is):", ""]
    entries = sorted([p for p in folder.iterdir() if p.name != "README.txt"],
                     key=lambda p: (not p.is_dir(), p.name.lower()))
    if not entries:
        lines.append("  (empty)")
    for p in entries:
        tag = p.name + ("/" if p.is_dir() else "")
        lines.append(f"  {tag:<34} {_describe(p.name, p.is_dir(), files_map)}")

    # COLUMN DICTIONARY — only for files that are actually here AND have a registered
    # column list, so a folder of figures does not grow a dictionary section.
    documented = [p for p in entries if _columns_of(p.name, p.is_dir(), files_map)]
    if documented:
        lines += ["", "-" * 70, "COLUMN DICTIONARY", "-" * 70]
        for p in documented:
            cols = _columns_of(p.name, p.is_dir(), files_map)
            lines += ["", f"{p.name}"]
            for col, meaning in cols.items():
                wrapped = _wrap(meaning, width=96, indent=" " * 34)
                if len(col) >= 30:
                    # a long column name would run straight into its description with no
                    # gap; give it its own line instead of producing an unreadable join
                    lines.append(f"  {col}")
                    lines += wrapped
                else:
                    lines.append(f"  {col:<30}" + wrapped[0].lstrip())
                    lines += wrapped[1:]
            # Cross-check BOTH ways. Reporting only present-but-undocumented is how three
            # entries shipped describing columns their file does not have (13 phantom
            # columns on overprinted_loci.csv, 4 on overprinting_summary.csv): a dictionary
            # entry for a column that does not exist is never contradicted by the file.
            actual = _actual_columns(p)
            extra = [c for c in actual if c not in cols]
            if extra:
                lines += _wrap("NOTE: columns present in the file but not described here: "
                               + ", ".join(extra), width=96, indent="  ")
            # Guarded on `actual`: a figure, a folder or an unreadable table returns [] and
            # must not be reported as "every documented column is missing".
            absent = [c for c in cols if c not in actual] if actual else []
            if absent:
                lines += _wrap("NOTE: columns described here but NOT present in the file: "
                               + ", ".join(absent), width=96, indent="  ")
    lines += ["", "Generated by hmm-homologue-finder. See ../README.txt for the whole package,",
              f"or ../{DATA_DICTIONARY} for every table's columns in one place."]
    return "\n".join(lines) + "\n"


def _top_readme(pkg: Path) -> str:
    order = [("README.txt", "this file"),
             (DATA_DICTIONARY, "every table in this package, column by column"),
             ("METHODS.md", "how this run was produced (methods + citations)"),
             ("run_manifest.json", "machine-readable provenance (parameters, tool versions, calibration, seed recovery)")]
    lines = ["PACKAGE/  —  self-contained, shareable results of an hmm-homologue-finder run",
             "=" * 74, "",
             "Open  ../report.html  for a one-page visual summary of this run.",
             "", "At this level:", ""]
    for n, d in order:
        if (pkg / n).exists() or n in ("README.txt", DATA_DICTIONARY):
            lines.append(f"  {n:<22} {d}")
    lines += ["", "Folders (in reading order):", ""]
    for key in _ORDER:
        d = pkg / DIRS[key]
        if d.exists():
            lines.append(f"  {DIRS[key] + '/':<26} {_REGISTRY[DIRS[key]]['purpose']}")
    lines += ["", "Every folder has its own README.txt listing each file's purpose."]
    return "\n".join(lines) + "\n"


def _dictionary_md(pkg: Path) -> str:
    """Render PACKAGE/DATA_DICTIONARY.md from the SAME registry the folder READMEs use, so
    the two can never disagree about what a column means. Only files present in this package
    are listed — a dictionary describing files that did not ship is misinformation."""
    lines = ["# Data dictionary", "",
             "Every table in this package, column by column. Generated by "
             "hmm-homologue-finder from `scripts/package_layout.py`, which is also what "
             "writes each folder's `README.txt` — the two are the same source.", "",
             "Files that this run did not produce are not listed.", ""]
    for key in _ORDER:
        folder = pkg / DIRS[key]
        if not folder.is_dir():
            continue
        files_map = _REGISTRY.get(DIRS[key], {}).get("files", {})
        present = sorted((p for p in folder.iterdir() if p.name != "README.txt"),
                         key=lambda p: p.name.lower())
        # Group by COLUMN LIST identity: the ten stage tables deliberately share one
        # `columns` object, and printing the same column list ten times would bury the
        # tables that really are different. Grouping is on `id(cols)`, so only tables that
        # literally share the one object are folded together — a table that merely looks
        # similar keeps its own dict and its own section.
        groups, order = {}, []
        for p in present:
            cols = _columns_of(p.name, p.is_dir(), files_map)
            if not cols:
                continue
            e = _entry(p.name, p.is_dir(), files_map)
            k = id(cols)
            if k not in groups:
                groups[k] = {"entry": e, "names": [], "paths": []}
                order.append(k)
            groups[k]["names"].append(p.name)
            groups[k]["paths"].append(p)
        if not order:
            continue
        lines += [f"## `{DIRS[key]}/`", "", _REGISTRY[DIRS[key]]["purpose"], ""]
        for k in order:
            g = groups[k]
            lines += [f"### {', '.join('`%s`' % n for n in g['names'])}", "",
                      str(g["entry"].get("desc", "")), "",
                      "| column | meaning |", "| --- | --- |"]
            cols = g["entry"].get("columns", {})
            for col, meaning in cols.items():
                lines.append(f"| `{col}` | {str(meaning).replace('|', chr(92) + '|')} |")
            lines.append("")
            # Cross-check BOTH directions — see _folder_readme. `actual` is the UNION over
            # the grouped files, so a column is only called absent when NO file in the group
            # has it; and the whole check is skipped when the union is empty (nothing in the
            # group was a readable delimited file), which would otherwise flag every
            # documented column at once.
            actual = {c for p in g["paths"] for c in _actual_columns(p)}
            undocumented = sorted(actual - set(cols))
            if undocumented:
                lines += ["> Columns present in the file(s) but not described above: "
                          + ", ".join(f"`{c}`" for c in undocumented), ""]
            absent = [c for c in cols if c not in actual] if actual else []
            if absent:
                lines += ["> Columns described above but NOT present in the file(s): "
                          + ", ".join(f"`{c}`" for c in absent), ""]
    lines += ["---", "",
              "Counting note that applies to every homolog table: a homolog's identity is its "
              "genomic **locus** (organism + strand + overlapping ORF interval), not its "
              "amino-acid string. The sequence is the HMM envelope slice and is re-trimmed "
              "whenever the model is refined, so counting distinct sequences in "
              "`all_runs_hits.csv` over-counts the family.", ""]
    return "\n".join(lines)


def write_data_dictionary(pkg: Path) -> Path | None:
    """Write PACKAGE/DATA_DICTIONARY.md. Returns the path, or None if it could not be
    written. Never raises — a missing data dictionary must not cost the user the package."""
    try:
        pkg = Path(pkg)
        if not pkg.is_dir():
            return None
        p = pkg / DATA_DICTIONARY
        p.write_text(_dictionary_md(pkg), encoding="utf-8")
        return p
    except Exception:
        return None


def write_readmes(pkg: Path, log=None) -> None:
    """Write README.txt into the package root and every known subfolder, plus the
    package-wide DATA_DICTIONARY.md. Never raises."""
    pkg = Path(pkg)
    try:
        if not pkg.exists():
            return
        # subfolder READMEs (incl. the nested sequences/per_run/runN folders)
        for key in DIRS.values():
            folder = pkg / key
            if folder.is_dir():
                (folder / "README.txt").write_text(_folder_readme(folder, key), encoding="utf-8")
        seqdir = pkg / DIRS["sequences"] / PER_RUN
        if seqdir.is_dir():
            (seqdir / "README.txt").write_text(_folder_readme(seqdir, PER_RUN), encoding="utf-8")
            for run in sorted(p for p in seqdir.iterdir() if p.is_dir()):
                (run / "README.txt").write_text(_folder_readme(run, PER_RUN), encoding="utf-8")
        # Written BEFORE the top-level README so the README's "at this level" listing sees it.
        dd = write_data_dictionary(pkg)
        (pkg / "README.txt").write_text(_top_readme(pkg), encoding="utf-8")
        if log:
            log(f"  wrote README.txt in PACKAGE/ and {sum(1 for _ in pkg.rglob('README.txt'))-1} subfolder(s)"
                + (f" + {DATA_DICTIONARY}" if dd else ""))
    except Exception as e:
        if log:
            log(f"  (package READMEs skipped: {e})")
