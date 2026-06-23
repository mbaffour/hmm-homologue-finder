# Outputs reference

A run writes to `<out-dir>/` (default `<fasta>_discovery/`), with the
clean, shareable results assembled under `PACKAGE/`. Per-run working data lives
in `run1/`, `run2/`, … and `downstream/`.

## PACKAGE/ — the shareable result

Every folder contains a `README.txt` describing each file. Open `../report.html`
for the visual summary. The folders are numbered in reading order:
```
PACKAGE/
├── README.txt                        guide to the whole package
├── METHODS.md                        how this run was produced (methods + citations)
├── run_manifest.json                 machine-readable provenance (params, versions, calibration, seed recovery)
├── 01_summary_tables/                headline tables + database-hit bar chart  (START HERE)
│     paper_main_table.csv            one row per unique homolog (the main result)
│     hits_deduplicated.csv, hit_summary.csv, database_hit_summary.csv (+ database_hits.png/svg/pdf),
│     genome_metadata.csv, homolog_stats.csv, all_runs_hits.csv, database_summary.csv
├── 02_sequences/                     all_hits_aa.faa / all_hits_nt.fna, unique_homologs_aa.faa
│     └── per_run/runN/               hits.tsv (evidence table), hits.gff3, hits_aa.faa/hits_nt.fna,
│                                     orfs_aa.faa/orfs_nt.fna, hits_unique_aa.faa
├── 03_hmm_profile/profile.hmm        the calibrated profile HMM (submit to Pfam/CDD/VOGDB)
├── 04_alignment_phylogeny/           MSA (hits.aln.faa + stats + figure) and the ML tree (seeds marked)
├── 05_synteny/                       clinker/ (interactive), publication_figures/, genbank_with_sequence/
├── 06_database_summaries/runN_summary.tsv   per-database hit counts + provenance
├── 07_seed_qc/                       seed_recovery.csv (per-seed before/after) + seed alignment & QC tree
└── 08_scripts/                       a copy of the scripts that produced this run
```

## `hits.tsv` — column schema (one row per hit)

**Identity & provenance**
| column | meaning |
|--------|---------|
| `hit_id` | unique hit identifier (six-frame ORF id, or protein accession) |
| `genome_id` / `contig` | source genome / contig |
| `organism` | phage/organism name (NCBI) or "uncultured virus (db)" |
| `db_name` / `db_type` | database searched / nucleotide-or-protein |
| `run_label` | which iteration produced it |
| `source_url`, `source_sha256`, `accessed_at` | download provenance |
| `source_type` | `six_frame_orf` (genome hit) or `annotated_protein` (protein-DB hit) |

**Genomic location** (six-frame hits)
| `nt_start`, `nt_end`, `strand`, `frame` · `orf_nt_start`, `orf_nt_end` |

**ORF validation** (the "is it a real gene?" evidence)
| `orf_aa_len`, `domain_aa_len`, `domain_coverage`, `has_start_M`, `ends_at_stop`,
`internal_stops` (must be 0), `prodigal_concordant`, `prodigal_same_strand_pct`,
`in_coding_locus`, `prodigal_any_strand_pct`, `passes_orf_filter` |

**HMM statistics**
| `evalue`, `bit_score`, `bias_score`, `env_from`, `env_to`, `confidence_tier`, `qc_flags` |

**Sequences**
| `aa_sequence` (amino-acid domain), `nt_sequence` (matching DNA; blank for protein-DB hits) |

## Which tool opens what
| File | Open in |
|------|---------|
| `*.tsv` | Excel, R, pandas |
| `*.faa` / `*.fna` | Jalview, MEGA, AliView, BLAST, any aligner |
| `*.gff3` | IGV, JBrowse, Artemis (load with a genome FASTA) |
| `*.gbk` (GenBank) | Artemis, Geneious, UGENE, clinker, pyGenomeViz (sequence + features in one file) |
| clinker `*.html` | any web browser (interactive; "Save SVG" for figures) |
| `*.treefile` (Newick) | iTOL, FigTree, ggtree, Dendroscope |
| `*.hmm` | HMMER; submit to Pfam / NCBI CDD / VOGDB |

## Reading the result
- **Converged?** Compare `06_database_summaries/run*_summary.tsv` across rounds —
  stable counts mean the family is fully captured.
- **Novel & specific?** Zero hits in SwissProt / Pfam / VOGDB across rounds.
- **Score threshold trustworthy?** `controls/control_report.json` +
  `controls/roc_curve.svg` give sensitivity / specificity / FPR and the ROC **AUC**
  (1.0 = perfect separation) plus the Youden-optimal cutoff (advisory; the fixed
  strict 45 is kept for tiering). Summary in `run_manifest.json` → `threshold_calibration`.
- **Every hit is a real ORF** — see the ORF-validation columns; `passes_orf_filter`
  is the keep/flag decision.
- **Did every input seed come back?** `07_seed_qc/seed_recovery.csv` lists each seed
  with its best score against the initial vs final model and a recovered flag; the
  aggregate counts are in `run_manifest.json` → `seed_recovery_qc` and `METHODS.md`.
  A seed not recovered by the final model is usually a divergent outlier — consider
  dropping it or treating it as a separate sub-family.
