# Outputs reference

A run writes to `<out-dir>/` (default `<fasta>_discovery/`), with the
clean, shareable results assembled under `PACKAGE/`. Per-run working data lives
in `run1/`, `run2/`, … and `downstream/`.

> For the **complete per-file, per-column catalog** (every column in `hits.tsv`,
> `paper_main_table.csv` — `database_records`/`n_genomes`/`n_organisms` —, `genome_metadata.csv`,
> the interrupted-homolog tables, the scan outputs, `control_report.json`, `run_manifest.json`, …),
> see **[REFERENCE.md → Outputs catalog](REFERENCE.md)**. This page is the orientation overview.

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
│     interrupted_homologs.tsv        (--find-interrupted only) stop-interrupted / overprinted homologs (schema below)
├── 02_sequences/                     all_hits_aa.faa / all_hits_nt.fna, unique_homologs_aa.faa
│     interrupted_homologs_domain_aa.faa, _full_orf_aa.faa, _full_orf_nt.fna  (--find-interrupted only)
│     └── per_run/runN/               hits.tsv (evidence table), hits.gff3, hits_aa.faa/hits_nt.fna,
│                                     orfs_aa.faa/orfs_nt.fna, hits_unique_aa.faa
├── 03_hmm_profile/profile.hmm        the calibrated profile HMM (submit to Pfam/CDD/VOGDB)
├── 04_alignment_phylogeny/           MSA (hits.aln.faa + stats + figure), per-hit HMM alignment
│                                     (hits_hmmalign.sto/.a2m), and the ML tree (seeds marked)
├── 05_synteny/                       clinker/ (interactive .html + static cluster_*.png), publication_figures/ (PNG/SVG/PDF —
│                                     per-cluster cluster_<id>_synteny.* plus cluster_all_homologs_synteny.* = ALL loci in one
│                                     family-anchored panel; + neighbour_gene_annotations.csv = ordered table), genbank_with_sequence/
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

## `interrupted_homologs.tsv` — column schema (`--find-interrupted`)
One row per stop-interrupted / overprinted candidate (the read-through scan).
| column | meaning |
|--------|---------|
| `contig`, `organism` | source contig accession, and its organism / phage name (joined offline from `genome_metadata.csv`; falls back to the accession if unknown) |
| `strand`, `frame` | the read-through reading frame (strand + 0/1/2) |
| `domain_nt_start`, `domain_nt_end` | genome coordinates of the matched domain (forward strand, 1-based) |
| `domain_aa_len`, `internal_stops` | domain length (aa); number of premature internal stops in it |
| `stop_nt_positions`, `stop_aa_positions` | per-stop genome coordinate(s) and aa position(s) (`;`-separated) |
| `overprinting_support` | **strong** = stop synonymous in a fully-open overlapping antisense ORF (overprinting evidence); **partial** / **none** |
| `antisense_open_frame`, `antisense_open_stops` | the antisense frame with the fewest stops over the domain, and that count (0 = fully open overlapping ORF) |
| `stop_silent_antisense` | per-stop 1/0 — the premature stop is synonymous in that antisense frame |
| `domain_bit_score`, `i_evalue` | HMM domain score / independent E-value |
| `orf_aa_len`, `aa_before_first_stop`, `aa_after_last_stop` | full read-through ORF length; intact residues before the first / after the last stop |
| `orf_nt_start`, `orf_nt_end`, `natural_stop_nt` | full-ORF genome bounds; genome coordinate of the **actual** (natural) stop codon |
| `domain_nt`, `domain_aa_with_stops` | matched-domain DNA; matched-domain protein (internal stop shown as `*`) |
| `full_orf_aa`, `full_orf_nt` | full read-through ORF — protein (premature stops `*`, terminal `*` = gene end) and nucleotide (5'→3', ends in the actual stop codon; translates back to `full_orf_aa`) |

The three `interrupted_homologs_*.faa/.fna` files carry these sequences as FASTA.

## `neighbour_gene_annotations.csv` — the ordered gene-neighbourhood table (synteny)
One row per neighbouring gene per locus, **anchored on your gene of interest** — so you
can read off and describe/label the genes bordering it. Sort by `genome_id` then
`pos_index` to walk each neighbourhood in order.
| column | meaning |
|--------|---------|
| `cluster`, `genome_id`, `organism` | which synteny cluster / source locus the gene belongs to |
| `pos_index` | **gene order relative to your gene**: `0` = your gene, `-1/-2…` upstream, `+1/+2…` downstream |
| `is_anchor` | `1` for your gene of interest (the family homolog), else `0` |
| `rel_start`, `rel_end` | gene coordinates **relative to your gene** (bp; your gene sits at 0), strand-normalised |
| `strand_vs_gene` | `+` = same orientation as your gene, `-` = opposite |
| `length_bp`, `distance_to_anchor_bp` | gene length; signed gap to your gene (− upstream, + downstream) |
| `orthogroup`, `category`, `vfam`, `function` | cross-locus orthogroup, functional category, VOGDB VFAM, and product/function |

## Which tool opens what
| File | Open in |
|------|---------|
| `*.tsv` | Excel, R, pandas |
| `*.sto` / `*.a2m` | Belvu, Jalview, HMMER (`esl-reformat`), any alignment viewer |
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
