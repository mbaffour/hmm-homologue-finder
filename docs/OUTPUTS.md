# Outputs reference

A run writes to `<out-dir>/` (default `<fasta>_discovery/`), with the
clean, shareable results assembled under `PACKAGE/`. Per-run working data lives
in `run1/`, `run2/`, … and `downstream/`.

## PACKAGE/ — the shareable result
```
PACKAGE/
├── 01_hmm_profile/profile.hmm        the profile HMM (submit to Pfam/CDD/VOGDB)
├── 02_sequences_per_run/runN/
│     hits.tsv                          per-hit evidence table (see schema below)
│     hits.gff3                         genome-browser track of every hit
│     hits_aa.faa / hits_nt.fna    homologue domain: protein / DNA
│     orfs_aa.faa / orfs_nt.fna    full ORF context: protein / DNA
│     hits_unique_aa.faa                deduplicated domains (the next-round seed)
├── 03_database_summaries/runN_summary.tsv   per-database hit counts + provenance
├── 04_synteny_clinker/                       clinker figures + named GenBank neighbourhoods
├── 05_phylogeny/                             ML tree of the homologues (Newick + PNG/SVG)
└── 06_scripts/                               a copy of the scripts that produced this run
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
- **Converged?** Compare `03_database_summaries/run*_summary.tsv` across rounds —
  stable counts mean the family is fully captured.
- **Novel & specific?** Zero hits in SwissProt / Pfam / VOGDB across rounds.
- **Every hit is a real ORF** — see the ORF-validation columns; `passes_orf_filter`
  is the keep/flag decision.
