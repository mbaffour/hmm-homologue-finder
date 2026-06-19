# Methodology

This describes the method the pipeline implements, generically. Replace the
example seed set with your own protein family; the workflow is identical.

## 1. Overview
Distant homologues of a protein family are discovered by iterative profile-HMM
search of public phage/viral sequence databases. Genome databases are translated
in all six reading frames before search, so homologues encoded by genes that
standard annotation does not predict can be recovered. For every hit, the matched
open reading frame (ORF) is reconstructed from its genomic coordinates, validated
as a genuine ORF, and exported as both nucleotide and amino-acid sequence. The
discovered, deduplicated sequences re-seed further search rounds to test for
convergence.

## 2. Seed sequences
A curated set of family protein sequences (the only manual input). Quality
matters more than quantity: full-length, genuinely related members give a sharper
profile than many fragments.

## 3. Profile HMM construction
Seeds are aligned with **MAFFT** (`--auto`) and trimmed with **trimAl**
(`automated1`). A profile HMM is built with **hmmbuild** (HMMER 3.4). The profile
is validated by self-search against the seeds; a run proceeds only if seed
recovery exceeds a threshold (default 70%).

## 4. Database search
The profile is searched (**hmmsearch**, E ≤ 1 × 10⁻⁵) against the selected
databases. Each download is integrity-checked (SHA-256) with its URL and access
date recorded. Nucleotide databases are translated in all six frames into
stop-to-stop ORFs (minimum length 30 aa) and searched with the protein profile.
Large databases are chunked with **seqkit** and processed in parallel.

## 5. ORF-validated sequence extraction
For each hit the ORF is reconstructed directly from the genomic coordinates
(1-based, strand/frame-correct) and translated; the family domain within it is
delimited by the HMM envelope (`hmmsearch --domtblout`). Recorded per hit:
full-ORF length, domain length and coverage, internal-stop count (required to be
0), and overlap with **Prodigal** gene predictions (same-strand and any-strand).
A hit passes the ORF filter if it has no internal stop codons and sits within a
genuine coding locus. Hits in annotated protein databases are captured by
accession and marked accordingly. Both nucleotide and amino-acid sequences are
exported with a per-hit evidence table.

## 6. Iterative refinement and convergence
Unique, ORF-validated domains from one round seed the next; identical databases,
parameters, and extraction are applied each round. Convergence is assessed by
comparing per-database hit counts across rounds — stable counts indicate the
detectable family has been fully recovered. (Identical hit sequences recurring
across many genomes are deduplicated before re-seeding, so the unique-sequence
count is the meaningful measure of family diversity.)

## 7. Downstream characterisation
- **Clustering** — CD-HIT (40% identity, 80% coverage).
- **Synteny** — Prodigal gene calls provide flanking-gene context; neighbourhoods
  are compared per cluster with clinker; real-sequence GenBank files are written.
- **Phylogenetics** — MAFFT alignment -> trimAl (`-gt 0.5`) -> IQ-TREE
  (ModelFinder; 1000 ultrafast bootstrap).
- **Motifs** — MEME (<=3 motifs, width 6-30 aa); scanned with FIMO.

## 8. Reproducibility
The whole workflow runs from a single command requiring only a seed FASTA. Tool
versions, database URLs, access dates, and checksums are recorded in the run's
`reproducibility.json`.

Software: HMMER 3.4, MAFFT v7.526, trimAl v1.5, Prodigal V2.6.3, seqkit v2.13.0,
CD-HIT 4.8.1, IQ-TREE 3.1.2, MEME/FIMO 5.5.9, clinker v0.0.32; genome retrieval
via NCBI Entrez and direct catalogue streaming.

## 9. Interpreting results
- **Converged** — hit counts stop growing between rounds.
- **Novel & specific** — zero hits in reviewed-protein and domain databases
  (SwissProt, Pfam, VOGDB) across rounds, with hits found only via six-frame
  translation of genome databases.
- **Validated** — every reported hit is a real ORF (no internal stops; in a
  coding locus), with both DNA and protein sequence recorded.
