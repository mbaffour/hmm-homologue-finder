"""
databases/builtin.py — Pre-configured database entries.

Every database with a direct FASTA URL streams automatically — no manual
download. The search handler pipes data through:

    curl -sL <url> | gunzip | hmmsearch --tblout results.tbl <hmm> -

Nucleotide databases add inline 6-frame translation in the pipe.
Multi-file databases (RefSeq wildcard URLs) expand via NCBI FTP listing.

Databases with special formats (HMM libraries, tar archives) auto-download
and prepare themselves on first use.
"""

from typing import Optional


def _db(
    name: str,
    db_type: str,
    streaming: bool,
    optional: bool,
    size_hint: str,
    notes: str = "",
    relevance: str = "",
    path: Optional[str] = None,
    url: Optional[str] = None,
    download_url: Optional[str] = None,
    enabled: bool = True,
    est_time: str = "",
    search_mode: str = "hmmsearch",
    setup_handler: Optional[str] = None,
    annotation_url: Optional[str] = None,
    release: str = "",
    checksum_md5: str = "",
    annotation_checksum_md5: str = "",
) -> dict:
    return {
        "name": name,
        "type": db_type,
        "path": path,
        "url": url,
        "download_url": download_url,
        "streaming": streaming,
        "enabled": enabled,
        "optional": optional,
        "size_hint": size_hint,
        "est_time": est_time,
        "notes": notes,
        "relevance": relevance,
        "search_mode": search_mode,
        "setup_handler": setup_handler,
        "annotation_url": annotation_url,
        "release": release,
        "checksum_md5": checksum_md5,
        "annotation_checksum_md5": annotation_checksum_md5,
        "last_searched": None,
        "last_hit_count": None,
    }


BUILTIN_DATABASES: list[dict] = [

    # ══════════════════════════════════════════════════════════════════════
    # TIER 1 — Fast, small, streams instantly. Pre-ticked by default.
    # ══════════════════════════════════════════════════════════════════════

    _db(
        name="INPHARED genomes",
        db_type="nucleotide",
        streaming=True,
        optional=False,
        size_hint="~497 MB",
        est_time="8-12 min",
        download_url="https://millardlab-inphared.s3.climb.ac.uk/1Jan2024_genomes.fa.gz",
        notes=(
            "INPHARED complete phage genome collection (~20K genomes). "
            "Nucleotide — 6-frame translated on the fly. "
            "Best source for discovering novel phage gene family members."
        ),
        relevance=(
            "Core reference-phage discovery set. Use when reviewers need curated phage genome evidence "
            "and when unannotated ORFs may matter."
        ),
    ),
    _db(
        name="INPHARED proteins",
        db_type="protein",
        streaming=True,
        optional=True,
        size_hint="~215 MB",
        est_time="~15 sec",
        download_url="https://millardlab-inphared.s3.climb.ac.uk/1Jan2024_vConTACT2_proteins.faa.gz",
        notes=(
            "INPHARED vConTACT2 representative proteins (~1.9M). "
            "Fast but only a subset — use INPHARED genomes for comprehensive discovery."
        ),
        relevance=(
            "Fast curated protein check. Useful as a quick proteome baseline, but less exhaustive than genome six-frame search."
        ),
    ),
    _db(
        name="SwissProt",
        db_type="protein",
        streaming=True,
        optional=False,
        size_hint="~90 MB",
        est_time="1-3 min",
        download_url=(
            "https://ftp.uniprot.org/pub/databases/uniprot/"
            "current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
        ),
        notes="~570K manually reviewed proteins from UniProt. Fast to stream.",
        relevance=(
            "Specificity/control database. Helps show whether hits are viral-family specific or also appear in reviewed general proteins."
        ),
    ),
    _db(
        name="RefSeq viral proteins",
        db_type="protein",
        streaming=True,
        optional=False,
        size_hint="~300 MB",
        est_time="2-5 min",
        url="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/",
        download_url="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.*.protein.faa.gz",
        notes="All NCBI viral proteins. Multi-file, each part streamed separately.",
        relevance=(
            "Reference viral proteome search. Reviewers expect this for annotated RefSeq viral protein evidence."
        ),
    ),

    # ══════════════════════════════════════════════════════════════════════
    # TIER 2 — Streamable, larger, takes longer. Enabled but not pre-ticked.
    # ══════════════════════════════════════════════════════════════════════

    _db(
        name="RefSeq viral genomes",
        db_type="nucleotide",
        streaming=True,
        optional=False,
        size_hint="~5 GB",
        est_time="10-30 min",
        url="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/",
        download_url="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.*.genomic.fna.gz",
        notes=(
            "All NCBI viral genomes (nucleotide). 6-frame translated on the fly. "
            "Slower than protein DB but finds unannotated ORFs."
        ),
        relevance=(
            "Reference viral genome search. Strong for finding missed or unannotated ORFs in NCBI viral genomes."
        ),
    ),
    _db(
        name="Gut Phage Database (GPD)",
        db_type="nucleotide",
        streaming=True,
        optional=True,
        size_hint="~1.5 GB",
        est_time="15-30 min",
        download_url="https://zenodo.org/records/6503062/files/GPD_sequences.fa.gz",
        notes="Gut Phage Database — ~142K gut phage genomes. 6-frame translated on the fly.",
        relevance=(
            "Gut-phage breadth set. Use when claims involve gut phages or diversity beyond curated reference genomes."
        ),
    ),
    _db(
        name="GVD-AVrC",
        db_type="nucleotide",
        streaming=True,
        optional=True,
        size_hint="~5 GB",
        est_time="30-90 min",
        download_url="https://zenodo.org/records/11426065/files/AVrC_allrepresentatives.fasta.gz",
        notes=(
            "Aggregated Gut Viral Catalogue (AVrC) — successor to MGV. "
            "~300K gut viral representative genomes. Large, stream with patience."
        ),
        relevance=(
            "Large gut/environmental viral catalogue. Use to strengthen broad diversity claims and address reviewer concerns about viral dark matter."
        ),
    ),

    # ══════════════════════════════════════════════════════════════════════
    # TIER 3 — Very large, streamable but slow.
    # ══════════════════════════════════════════════════════════════════════

    _db(
        name="RefSeq bacterial proteins",
        db_type="protein",
        streaming=True,
        optional=True,
        size_hint="~80 GB total",
        est_time="15-30 min",
        url="https://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/",
        download_url="https://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.*.protein.faa.gz",
        notes="All NCBI bacterial proteins (~973 files). Streams sequentially, ~1s per file.",
        relevance=(
            "Host/background specificity check. Useful for showing whether candidate homologs are phage-enriched rather than broadly bacterial."
        ),
    ),

    # ══════════════════════════════════════════════════════════════════════
    # TIER 4 — Auto-download on first use. Need prep but the app handles it.
    # ══════════════════════════════════════════════════════════════════════

    _db(
        name="Pfam (sequences)",
        db_type="protein",
        streaming=True,
        optional=True,
        size_hint="~6 GB",
        est_time="20-60 min",
        download_url="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.fasta.gz",
        notes=(
            "Pfam seed+full sequences as FASTA — streamed through hmmsearch. "
            "Finds which Pfam families your protein matches."
        ),
        relevance=(
            "Domain-family context. Useful for checking whether the HMM retrieves known Pfam-associated proteins or unexpected relatives."
        ),
    ),
    _db(
        name="Pfam (domain scan)",
        db_type="protein",
        streaming=False,
        optional=True,
        size_hint="~300 MB",
        est_time="~5 min setup, then seconds per search",
        download_url="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
        search_mode="hmmscan",
        setup_handler="pfam_hmmscan",
        notes=(
            "Pfam HMM library — downloaded once, then scans your hits for known domains. "
            "Auto-downloads, decompresses, and runs hmmpress on first use."
        ),
        relevance=(
            "Annotation support, not primary discovery. Use to describe conserved domains in hits for papers and supplements."
        ),
    ),
    _db(
        name="VOGDB VFAM (annotation)",
        db_type="protein",
        streaming=False,
        optional=True,
        size_hint="~300-600 MB setup",
        est_time="~5-15 min setup, then seconds to minutes",
        download_url="https://fileshare.csb.univie.ac.at/vog/vog230/vfam.hmm.tar.gz",
        annotation_url="https://fileshare.csb.univie.ac.at/vog/vog230/vfam.annotations.tsv.gz",
        release="VOGDB release 230 / RefSeq release 230; 39,585 VFAMs",
        search_mode="hmmscan",
        setup_handler="vogdb_hmmscan",
        notes=(
            "Preferred viral ortholog/family annotation layer. VOGDB VFAM HMMs "
            "are downloaded once, concatenated when needed, indexed with hmmpress, "
            "and scanned with HMMER hmmscan. Optional annotation support, not a "
            "required discovery database."
        ),
        relevance=(
            "Preferred viral ortholog annotation. Use after discovery hits are found "
            "to attach VOG/VFAM family, function, and category context using stable HMMER."
        ),
    ),
    _db(
        name="PHROGs (annotation)",
        db_type="protein",
        streaming=False,
        optional=True,
        size_hint="~656 MB bundle",
        est_time="~10 min download + ~1 min hmmpress, then seconds to minutes",
        download_url="https://zenodo.org/records/17110353/files/pharokka_v1.8.0_databases.tar.gz?download=1",
        release="PHROGs v4 (via Pharokka v1.8.0 DB; 38,880 families)",
        search_mode="hmmscan",
        setup_handler="phrogs_hmmscan",
        notes=(
            "PHROGs prokaryotic-virus protein HMMs (38,880 families) from the "
            "Pharokka v1.8.0 bundle. Ships a pre-pressed all_phrogs.h3m; the hmmscan "
            "handler extracts it and re-presses the aux files. Optional annotation "
            "layer like VOGDB, not a required discovery database."
        ),
        relevance=(
            "Viral protein-family annotation. Assigns PHROG family / function / "
            "category context to discovery hits using stable HMMER hmmscan."
        ),
    ),
]

PHAGE_MODE_DEFAULTS: list[str] = [
    "INPHARED genomes",
    "RefSeq viral proteins",
    "SwissProt",
    "RefSeq viral genomes",
]

BACTERIAL_MODE_DEFAULTS: list[str] = [
    "RefSeq bacterial proteins",
    "SwissProt",
]

GENERIC_MODE_DEFAULTS: list[str] = [
    "INPHARED genomes",
    "RefSeq viral proteins",
    "SwissProt",
]
