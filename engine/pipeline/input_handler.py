"""
pipeline/input_handler.py — Parse input sequences; auto-detect type.

Handles FASTA (.fasta, .fa, .faa, .fna), GenBank (.gb, .gbk), gzip-compressed
variants, and folders of mixed files. Auto-detects protein vs nucleotide.
"""
import gzip
import sys
from pathlib import Path
from typing import Union

# ---------------------------------------------------------------------------
# Optional Biopython — fail gracefully at import time rather than at call time
# ---------------------------------------------------------------------------
try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    _BIO_AVAILABLE = True
except ImportError:
    _BIO_AVAILABLE = False
    SeqRecord = object  # type: ignore

# Characters that are exclusively nucleotide (IUPAC + gap)
_NT_CHARS = set("ACGTNUWSMKRYBDHVacgtnuwsmkrybdhv-.")

# Characters that ONLY appear in protein sequences (not in any IUPAC nt code).
# IUPAC nt covers: A C G T U W S M K R Y B D H V N — so amino acids
# exclusively NOT in nt are: E F I L O P Q Z J X (standard + extended).
# Using a low threshold (≥5% protein-only chars) is sufficient to call protein
# because real protein sequences are full of L, I, F, P, Q, E.
_PROTEIN_ONLY_CHARS = set("EFILPQZJXefilpqzjx")

# Extensions recognised as sequence files
_FASTA_EXTS  = {".fasta", ".fa", ".faa", ".fna"}
_GENBANK_EXTS = {".gb", ".gbk"}
_ALL_EXTS    = _FASTA_EXTS | _GENBANK_EXTS


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def detect_file_type(path: Path) -> str:
    """Return "fasta", "genbank", or "folder" based on path suffix / type.

    Parameters
    ----------
    path : Path
        File or directory to inspect.

    Returns
    -------
    str
        "fasta", "genbank", or "folder".
    """
    path = Path(path)
    if path.is_dir():
        return "folder"
    suffix = path.suffix.lower()
    # Strip .gz to inspect the inner suffix
    if suffix == ".gz":
        suffix = Path(path.stem).suffix.lower()
    if suffix in _GENBANK_EXTS:
        return "genbank"
    return "fasta"  # Default — covers .fa, .faa, .fna, .fasta, unknown


def detect_seq_type(path: Path) -> str:
    """Return "protein" or "nucleotide" by inspecting the first 50 sequences.

    Rule: if >50 % of characters in sampled sequences are outside the
    IUPAC nucleotide alphabet → "protein".

    Parameters
    ----------
    path : Path
        Path to a sequence file (FASTA or GenBank).

    Returns
    -------
    str
        "protein" or "nucleotide".
    """
    path = Path(path)

    # Fast path: extension is unambiguous
    suffix = path.suffix.lower()
    if suffix == ".gz":
        suffix = Path(path.stem).suffix.lower()
    if suffix == ".faa":
        return "protein"
    if suffix == ".fna":
        return "nucleotide"

    # Character-inspection: count chars that can ONLY appear in proteins
    # (E, F, I, L, P, Q, Z, J — not present in any IUPAC nt code).
    # A genuine protein sequence typically has ≥5% such chars.
    records = load_sequences(path)[:50]
    if not records:
        return "nucleotide"

    total = 0
    protein_only = 0
    for rec in records:
        seq_str = str(rec.seq).upper()
        for ch in seq_str:
            if ch in ("-", ".", "*"):
                continue
            total += 1
            if ch in _PROTEIN_ONLY_CHARS:
                protein_only += 1

    if total == 0:
        return "nucleotide"
    return "protein" if (protein_only / total) > 0.05 else "nucleotide"


def load_sequences(path: Path) -> "list[SeqRecord]":
    """Load sequences from a FASTA or GenBank file (gzip transparent).

    For GenBank files, protein translations are extracted from CDS features
    when available; if none are annotated the raw CDS nucleotide sequence is
    returned instead.

    Parameters
    ----------
    path : Path
        File to parse (.fasta/.fa/.faa/.fna/.gb/.gbk, optionally .gz).

    Returns
    -------
    list[SeqRecord]
        Parsed records, or empty list on any error.
    """
    if not _BIO_AVAILABLE:
        print("ERROR: Biopython not installed — cannot parse sequences.", file=sys.stderr)
        return []

    path = Path(path)
    if not path.exists():
        print(f"ERROR: File not found: {path}", file=sys.stderr)
        return []

    file_type = detect_file_type(path)
    if file_type == "folder":
        records: list = []
        for fp in find_files_in_folder(path):
            records.extend(load_sequences(fp))
        return records

    # Determine format string for SeqIO
    fmt = "genbank" if file_type == "genbank" else "fasta"

    try:
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        mode = "rt"
        with opener(path, mode) as fh:  # type: ignore[call-overload]
            if fmt == "genbank":
                return _load_genbank_records(fh)
            else:
                return list(SeqIO.parse(fh, "fasta"))
    except Exception as exc:
        print(f"ERROR: Could not parse {path}: {exc}", file=sys.stderr)
        return []


def input_summary(path: Path) -> dict:
    """Return a summary dict for the sequences at *path*.

    Returns
    -------
    dict with keys:
        seq_count, avg_len, min_len, max_len, seq_type, duplicates, alphabet
    """
    path = Path(path)
    records = load_sequences(path)

    if not records:
        return {
            "seq_count": 0,
            "avg_len": 0,
            "min_len": 0,
            "max_len": 0,
            "seq_type": "unknown",
            "duplicates": 0,
            "alphabet": "unknown",
        }

    lengths = [len(r.seq) for r in records]
    ids = [r.id for r in records]
    seq_type = detect_seq_type(path)

    # Detect alphabet characters actually present
    all_chars: set = set()
    for rec in records:
        all_chars.update(str(rec.seq).upper())
    alphabet = "".join(sorted(all_chars))

    return {
        "seq_count": len(records),
        "avg_len": round(sum(lengths) / len(lengths), 1),
        "min_len": min(lengths),
        "max_len": max(lengths),
        "seq_type": seq_type,
        "duplicates": len(ids) - len(set(ids)),
        "alphabet": alphabet[:80],  # truncate for display
    }


def find_files_in_folder(path: Path) -> "list[Path]":
    """Recursively find all recognised sequence files under *path*.

    Parameters
    ----------
    path : Path
        Directory to search.

    Returns
    -------
    list[Path]
        Sorted list of matching file paths.
    """
    path = Path(path)
    if not path.is_dir():
        print(f"WARNING: {path} is not a directory.", file=sys.stderr)
        return []

    found: list = []
    for ext in _ALL_EXTS:
        found.extend(path.rglob(f"*{ext}"))
        found.extend(path.rglob(f"*{ext}.gz"))

    return sorted(set(found))


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_genbank_records(fh) -> "list[SeqRecord]":
    """Parse GenBank file handle; prefer protein translations from CDS features."""
    if not _BIO_AVAILABLE:
        return []

    records: list = []
    for gb_rec in SeqIO.parse(fh, "genbank"):
        has_translation = False
        for feature in gb_rec.features:
            if feature.type != "CDS":
                continue
            if "translation" in feature.qualifiers:
                aa_seq = feature.qualifiers["translation"][0]
                gene_name = (
                    feature.qualifiers.get("gene", [""])
                    or feature.qualifiers.get("product", [""])
                )[0]
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                prot_id = feature.qualifiers.get("protein_id", [locus_tag or gb_rec.id])[0]
                from Bio.SeqRecord import SeqRecord as SR
                from Bio.Seq import Seq
                prot_rec = SR(
                    Seq(aa_seq),
                    id=prot_id,
                    description=f"{gene_name} [{gb_rec.description}]",
                )
                records.append(prot_rec)
                has_translation = True

        if not has_translation:
            # Fall back to the raw nucleotide record
            records.append(gb_rec)

    return records
