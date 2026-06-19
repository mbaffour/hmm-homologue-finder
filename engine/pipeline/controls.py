"""
pipeline/controls.py — Built-in positive and negative controls for HMM validation.

Controls are used in Step 5 (Validation & Calibration) to measure:
  - Sensitivity  : how well the HMM recovers known positives
  - Specificity  : how few false positives appear in negative sets
  - ROC-like data: score distributions for both sets

Each biology mode has domain-appropriate control sets:
  - generic  : seeds (pos) + shuffled sequences (neg)
  - phage    : seeds (pos) + fungal / mammalian proteins (neg)
  - bacterial: seeds (pos) + unrelated viral proteins (neg)

Control FASTA files are stored in www/controls/ and ship with the app.
Users can register their own controls via the UI.
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

from .utils import find_tool, run_cmd

# ---------------------------------------------------------------------------
# Built-in control catalogue
# ---------------------------------------------------------------------------

# Absolute path to the controls data directory (www/controls/)
_APP_DIR  = Path(__file__).parent.parent
_CTRL_DIR = _APP_DIR / "www" / "controls"

BUILTIN_CONTROLS: list[dict] = [
    # ── Positive controls ──────────────────────────────────────────────
    {
        "name":        "seeds_self_test",
        "role":        "positive",
        "desc":        "Seed sequences — 100% recovery expected",
        "modes":       ["generic", "phage", "bacterial"],
        "file":        None,            # set dynamically from seed_faa
        "source":      "user_seed",
        "expect_hits": True,
        "n_seqs":      None,
    },
    # ── Negative controls — phage mode ─────────────────────────────────
    {
        "name":    "fungi_proteome_sample",
        "role":    "negative",
        "desc":    "500 fungal proteins (Candida, Aspergillus, Saccharomyces)",
        "modes":   ["phage", "generic"],
        "file":    "fungi_500.faa",
        "source":  "uniprot_fungi",
        "expect_hits": False,
        "n_seqs":  500,
    },
    {
        "name":    "mammalian_housekeeping",
        "role":    "negative",
        "desc":    "200 human/mouse housekeeping proteins (ribosomal, metabolic)",
        "modes":   ["phage"],
        "file":    "mammalian_200.faa",
        "source":  "uniprot_human",
        "expect_hits": False,
        "n_seqs":  200,
    },
    {
        "name":    "archaeal_proteins",
        "role":    "negative",
        "desc":    "300 archaeal proteins (Methanobacterium, Haloarcula)",
        "modes":   ["phage", "generic"],
        "file":    "archaea_300.faa",
        "source":  "uniprot_archaea",
        "expect_hits": False,
        "n_seqs":  300,
    },
    # ── Negative controls — bacterial mode ─────────────────────────────
    {
        "name":    "viral_structural_proteins",
        "role":    "negative",
        "desc":    "200 eukaryotic virus structural proteins (HIV, influenza)",
        "modes":   ["bacterial"],
        "file":    "euk_viral_200.faa",
        "source":  "uniprot_viral",
        "expect_hits": False,
        "n_seqs":  200,
    },
    {
        "name":    "plant_proteins",
        "role":    "negative",
        "desc":    "400 plant proteins (Arabidopsis thaliana)",
        "modes":   ["bacterial", "generic"],
        "file":    "plant_400.faa",
        "source":  "uniprot_plant",
        "expect_hits": False,
        "n_seqs":  400,
    },
    # ── Universal negative control ──────────────────────────────────────
    {
        "name":    "shuffled_seeds",
        "role":    "negative",
        "desc":    "Amino-acid-shuffled seeds — same composition, random order",
        "modes":   ["generic", "phage", "bacterial"],
        "file":    None,        # generated on the fly from seed_faa
        "source":  "shuffled",
        "expect_hits": False,
        "n_seqs":  None,
    },
]


# ---------------------------------------------------------------------------
# Control file management
# ---------------------------------------------------------------------------

def available_controls(mode: str = "generic", app_dir: Optional[Path] = None) -> list[dict]:
    """Return controls available for the given biology mode.

    A control is "available" when either:
    - its ``file`` is present in ``www/controls/``, or
    - it is generated dynamically (file=None).

    Parameters
    ----------
    mode : str
        One of ``"generic"``, ``"phage"``, ``"bacterial"``.
    app_dir : Path, optional
        Root of the app (default: this file's parent.parent).

    Returns
    -------
    list[dict]
        Control descriptors with an added ``path`` key (Path | None).
    """
    base = (app_dir or _APP_DIR) / "www" / "controls"
    result = []
    for ctrl in BUILTIN_CONTROLS:
        if mode not in ctrl["modes"]:
            continue
        entry = dict(ctrl)
        if ctrl["file"]:
            p = base / ctrl["file"]
            entry["path"] = p if p.exists() else None
            entry["available"] = p.exists()
        else:
            # Dynamic controls (seeds_self_test, shuffled_seeds) are always
            # "available" once the seed FASTA is known
            entry["path"] = None
            entry["available"] = True  # will be created at run time
        result.append(entry)
    return result


def generate_shuffled_control(
    seed_faa: Path,
    out_path: Optional[Path] = None,
    n: int = 100,
    seed_rng: int = 42,
) -> Path:
    """Generate a shuffled-sequence negative control from the seed FASTA.

    Each sequence has its amino acids randomly permuted so composition is
    identical to the real sequences but primary structure is destroyed.

    Parameters
    ----------
    seed_faa : Path
        Input seed FASTA.
    out_path : Path, optional
        Where to write the shuffled FASTA (default: sibling ``shuffled_ctrl.faa``).
    n : int
        Maximum number of shuffled sequences to generate.
    seed_rng : int
        Random seed for reproducibility.

    Returns
    -------
    Path
        Path to the generated FASTA.
    """
    seed_faa = Path(seed_faa)
    if out_path is None:
        out_path = seed_faa.parent / "shuffled_ctrl.faa"
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
    except ImportError:
        print("ERROR: Biopython required for shuffled control generation.", file=sys.stderr)
        return Path()

    random.seed(seed_rng)
    records = list(SeqIO.parse(str(seed_faa), "fasta"))
    if not records:
        print(f"ERROR: No sequences found in {seed_faa}", file=sys.stderr)
        return Path()

    random.shuffle(records)
    records = records[:n]

    shuffled = []
    for i, rec in enumerate(records):
        aa_list = list(str(rec.seq).upper())
        random.shuffle(aa_list)
        shuffled.append(
            SeqRecord(
                Seq("".join(aa_list)),
                id=f"shuffled_{i+1:04d}",
                description=f"shuffled|original={rec.id}",
            )
        )

    SeqIO.write(shuffled, str(out_path), "fasta")
    return out_path


# ---------------------------------------------------------------------------
# Control runner
# ---------------------------------------------------------------------------

def run_control_search(
    hmm_path: Path,
    control_faa: Path,
    out_dir: Path,
    control_name: str,
    evalue: float = 0.001,
    cpu: int = 4,
) -> dict:
    """Run hmmsearch against a control FASTA and return metrics.

    Parameters
    ----------
    hmm_path : Path
        Profile HMM file.
    control_faa : Path
        Control sequence FASTA.
    out_dir : Path
        Directory for hmmsearch output.
    control_name : str
        Label used in output file names.
    evalue : float
        E-value threshold (use a lenient value to see score distributions).
    cpu : int
        Threads for hmmsearch.

    Returns
    -------
    dict with keys:
        control_name, role, n_seqs, n_hits, hit_rate,
        score_distribution (list of floats), fp_rate,
        sensitivity (only for positive controls),
        result_file
    """
    hmm_path    = Path(hmm_path)
    control_faa = Path(control_faa)
    out_dir     = Path(out_dir)

    if not hmm_path.exists():
        return {"error": f"HMM not found: {hmm_path}"}
    if not control_faa.exists():
        return {"error": f"Control FASTA not found: {control_faa}"}

    out_dir.mkdir(parents=True, exist_ok=True)
    safe_name = control_name.replace(" ", "_").lower()
    tbl_out   = out_dir / f"ctrl_{safe_name}.tbl"

    hmmsearch_bin = find_tool("hmmsearch") or "hmmsearch"
    cmd = [
        hmmsearch_bin,
        "--tblout", str(tbl_out),
        "--noali",
        "-E", str(evalue),
        "--cpu", str(cpu),
        str(hmm_path),
        str(control_faa),
    ]

    result = run_cmd(cmd)
    if result.returncode != 0:
        return {"error": f"hmmsearch failed: {result.stderr[-300:]}"}

    # Parse tblout
    n_seqs, hits, scores = _parse_tbl_simple(tbl_out)
    if n_seqs == 0:
        # Count from FASTA
        try:
            from Bio import SeqIO
            n_seqs = sum(1 for _ in SeqIO.parse(str(control_faa), "fasta"))
        except Exception:
            n_seqs = 1  # avoid div-by-zero

    hit_rate = len(hits) / n_seqs if n_seqs > 0 else 0.0

    return {
        "control_name":        control_name,
        "n_seqs":              n_seqs,
        "n_hits":              len(hits),
        "hit_rate":            round(hit_rate, 4),
        "hit_rate_pct":        round(hit_rate * 100, 2),
        "fp_rate":             round(hit_rate, 4),   # alias for negatives
        "scores":              scores,
        "min_score":           min(scores) if scores else 0.0,
        "max_score":           max(scores) if scores else 0.0,
        "mean_score":          round(sum(scores) / len(scores), 2) if scores else 0.0,
        "result_file":         str(tbl_out),
        "hit_ids":             hits,
    }


def run_all_controls(
    hmm_path: Path,
    seed_faa: Path,
    out_dir: Path,
    mode: str = "generic",
    strict_threshold: float = 45.0,
    moderate_threshold: float = 30.0,
    cpu: int = 4,
    app_dir: Optional[Path] = None,
) -> "ControlReport":
    """Run all available controls and return a :class:`ControlReport`.

    Parameters
    ----------
    hmm_path : Path
        Trained profile HMM.
    seed_faa : Path
        Seed FASTA used to build the HMM.
    out_dir : Path
        Directory for all control output.
    mode : str
        Biology mode (``"generic"``, ``"phage"``, ``"bacterial"``).
    strict_threshold : float
        Strict bit-score threshold.
    moderate_threshold : float
        Moderate bit-score threshold.
    cpu : int
        Threads.
    app_dir : Path, optional
        App root for locating control files.

    Returns
    -------
    ControlReport
    """
    out_dir  = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seed_faa = Path(seed_faa)
    controls = available_controls(mode=mode, app_dir=app_dir)

    results: list[dict] = []

    for ctrl in controls:
        print(f"INFO: Running control: {ctrl['name']} ({ctrl['role']})", file=sys.stderr)

        # Resolve FASTA path
        if ctrl["source"] == "user_seed":
            faa = seed_faa
        elif ctrl["source"] == "shuffled":
            shuffled_path = out_dir / "shuffled_ctrl.faa"
            faa = generate_shuffled_control(seed_faa, shuffled_path)
            if not faa.exists():
                print(f"WARNING: Could not generate shuffled control.", file=sys.stderr)
                continue
        elif ctrl["path"] and ctrl["available"]:
            faa = ctrl["path"]
        else:
            print(
                f"INFO: Control '{ctrl['name']}' FASTA not found — skipping. "
                f"(Run setup_controls.py to download control sequences.)",
                file=sys.stderr,
            )
            continue

        res = run_control_search(
            hmm_path     = hmm_path,
            control_faa  = faa,
            out_dir      = out_dir,
            control_name = ctrl["name"],
            evalue       = 1.0,  # very lenient for score distribution
            cpu          = cpu,
        )
        res["role"]         = ctrl["role"]
        res["desc"]         = ctrl["desc"]
        res["expect_hits"]  = ctrl["expect_hits"]
        res["strict_hits"]  = sum(1 for s in res.get("scores", []) if s >= strict_threshold)
        res["moderate_hits"]= sum(1 for s in res.get("scores", []) if s >= moderate_threshold)
        results.append(res)

    return ControlReport(results, strict_threshold, moderate_threshold)


# ---------------------------------------------------------------------------
# ControlReport — aggregated metrics
# ---------------------------------------------------------------------------

class ControlReport:
    """Aggregated control results with sensitivity/specificity metrics."""

    def __init__(
        self,
        results: list[dict],
        strict_threshold: float = 45.0,
        moderate_threshold: float = 30.0,
    ):
        self.results           = results
        self.strict_threshold  = strict_threshold
        self.moderate_threshold= moderate_threshold

    # ---- Per-role aggregation -----------------------------------------------

    def positive_results(self) -> list[dict]:
        return [r for r in self.results if r.get("role") == "positive"]

    def negative_results(self) -> list[dict]:
        return [r for r in self.results if r.get("role") == "negative"]

    def sensitivity(self, threshold: Optional[float] = None) -> float:
        """Fraction of positive control sequences recovered above threshold."""
        t = threshold if threshold is not None else self.strict_threshold
        total = sum(r.get("n_seqs", 0) for r in self.positive_results())
        hits  = sum(
            sum(1 for s in r.get("scores", []) if s >= t)
            for r in self.positive_results()
        )
        return round(hits / total, 4) if total > 0 else 0.0

    def specificity(self, threshold: Optional[float] = None) -> float:
        """Fraction of negative control sequences correctly rejected (1 - FPR)."""
        t = threshold if threshold is not None else self.strict_threshold
        total = sum(r.get("n_seqs", 0) for r in self.negative_results())
        fps   = sum(
            sum(1 for s in r.get("scores", []) if s >= t)
            for r in self.negative_results()
        )
        tn = total - fps
        return round(tn / total, 4) if total > 0 else 1.0

    def false_positive_rate(self, threshold: Optional[float] = None) -> float:
        return round(1.0 - self.specificity(threshold), 4)

    def summary(self) -> dict:
        """Return a compact summary dict for the UI.

        Keys follow the names used by step_05_validate.py:
          sensitivity, specificity, false_positive_rate,
          total_positives, true_positives, total_negatives, false_positives
        """
        t = self.strict_threshold
        total_pos = sum(r.get("n_seqs", 0) for r in self.positive_results())
        tp        = sum(
            sum(1 for s in r.get("scores", []) if s >= t)
            for r in self.positive_results()
        )
        total_neg = sum(r.get("n_seqs", 0) for r in self.negative_results())
        fp        = sum(
            sum(1 for s in r.get("scores", []) if s >= t)
            for r in self.negative_results()
        )
        sens = round(tp / total_pos, 4) if total_pos > 0 else 0.0
        spec = round((total_neg - fp) / total_neg, 4) if total_neg > 0 else 1.0
        fpr  = round(fp / total_neg, 4) if total_neg > 0 else 0.0
        return {
            # Short canonical keys used by UI
            "sensitivity":         sens,
            "specificity":         spec,
            "false_positive_rate": fpr,
            "total_positives":     total_pos,
            "true_positives":      tp,
            "total_negatives":     total_neg,
            "false_positives":     fp,
            # Verbose aliases kept for backwards compat
            "n_controls":           len(self.results),
            "n_positive_controls":  len(self.positive_results()),
            "n_negative_controls":  len(self.negative_results()),
            "sensitivity_strict":   sens,
            "sensitivity_moderate": self.sensitivity(self.moderate_threshold),
            "specificity_strict":   spec,
            "results":              self.results,
        }

    def to_dataframe(self) -> pd.DataFrame:
        """Convert results to a tidy DataFrame for display.

        Columns match the names used by step_05_validate.py:
          name, role, n_sequences, n_hits_strict, min_score, max_score
        """
        rows = []
        for r in self.results:
            rows.append({
                "name":         r.get("control_name", r.get("name", "")),
                "role":         r.get("role", ""),
                "desc":         r.get("desc", ""),
                "n_sequences":  r.get("n_seqs", 0),
                "n_hits":       r.get("n_hits", 0),
                "n_hits_strict":r.get("strict_hits", r.get("n_hits", 0)),
                "hit_rate_pct": r.get("hit_rate_pct", 0.0),
                "min_score":    r.get("min_score", 0.0),
                "max_score":    r.get("max_score", 0.0),
                "mean_score":   r.get("mean_score", 0.0),
                "pass": (
                    r.get("n_hits", 0) > 0 if r.get("role") == "positive"
                    else r.get("n_hits", 0) == 0
                ),
            })
        return pd.DataFrame(rows)

    def to_json(self) -> str:
        return json.dumps(self.summary(), indent=2, default=str)


# ---------------------------------------------------------------------------
# Control sequence downloader / bundler
# ---------------------------------------------------------------------------

def download_control_sequences(
    out_dir: Optional[Path] = None,
    n_per_taxon: int = 100,
) -> dict:
    """
    Download small curated control sequence sets from UniProt REST API.

    Fetches reviewed (Swiss-Prot) entries for:
      - Fungi         (tax:4751)  → fungi_500.faa
      - Mammalia      (tax:40674) → mammalian_200.faa
      - Archaea       (tax:2157)  → archaea_300.faa
      - Eukaryotic viruses (tax:10239 NOT tax:10244) → euk_viral_200.faa
      - Embryophyta   (tax:3193)  → plant_400.faa

    Parameters
    ----------
    out_dir : Path, optional
        Where to write FASTAs (default: www/controls/).
    n_per_taxon : int
        Number of sequences to fetch per taxonomic group.

    Returns
    -------
    dict
        {filename: n_sequences} for each group downloaded.
    """
    dest = (out_dir or _CTRL_DIR)
    dest.mkdir(parents=True, exist_ok=True)

    import urllib.request
    import time

    UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb/stream"

    groups = [
        ("fungi_500.faa",      "taxonomy_id:4751 AND reviewed:true",           500),
        ("mammalian_200.faa",  "taxonomy_id:40674 AND reviewed:true",          200),
        ("archaea_300.faa",    "taxonomy_id:2157 AND reviewed:true",           300),
        ("euk_viral_200.faa",  "taxonomy_id:10239 AND reviewed:true",          200),
        ("plant_400.faa",      "taxonomy_id:3193 AND reviewed:true",           400),
    ]

    downloaded = {}
    for fname, query, n in groups:
        out_faa = dest / fname
        if out_faa.exists() and out_faa.stat().st_size > 1000:
            print(f"  {fname}: already exists, skipping.", file=sys.stderr)
            downloaded[fname] = sum(1 for line in out_faa.open() if line.startswith(">"))
            continue

        try:
            import urllib.parse
            url = (
                f"{UNIPROT_BASE}?query={urllib.parse.quote(query)}"
                f"&format=fasta&size={n}"
            )
            print(f"  Downloading {fname} ...", file=sys.stderr)
            with urllib.request.urlopen(url, timeout=60) as resp:
                content = resp.read().decode("utf-8", errors="replace")
            out_faa.write_text(content)
            n_seqs = content.count(">")
            downloaded[fname] = n_seqs
            print(f"  {fname}: {n_seqs} sequences", file=sys.stderr)
            time.sleep(1)  # be polite to UniProt
        except Exception as exc:
            print(f"  WARNING: Could not download {fname}: {exc}", file=sys.stderr)
            downloaded[fname] = 0

    return downloaded


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_tbl_simple(tbl_path: Path) -> tuple[int, list[str], list[float]]:
    """Parse hmmsearch tblout; return (n_searched, hit_ids, bit_scores).

    Note: tblout doesn't tell us how many sequences were searched; n_searched
    is determined by counting FASTA sequences separately.
    """
    hit_ids: list[str] = []
    scores:  list[float] = []

    if not tbl_path.exists():
        return 0, hit_ids, scores

    with open(tbl_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                hit_ids.append(parts[0])
                scores.append(float(parts[5]))  # full-sequence bit score
            except (IndexError, ValueError):
                continue

    return 0, hit_ids, scores
