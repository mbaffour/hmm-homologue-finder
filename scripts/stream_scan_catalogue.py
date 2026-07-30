#!/usr/bin/env python3
"""stream_scan_catalogue.py — six-frame scan a profile HMM through a LARGE gzipped sequence
catalogue that is streamed, chunked, and discarded as it goes.

WHY THIS EXISTS SEPARATELY FROM scan_genome_collection.sh. That script treats each line of its
list as one genome and materialises a whole batch of them; it is right for thousands of small
per-genome URLs. The metagenome catalogues are the opposite shape — ONE file of ~1.5 GB (GPD)
or ~5 GB (GVD-AVrC) holding 142k / 300k contigs, which decompress to roughly 7 GB and 25 GB.
Feeding those to a per-file batcher would put the entire catalogue on disk at once and lose the
bounded-disk property that makes this affordable.

So: decompress the HTTP stream on the fly, cut it into contig batches, scan each batch, delete
it, and keep going. Peak disk is ONE batch (tens of MB) no matter how large the catalogue, and
nothing is ever added to the database cache — which is the whole point, since the user's
constraint is to search these without downloading them.

The scan itself reuses scan_genome.py, so a hit found here has been through exactly the same
six-frame translation, read-through windowing, ORF validation and overprinting analysis as a
hit from the main pipeline. Nothing about the evidence is weaker for having been streamed.

Progress is checkpointed after every batch. A multi-hour run that dies resumes from the last
completed batch instead of starting over — the stream is in deterministic order, so skipping N
already-scanned contigs is exact.

    python3 stream_scan_catalogue.py --hmm profile.hmm --catalogue gpd --out gpd_scan
    python3 stream_scan_catalogue.py --hmm profile.hmm --url https://.../x.fa.gz --out x_scan
    ... --max-contigs 5000        # bounded smoke test
    ... --resume                  # continue an interrupted run
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent

# The two metagenome catalogues in the database catalog that a default phage run never touches.
# Kept in step with engine/databases/builtin.py.
CATALOGUES = {
    "gpd": {
        "name": "Gut Phage Database (GPD)",
        "url": "https://zenodo.org/records/6503062/files/GPD_sequences.fa.gz",
        "note": "~142k gut phage genomes, ~1.5 GB compressed",
    },
    "gvd": {
        "name": "GVD-AVrC",
        "url": "https://zenodo.org/records/11426065/files/AVrC_allrepresentatives.fasta.gz",
        "note": "~300k viral representatives, ~5 GB compressed",
    },
}

BATCH_CONTIGS = 2000        # contigs per scan; keeps a batch file in the tens of MB
BATCH_MAX_BASES = 40_000_000
UA = "hmm-homologue-finder (streaming catalogue scan)"


def fetch_compressed(url: str, dst: Path, log=print) -> Path:
    """Fetch the COMPRESSED catalogue to a temporary file, resumably, and return its path.

    WHY NOT DECOMPRESS THE SOCKET DIRECTLY. That was the first design and it fails in practice:
    scanning a batch takes ~1 minute, during which nothing reads the socket, and the server
    closes the idle connection — the stream died after 956 of 142,000 contigs with
    "Compressed file ended before the end-of-stream marker was reached". A bounded read-ahead
    queue does not fix it either, because blocking on a full queue is still an idle socket, and
    reconnecting means re-downloading gigabytes for every gap.

    So the transfer and the scan are decoupled. Peak disk becomes the COMPRESSED file (~1.5 GB
    GPD, ~5 GB GVD) plus one batch, instead of one batch alone — and it is deleted when the scan
    finishes. That is a transient working file, not a cache: nothing is added to
    ~/.cache/hmm-homologue-finder, nothing survives the run, and a later run re-fetches. The
    constraint was to avoid accumulating databases on disk, which this still honours.

    `curl -C -` resumes a partial transfer, so an interrupted download costs only the remainder.
    """
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    exe = shutil.which("curl")
    if exe:
        have = dst.stat().st_size if dst.exists() else 0
        if have:
            log(f"  resuming download from {have / 1e9:.2f} GB")
        cmd = [exe, "-sSL", "--retry", "5", "--retry-delay", "5", "-C", "-",
               "-A", UA, "-o", str(dst), url]
        subprocess.run(cmd, check=False)
    else:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=300) as r, dst.open("wb") as fh:
            shutil.copyfileobj(r, fh, 1 << 20)
    if not dst.exists() or dst.stat().st_size == 0:
        raise RuntimeError(f"download produced nothing for {url}")
    log(f"  fetched {dst.stat().st_size / 1e9:.2f} GB (temporary — deleted after the scan)")
    return dst


def _iter_contigs(fh):
    """Yield (header, sequence) from a streaming FASTA handle, one contig at a time."""
    name, chunks = None, []
    for raw in fh:
        line = raw.decode("utf-8", "replace") if isinstance(raw, bytes) else raw
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name, chunks = line[1:].strip(), []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def _scan_batch(hmm: Path, batch: Path, out: Path, min_bit: float, cpu: int, table: int) -> list:
    """Six-frame scan one batch with scan_genome.py; return its hit rows (header + data)."""
    bout = out / "_bout"
    shutil.rmtree(bout, ignore_errors=True)
    cmd = [sys.executable, str(HERE / "scan_genome.py"), "--hmm", str(hmm),
           "--genome", str(batch), "--out", str(bout), "--find-interrupted",
           "--no-neighbours", "--min-bit", str(min_bit), "--cpu", str(cpu),
           "--trans-table", str(table)]
    # Return (rows, aa, ok). `ok` is load-bearing: a batch whose scan COULD NOT RUN — missing
    # hmmsearch, an unreadable or 0-byte HMM, a crash — used to be indistinguishable from a
    # batch that legitimately contained nothing, because the output was sent to DEVNULL with
    # check=False and a missing scan_hits.tsv simply read as "no hits". The run then reported
    # status ok / hits 0 / complete true, and coverage_report turned that into "0 hits = the
    # family is absent from this catalogue". A scan that did not happen must never become a
    # biological statement.
    err = ""
    try:
        pr = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                            check=False, timeout=7200)
        rc = pr.returncode
        err = (pr.stderr or b"").decode("utf-8", "replace")[-400:]
    except subprocess.TimeoutExpired:
        shutil.rmtree(bout, ignore_errors=True)
        return [], "", False
    except Exception:
        shutil.rmtree(bout, ignore_errors=True)
        return [], "", False
    tsv = bout / "scan_hits.tsv"
    # scan_genome.py EXIT 1 MEANS "gene not detected" — a normal, expected result for most
    # batches, not a crash. Treating any non-zero exit as a failure marked 125 of 127 GVD
    # batches failed when they had simply found nothing, and reported a completed scan as
    # incomplete. Exit 1 is overloaded (argument and fetch errors also use it), so the
    # discriminator is whether the scan actually produced its table: a real "absent" result
    # still writes scan_hits.tsv, whereas a setup failure exits before the scan runs.
    ok = (rc in (0, 1)) and tsv.exists()
    rows = tsv.read_text(encoding="utf-8", errors="replace").splitlines() if tsv.exists() else []
    faa = bout / "scan_hits_aa.faa"
    aa = faa.read_text(encoding="utf-8", errors="replace") if faa.exists() else ""
    shutil.rmtree(bout, ignore_errors=True)
    if not ok and err:
        print(f"    batch scan failed: {err.strip().splitlines()[-1][:160]}")
    return rows, aa, ok


def run(hmm: Path, url: str, out: Path, label: str = "", min_bit: float = 25.0, cpu: int = 4,
        batch_contigs: int = BATCH_CONTIGS, max_contigs: int = 0, resume: bool = False,
        table: int = 11, keep_download: bool = False, log=print) -> dict:
    """Stream the catalogue, scanning and discarding batch by batch."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    agg = out / "collection_hits.tsv"
    aggaa = out / "collection_hits_aa.faa"
    ckpt = out / "_progress.json"

    done_contigs, hits, batches = 0, 0, 0
    if resume and ckpt.exists():
        try:
            st = json.loads(ckpt.read_text())
            done_contigs = int(st.get("contigs_done", 0))
            hits = int(st.get("hits", 0))
            batches = int(st.get("batches", 0))
            log(f"  resuming after {done_contigs:,} contigs already scanned ({hits} hit(s))")
        except Exception:
            done_contigs = 0
    if not resume or not agg.exists():
        agg.write_text("")
        aggaa.write_text("")

    started = time.time()
    log(f"scanning {label or url}")
    log(f"  batches of {batch_contigs} contigs, each scanned then deleted. The compressed "
        f"catalogue is a TEMPORARY working file, removed when the scan ends — nothing is added "
        f"to the database cache.")
    dl = out / ("_catalogue" + (".fa.gz" if url.endswith(".gz") else ".dat"))
    try:
        fetch_compressed(url, dl, log=log)
        gz = gzip.open(dl, "rb")
        resp = None
    except Exception as e:
        log(f"  ERROR fetching the catalogue: {e}")
        return {"status": "error", "error": str(e)[:200]}

    batch = out / "_batch.fna"
    seen = 0
    # ended_early / failed_batches are what stop a scan that DID NOT FINISH from being written
    # out as a completed one. The stream can die mid-catalogue (a truncated download, a network
    # drop) and the old code logged "re-run with --resume" and then fell through to the same
    # final write as a clean finish — status ok, complete true, exit 0 — so a partial scan was
    # published as full coverage, and at zero hits as an absence.
    ended_early, failed_batches, scanned_batches = False, 0, 0
    have_header = agg.stat().st_size > 0 if agg.exists() else False
    try:
        fh = batch.open("w", encoding="utf-8")
        n_in_batch = bases = 0
        for name, seq in _iter_contigs(gz):
            seen += 1
            if seen <= done_contigs:            # already scanned in a previous attempt
                continue
            if not seq:
                continue
            fh.write(f">{name}\n{seq}\n")
            n_in_batch += 1
            bases += len(seq)
            if n_in_batch >= batch_contigs or bases >= BATCH_MAX_BASES:
                fh.close()
                rows, aa, ok = _scan_batch(hmm, batch, out, min_bit, cpu, table)
                scanned_batches += 1
                if not ok:
                    failed_batches += 1
                if rows:
                    with agg.open("a", encoding="utf-8") as af:
                        if not have_header:
                            af.write(rows[0] + "\n")
                            have_header = True
                        for r in rows[1:]:
                            af.write(r + "\n")
                            hits += 1
                if aa:
                    with aggaa.open("a", encoding="utf-8") as bf:
                        bf.write(aa)
                batches += 1
                done_contigs = seen
                ckpt.write_text(json.dumps({"contigs_done": done_contigs, "hits": hits,
                                            "batches": batches, "url": url}, indent=2))
                el = time.time() - started
                log(f"  batch {batches}: {done_contigs:,} contigs scanned, {hits} hit(s), "
                    f"{el/60:.1f} min elapsed")
                batch.unlink(missing_ok=True)
                fh = batch.open("w", encoding="utf-8")
                n_in_batch = bases = 0
            if max_contigs and seen >= max_contigs:
                break
        fh.close()
        # the final partial batch
        if n_in_batch:
            rows, aa, ok = _scan_batch(hmm, batch, out, min_bit, cpu, table)
            scanned_batches += 1
            if not ok:
                failed_batches += 1
            if rows:
                with agg.open("a", encoding="utf-8") as af:
                    if not have_header:
                        af.write(rows[0] + "\n")
                        have_header = True
                    for r in rows[1:]:
                        af.write(r + "\n")
                        hits += 1
            if aa:
                with aggaa.open("a", encoding="utf-8") as bf:
                    bf.write(aa)
            batches += 1
            done_contigs = seen
    except KeyboardInterrupt:
        ended_early = True
        log("  interrupted — progress is checkpointed; re-run with --resume")
    except Exception as e:
        ended_early = True
        log(f"  stream ended early ({e}) — progress checkpointed; re-run with --resume")
    finally:
        try:
            gz.close()
            if resp is not None:
                resp.close()
        except Exception:
            pass
        batch.unlink(missing_ok=True)
        shutil.rmtree(out / "_bout", ignore_errors=True)
        # Delete the temporary catalogue unless the run is resumable-incomplete, in which case
        # keeping it means --resume costs no re-download. Bounded either way, and never a cache.
        if keep_download:
            log(f"  keeping {dl.name} (--keep-download)")
        elif not max_contigs:
            dl.unlink(missing_ok=True)
            log("  deleted the temporary catalogue file")

    # A resume whose checkpoint claims more contigs than the catalogue holds scans NOTHING and
    # would otherwise report a clean, complete run of that fabricated count.
    if resume and scanned_batches == 0 and seen < done_contigs:
        ended_early = True
        log(f"  ERROR: the checkpoint claims {done_contigs:,} contigs but the catalogue yielded "
            f"only {seen:,} — nothing was scanned. Wrong --out directory, or the catalogue "
            f"changed. Not recording this as a completed scan.")

    complete = bool(not max_contigs and not ended_early and failed_batches == 0)
    ckpt.write_text(json.dumps({"contigs_done": done_contigs, "hits": hits,
                                "batches": batches, "url": url, "complete": complete},
                               indent=2))
    el = (time.time() - started) / 60
    status = "ok" if complete else ("bounded" if max_contigs and not ended_early
                                    and not failed_batches else "incomplete")
    res = {"status": status, "complete": complete, "catalogue": label or url,
           "contigs_scanned": done_contigs, "batches": batches, "hits": hits,
           "minutes": round(el, 1), "hits_tsv": str(agg),
           "bounded_test": bool(max_contigs), "ended_early": ended_early,
           "failed_batches": failed_batches, "scanned_batches": scanned_batches}
    if not complete:
        res["warning"] = (
            "THIS SCAN DID NOT COVER THE WHOLE CATALOGUE"
            + (" — the stream ended early" if ended_early else "")
            + (f" — {failed_batches} batch scan(s) FAILED" if failed_batches else "")
            + (" — bounded by --max-contigs" if max_contigs else "")
            + ". A hit count from it is a LOWER BOUND and a zero is NOT evidence of absence.")
        log(f"INCOMPLETE {label or url}: {done_contigs:,} contigs, {hits} hit(s) in {el:.1f} min")
        log(f"  {res['warning']}")
        log("  re-run with --resume to continue")
    else:
        log(f"DONE {label or url}: {done_contigs:,} contigs, {hits} hit(s) in {el:.1f} min -> {agg}")
    (out / "stream_scan_summary.json").write_text(json.dumps(res, indent=2))
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hmm", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--catalogue", choices=sorted(CATALOGUES))
    g.add_argument("--url")
    ap.add_argument("--min-bit", type=float, default=25.0)
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--batch-contigs", type=int, default=BATCH_CONTIGS)
    ap.add_argument("--max-contigs", type=int, default=0,
                    help="stop after N contigs (bounded smoke test)")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--keep-download", action="store_true",
                    help="keep the temporary compressed catalogue so --resume costs no "
                         "re-download (it is deleted by default)")
    ap.add_argument("--trans-table", type=int, default=11)
    a = ap.parse_args()
    if not a.hmm.exists():
        print(f"HMM not found: {a.hmm}", file=sys.stderr)
        return 2
    url = a.url or CATALOGUES[a.catalogue]["url"]
    label = CATALOGUES[a.catalogue]["name"] if a.catalogue else a.url
    r = run(a.hmm, url, a.out, label=label, min_bit=a.min_bit, cpu=a.cpu,
            batch_contigs=a.batch_contigs, max_contigs=a.max_contigs, resume=a.resume,
            table=a.trans_table, keep_download=a.keep_download)
    # Non-zero for an incomplete scan so a wrapper cannot mistake it for full coverage.
    return 0 if r.get("status") in ("ok", "bounded") else 1


if __name__ == "__main__":
    sys.exit(main())
