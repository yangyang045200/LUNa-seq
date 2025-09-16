#!/usr/bin/env python3
# consensus_pipeline.py  (abPOA; per-UMI top-N, thresholds control which FASTQs are used)

import os
import sys
import glob
import argparse
import logging
import tempfile
import yaml
import subprocess
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)

def count_reads(fastq_path: str) -> int:
    """Return number of reads in a FASTQ file."""
    return sum(1 for _ in SeqIO.parse(fastq_path, "fastq"))

def write_top_n_fastq_with_stats(fastq_path: str, n: int, out_fastq: str):
    """
    Select the top-N reads by length (descending). If fewer than N, take all.
    Write them to `out_fastq` and return a stats dict.
    """
    recs = list(SeqIO.parse(fastq_path, "fastq"))
    recs.sort(key=lambda r: len(r.seq), reverse=True)
    picked = recs[:n]
    SeqIO.write(picked, out_fastq, "fastq")
    total_reads = len(recs)
    picked_reads = len(picked)
    mean_len_all = (sum(len(r.seq) for r in recs) / total_reads) if total_reads else 0.0
    mean_len_picked = (sum(len(r.seq) for r in picked) / picked_reads) if picked_reads else 0.0
    return {
        "total_reads": total_reads,
        "picked_reads": picked_reads,
        "mean_len_all": mean_len_all,
        "mean_len_picked": mean_len_picked,
    }

def build_abpoa_cmd(infile: str, is_fastq: bool) -> list[str]:
    """Build fixed abPOA command (tuned for Nanopore + UMI consensus)."""
    cmd = [
        "abpoa",
        "-m", "0",        # global alignment
        "-a", "0",        # heaviest-bundle consensus
        "-O", "6,26",     # gap open (convex O1,O2)
        "-E", "2,1",      # gap extend (E1,E2)
        "-L",             # sort sequences by length (desc)
        "-p",             # progressive mode
        "-S", "-k", "15", "-w", "10",  # minimizer seeding
        "-n", "500",      # minimum POA window
        "-b", "10", "-f", "0.02",      # adaptive band
        "-d", "1",        # one consensus only
        "-r", "0",        # output FASTA to stdout
    ]
    if is_fastq:
        cmd.append("-Q")  # use FASTQ quality as weights
    cmd.append(infile)
    return cmd

def run_consensus_abpoa(in_fastq: str, out_fa: str, use_quality: bool = True):
    """Run abPOA and write a single consensus in FASTA to `out_fa`."""
    cmd = build_abpoa_cmd(in_fastq, is_fastq=use_quality)
    logging.info("Running: %s", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        logging.error(res.stderr)
        raise RuntimeError("abPOA failed")
    with open(out_fa, "w") as oh:
        oh.write(res.stdout)

def main():
    parser = argparse.ArgumentParser(
        description="Per-sample abPOA consensus: filter input FASTQs by read-count thresholds, "
                    "downsample per FASTQ to top-N-by-length, then call abPOA."
    )
    parser.add_argument("-c", "--config", required=True, help="Path to YAML configuration file")
    args = parser.parse_args()
    cfg = load_config(args.config)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    base_dir = os.getcwd()
    sample_root = os.path.join(base_dir, cfg.get("sample_root", "."))

    # Glob pattern(s) for per-UMI FASTQs (e.g., <sample>/umi_fastq_*_files/*.fastq)
    raw_pats = cfg["input_pattern"]
    pats = [raw_pats] if isinstance(raw_pats, str) else raw_pats

    # abPOA input cap per FASTQ (take the longest top-N reads)
    top_n = int(cfg.get("max_reads_per_file", 20))

    # Thresholds (minimum read count per FASTQ to include that file)
    thresholds = cfg.get("thresholds", [1])
    if isinstance(thresholds, int):
        thresholds = [thresholds]
    thresholds = sorted(set(int(t) for t in thresholds))

    # Find all FASTQs
    all_fastqs = []
    for p in pats:
        glob_pat = os.path.join(sample_root, "*", p)
        logging.info("Looking for FASTQs with pattern %s", glob_pat)
        all_fastqs.extend(glob.glob(glob_pat))

    if not all_fastqs:
        logging.error("No FASTQs found under any of %s", pats)
        sys.exit(1)

    # Group by sample (parent-of-parent: <sample>/umi_fastq_*_files/<umi>.fastq)
    sample_fastqs = defaultdict(list)
    for fq in all_fastqs:
        sample = os.path.basename(os.path.dirname(os.path.dirname(fq)))
        sample_fastqs[sample].append(fq)

    # Process each sample
    for sample, fastqs in sample_fastqs.items():
        samp_dir = os.path.join(sample_root, sample)
        logging.info("=== Processing sample '%s' (%d FASTQ files) ===", sample, len(fastqs))

        # Count reads per FASTQ once (used for threshold filtering)
        read_counts = {fq: count_reads(fq) for fq in fastqs}
        for fq, cnt in read_counts.items():
            logging.info("[%s] %s → %d reads", sample, os.path.basename(fq), cnt)

        # For each threshold, include only FASTQs with reads >= threshold
        for thr in thresholds:
            passing = [fq for fq, cnt in read_counts.items() if cnt >= thr]
            out_fname = f"consensus_reads_ge_{thr}.fasta"
            out_path  = os.path.join(samp_dir, out_fname)

            if not passing:
                logging.info("[%s] No FASTQs with reads >= %d; skip %s", sample, thr, out_fname)
                continue

            logging.info("[%s] Generating %s for %d FASTQs (min reads per FASTQ ≥ %d; topN=%d)",
                         sample, out_fname, len(passing), thr, top_n)

            with open(out_path, "w") as out_handle:
                for fq in passing:
                    # (1) Downsample to longest top-N reads and collect stats
                    tmp_fastq = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False).name
                    stats = write_top_n_fastq_with_stats(fq, top_n, tmp_fastq)

                    # (2) abPOA consensus
                    tmp_cons = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False).name
                    run_consensus_abpoa(tmp_fastq, tmp_cons, use_quality=True)

                    # (3) Read consensus and set informative header (includes threshold)
                    base = os.path.splitext(os.path.basename(fq))[0]
                    tag = (
                        f"{base}"
                        f"|min_reads_ge={thr}"
                        f"|reads={stats['total_reads']}"
                        f"|picked={stats['picked_reads']}"
                        f"|mean_len_all={stats['mean_len_all']:.1f}"
                        f"|mean_len_picked={stats['mean_len_picked']:.1f}"
                        f"|topN={top_n}"
                    )
                    for rec in SeqIO.parse(tmp_cons, "fasta"):
                        rec.id = tag
                        rec.description = ""
                        SeqIO.write(rec, out_handle, "fasta")

                    # cleanup temps
                    try: os.remove(tmp_fastq)
                    except Exception: pass
                    try: os.remove(tmp_cons)
                    except Exception: pass

            logging.info("[%s] Wrote %s", sample, out_fname)

    logging.info("All samples completed.")

if __name__ == "__main__":
    main()
