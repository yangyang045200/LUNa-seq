#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classify raw data FASTQ inputs per sample:
- Fully on-target proviral reads (contiguous-coverage thresholds)
- Iterative contaminant/recombinant subsets (seed via 5-prime Psi & 3-prime WPRE motifs)
- Remaining sequences as non-specific residual products

NOTE:
- Streaming design: low memory, no temporary SAM on disk.
- Robust FASTA reader: tolerant to BOM / leading comments.
- Duplicate FASTQ IDs tolerated via per-ID counters.
"""

import argparse
import csv
import io
import logging
import os
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import regex
import pysam
from Bio import SeqIO

# ------------------------------------------------------------------------------
# Globals
# ------------------------------------------------------------------------------
default_ref_length_cache = {}  # cache reference lengths to avoid re-reading

# ------------------------------------------------------------------------------
# Logging / config
# ------------------------------------------------------------------------------
def setup_logging(level=logging.INFO):
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s: %(message)s")


def load_config(config_path):
    with open(config_path, "r", encoding="utf-8") as f:
        return SeqIO.yaml.safe_load(f) if hasattr(SeqIO, "yaml") else __import__("yaml").safe_load(f)

# ------------------------------------------------------------------------------
# FASTA (reference) helpers
# ------------------------------------------------------------------------------
def read_ref_len_tolerant(ref_fa_path: str) -> int:
    """
    Robustly read a single-reference FASTA length.
    Tries 'fasta', then comment-tolerant variants; finally strips leading bytes before first '>'.
    """
    for fmt in ("fasta", "fasta-pearson", "fasta-blast"):
        it = SeqIO.parse(ref_fa_path, fmt)
        rec = next(it, None)
        if rec is not None:
            return len(rec.seq)

    with open(ref_fa_path, "r", encoding="utf-8-sig", errors="ignore") as fh:
        txt = fh.read().replace("\r\n", "\n").replace("\r", "\n")
    p = txt.find(">")
    if p == -1:
        raise ValueError(f"No FASTA header ('>') found in {ref_fa_path}")
    it = SeqIO.parse(io.StringIO(txt[p:]), "fasta")
    rec = next(it, None)
    if rec is None:
        raise ValueError(f"No FASTA records found in cleaned {ref_fa_path}")
    return len(rec.seq)


def cached_ref_len(ref_fa: Path) -> int:
    key = str(ref_fa)
    L = default_ref_length_cache.get(key)
    if L is None:
        L = read_ref_len_tolerant(key)
        default_ref_length_cache[key] = L
    return L

# ------------------------------------------------------------------------------
# File discovery
# ------------------------------------------------------------------------------
def find_ref_and_move(sample_dir: Path, ref_dir: Path):
    prefix = sample_dir.name
    fasta = ref_dir / f"{prefix}.fasta"
    if fasta.exists():
        dest = sample_dir / fasta.name
        shutil.copy(fasta, dest)
        logging.info(f"Copied {fasta} to {dest}")
        return dest
    logging.warning(f"Reference {fasta} not found.")
    return None

# ------------------------------------------------------------------------------
# CIGAR processing
# ------------------------------------------------------------------------------
def longest_contiguous_len(cigartuples, indel_tol: int) -> int:
    """
    Compute the longest contiguous aligned segment length on the reference,
    allowing small I/D (<= indel_tol) to be tolerated inside the segment.
    """
    if not cigartuples:
        return 0
    longest = 0
    current = 0
    for op, length in cigartuples:
        if op == 0:  # M
            current += length
        elif op in (1, 2) and length <= indel_tol:  # small I/D tolerated
            # Only D consumes reference; we keep the original behavior that treats both as contiguous
            current += length
        elif op in (4, 5):  # S/H clipping ignored
            continue
        else:
            current = 0
        if current > longest:
            longest = current
    return longest

# ------------------------------------------------------------------------------
# Core alignment + filtering (streaming, no SAM on disk)
# ------------------------------------------------------------------------------
def align_and_filter_counts(
    reads_fq: Path,
    ref_fa: Path,
    out_fq: Path,
    min_self: float,
    min_ref: float,
    indel_tol: int,
    threads: int = 4,
):
    """
    Stream minimap2 SAM via pipe; reconstruct FASTQ for passing reads only.
    Returns:
      passed_count           : int
      total_primary          : int
      longest_passed_lengths : list[int]  (per passing read)
      passed_id_counter      : Counter[str] (qname occurrences that passed)
    """
    full_ref_len = cached_ref_len(ref_fa)
    cmd = ["minimap2", "-x", "map-ont", "-a", "-t", str(threads), str(ref_fa), str(reads_fq)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    sf = pysam.AlignmentFile(proc.stdout, "r")  # read SAM from pipe

    passed = 0
    total_primary = 0
    longest_passed_lengths = []
    passed_id_counter: Counter[str] = Counter()

    with open(out_fq, "w", encoding="utf-8") as fout:
        for aln in sf.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            qname = aln.query_name
            total_primary += 1

            L = longest_contiguous_len(aln.cigartuples or [], indel_tol)
            qlen = aln.query_length or (len(aln.query_sequence or ""))

            pass_self = (qlen > 0) and (L / qlen >= min_self)
            pass_ref = (full_ref_len > 0) and (L / full_ref_len >= min_ref)

            if pass_self or pass_ref:
                seq = aln.query_sequence or ""
                qual = aln.qual if aln.qual is not None else ("I" * len(seq))
                fout.write(f"@{qname}\n{seq}\n+\n{qual}\n")
                passed += 1
                longest_passed_lengths.append(int(L))
                passed_id_counter[qname] += 1

    sf.close()
    proc.wait()
    if proc.returncode not in (0, None):
        raise subprocess.CalledProcessError(proc.returncode, cmd)

    return passed, total_primary, longest_passed_lengths, passed_id_counter

# ------------------------------------------------------------------------------
# Remaining-writer that tolerates duplicate IDs
# ------------------------------------------------------------------------------
def write_remaining_by_counter(in_fastq: Path, exclude_counter: Counter, out_fastq: Path):
    """
    Stream input FASTQ; for each rec.id, if exclude_counter[id] > 0, decrement and skip.
    Otherwise write to out_fastq. Handles duplicate IDs correctly.
    """
    cnt = exclude_counter.copy()
    written = 0
    with open(out_fastq, "w", encoding="utf-8") as out:
        for rec in SeqIO.parse(in_fastq, "fastq"):
            rid = rec.id
            if cnt.get(rid, 0) > 0:
                cnt[rid] -= 1
                continue
            SeqIO.write(rec, out, "fastq")
            written += 1
    return written

# ------------------------------------------------------------------------------
# Contaminant seed extraction (motif-based)
# ------------------------------------------------------------------------------
def extract_contamination_ref(candidates, ltr5, ltr3, fuzzy):
    """
    From candidate reads carrying both motifs (5-prime Psi, 3-prime WPRE),
    extract a seed reference spanning from the first 5' motif to the last 3' motif.
    """
    pat5 = regex.compile(f"({ltr5}){{e<={fuzzy}}}", flags=regex.IGNORECASE)
    pat3 = regex.compile(f"({ltr3}){{e<={fuzzy}}}", flags=regex.IGNORECASE)
    best = None
    best_len = 0
    for rec in candidates:
        seq = str(rec.seq)
        matches5 = list(pat5.finditer(seq))
        matches3 = list(pat3.finditer(seq))
        if not matches5 or not matches3:
            continue
        first5 = matches5[0]
        last3 = matches3[-1]
        if first5.end() < last3.start():
            inner_len = last3.start() - first5.end()
            if inner_len > best_len:
                best = seq[first5.start(): last3.end()]
                best_len = inner_len
    return best

def contamination_stage(reads_fq: Path, ltr5, ltr3, stage: int, work_dir: Path, params: dict, indel_tol: int):
    """
    One iteration:
      - find up to max_candidates seed reads with both motifs,
      - build a seed reference (FASTA),
      - align all reads to the seed,
      - move matched reads into a stage subset (FASTQ),
      - write remaining for the next stage using a per-ID counter.
    """
    fuzzy = params.get("mismatch_max", 5)
    scan_len = params["scan_length"]
    max_cand = params["max_candidates"]

    pat5 = regex.compile(f"({ltr5}){{e<={fuzzy}}}", flags=regex.IGNORECASE)
    pat3 = regex.compile(f"({ltr3}){{e<={fuzzy}}}", flags=regex.IGNORECASE)

    # Collect a small set of candidates via a single streaming pass
    candidates = []
    for rec in SeqIO.parse(reads_fq, "fastq"):
        seq = str(rec.seq)
        head = seq[:scan_len]
        tail = seq[-scan_len:]
        if pat5.search(head) and pat3.search(tail):
            candidates.append(rec)
            if len(candidates) >= max_cand:
                break
    if not candidates:
        return None, 0, reads_fq

    ref_seq = extract_contamination_ref(candidates, ltr5, ltr3, fuzzy)
    if not ref_seq:
        return None, 0, reads_fq

    ref_fa = work_dir / f"contaminant_or_recombinant_ref_{stage}.fasta"
    with open(ref_fa, "w", encoding="utf-8") as f:
        f.write(f">contaminant_or_recombinant_ref_{stage}\n{ref_seq}\n")

    out_fq = work_dir / f"contaminant_or_recombinant_stage_{stage}.fastq"
    passed, _, _, passed_ids = align_and_filter_counts(
        reads_fq, ref_fa, out_fq, params["min_self"], params["min_ref"], indel_tol
    )
    if passed == 0:
        return None, 0, reads_fq

    rem_fq = work_dir / f"remaining_after_stage_{stage}.fastq"
    write_remaining_by_counter(reads_fq, passed_ids, rem_fq)
    return ref_fa, passed, rem_fq

# ------------------------------------------------------------------------------
# Main per-sample processing
# ------------------------------------------------------------------------------
def process_sample(sample_dir: Path, config: dict):
    work_dir = Path(sample_dir)
    script_dir = Path(__file__).resolve().parent
    ref_dir = script_dir / config["ref_folder"]
    ref_fa = find_ref_and_move(work_dir, ref_dir)
    if not ref_fa:
        return

    # Input FASTQ is assumed to be <sample>/<sample>.fastq (kept for compatibility)
    reads_fq = work_dir / f"{work_dir.name}.fastq"
    # Count total reads (streaming)
    original_total = sum(1 for _ in SeqIO.parse(reads_fq, "fastq"))

    # Primary on-target filtering
    filt = config["filter1"]
    indel_tol = int(filt.get("indel_tol", 50))

    on_target_fq = work_dir / "fully_on_target_proviral_reads.fastq"
    on_target_count, _, longest_list, passed_ids_counter = align_and_filter_counts(
        reads_fq, ref_fa, on_target_fq, filt["min_self"], filt["min_ref"], indel_tol
    )

    # "Complete coverage within on-target" using LONGEST aligned segment lengths
    ref_len = cached_ref_len(ref_fa)
    complete_cov_ratio = float(filt.get("complete_cov_ratio", filt.get("hq_self_ratio", 0.95)))
    threshold = complete_cov_ratio * ref_len
    complete_cov_count = sum(1 for L in longest_list if L >= threshold)
    complete_cov_fraction = (complete_cov_count / on_target_count) if on_target_count else 0.0

    # Remaining after on-target extraction (streaming; duplicate IDs safe)
    rem_fq = work_dir / "remaining_after_on_target.fastq"
    write_remaining_by_counter(reads_fq, passed_ids_counter, rem_fq)

    # Iterative contaminant/recombinant harvesting
    ltr_cfg = config["ltr"]
    ltr5 = ltr_cfg.get("5") or ltr_cfg.get(5)
    ltr3 = ltr_cfg.get("3") or ltr_cfg.get(3)
    max_stages = int(config.get("contamination", {}).get("max_stages", 50))

    contamination_counts = []
    stage = 1
    while stage <= max_stages:
        ref_file, passed, new_rem = contamination_stage(
            rem_fq, ltr5, ltr3, stage, work_dir, config["contamination"], indel_tol
        )
        if ref_file is None or passed == 0:
            break
        contamination_counts.append(passed)
        try:
            rem_fq.unlink()
        except OSError:
            pass
        rem_fq = new_rem
        stage += 1

    # Final remainder -> non-specific residual products (stream copy + count)
    residual_fq = work_dir / "non_specific_residual_products.fastq"
    residual_count = 0
    with open(residual_fq, "w", encoding="utf-8") as out:
        for rec in SeqIO.parse(rem_fq, "fastq"):
            SeqIO.write(rec, out, "fastq")
            residual_count += 1
    try:
        rem_fq.unlink()
    except OSError:
        pass

    # Ratios CSV (keys unchanged except they are already your latest names)
    ratios = {
        "fully_on_target_proviral_reads": (on_target_count / original_total) if original_total else 0.0,
        "complete_coverage_within_on_target": complete_cov_fraction,
    }
    for idx, count in enumerate(contamination_counts, start=1):
        ratios[f"contaminant_or_recombinant_stage_{idx}"] = (count / original_total) if original_total else 0.0
    ratios["contaminant_or_recombinant_total"] = (sum(contamination_counts) / original_total) if original_total else 0.0
    ratios["non_specific_residual_products"] = (residual_count / original_total) if original_total else 0.0

    csv_file = work_dir / f"{work_dir.name}_ratios.csv"
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["name"] + list(ratios.keys()))
        writer.writerow([work_dir.name] + [ratios[k] for k in ratios])
    logging.info(f"Processed {work_dir.name}; results in {csv_file}")

# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Filter sequences: on-target, contaminant/recombinant, residual (streaming, low-memory)."
    )
    parser.add_argument("input_folder", help="Parent folder containing sample subdirectories")
    parser.add_argument("-c", "--config", default="config.yaml", help="YAML config file path")
    args = parser.parse_args()

    setup_logging()
    # Safe YAML load
    with open(args.config, "r", encoding="utf-8") as f:
        import yaml
        config = yaml.safe_load(f)

    base = Path(args.input_folder)
    for sample in base.iterdir():
        if sample.is_dir():
            process_sample(sample, config)


if __name__ == "__main__":
    main()
