#!/usr/bin/env python3
# consensus_classification.py  (no on-disk "remaining_after_*" files)
#
# Classify consensus FASTA inputs per sample:
# - Fully on-target proviral consensus (contiguous coverage thresholds)
# - Iterative contaminant/recombinant subsets (seed via 5'Ψ & 3' WPRE motifs)
# - Remaining sequences as non-specific residual products
#
# Multiple inputs per sample are supported. Every output file name is suffixed
# with the input file's stem, e.g. "__consensus_reads_ge_5".
#
# CHANGE: intermediate "remaining_after_*" files are NOT written anymore.
# UPDATE: complete_coverage_within_on_target now uses the "longest contiguous
#         aligned segment" (per-read) as numerator criterion, NOT read length.

import os
import shutil
import argparse
import yaml
import subprocess
import logging
from pathlib import Path
from typing import Tuple, List, Dict
from tempfile import NamedTemporaryFile

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam
import csv
import regex

# Cache for reference lengths to avoid repeated FASTA reads
default_ref_length_cache: Dict[str, int] = {}


def setup_logging(level=logging.INFO):
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s: %(message)s")


def load_config(config_path: str):
    with open(config_path) as f:
        return yaml.safe_load(f)


def resolve_reference(sample_dir: Path, sample_name: str, ref_folder: str, config_dir: Path) -> Path | None:
    """
    Prefer <sample_dir>/<sample>.fasta if present.
    Else look in <config_dir>/<ref_folder>/<sample>.fasta; if found, copy into sample_dir and use it.
    """
    local = sample_dir / f"{sample_name}.fasta"
    if local.exists():
        logging.info("[%s] Using local reference: %s", sample_name, local.name)
        return local

    ref_base = (config_dir / ref_folder).resolve()
    candidate = ref_base / f"{sample_name}.fasta"
    if candidate.exists():
        dest = sample_dir / candidate.name
        shutil.copy(candidate, dest)
        logging.info("[%s] Copied reference from %s to %s", sample_name, candidate, dest)
        return dest

    logging.warning("[%s] Reference not found: %s OR %s", sample_name, local, candidate)
    return None


def _ensure_a(args_list: List[str]) -> List[str]:
    """Ensure '-a' is present for SAM output."""
    return args_list if "-a" in args_list else ["-a", *args_list]


def _parse_all_from_path(p: Path, fmt: str) -> List[SeqRecord]:
    return list(SeqIO.parse(str(p), fmt))


def _write_all_to_path(records: List[SeqRecord], out_path: Path, fmt: str):
    with open(out_path, "w") as fout:
        SeqIO.write(records, fout, fmt)


def _write_temp_fasta(records: List[SeqRecord]) -> Path:
    """Write records to a temp FASTA file and return its path."""
    tf = NamedTemporaryFile(delete=False, suffix=".fasta", mode="w")
    try:
        SeqIO.write(records, tf, "fasta")
    finally:
        tf.close()
    return Path(tf.name)


def align_and_filter_counts(
    reads_path: Path,
    ref_fa: Path,
    out_path: Path,
    min_self: float,
    min_ref: float,
    indel_tol: int,
    minimap2_args: List[str],
) -> Tuple[int, int, Dict[str, int]]:
    """
    Align sequences to reference, filter by the longest contiguous alignment segment.
    Return (passed_count, total_count, longest_passed_map).

    - longest_passed_map: {read_id -> longest_contiguous_aligned_len_on_ref} for PASSED reads.
    - Written sequences use the same format as input (FASTA here).
    - Temporary SAM is cleaned up.
    """
    # consensus inputs are FASTA
    records = _parse_all_from_path(reads_path, "fasta")
    record_dict = {rec.id: rec for rec in records}
    total = len(record_dict)

    # Cache full reference length
    ref_key = str(ref_fa)
    full_ref_len = default_ref_length_cache.get(ref_key)
    if full_ref_len is None:
        ref_record = SeqIO.read(str(ref_fa), "fasta")
        full_ref_len = len(ref_record.seq)
        default_ref_length_cache[ref_key] = full_ref_len

    # minimap2 -> SAM
    sam_path = reads_path.with_suffix('.sam')
    args_list = _ensure_a(list(minimap2_args))
    with open(sam_path, 'w') as sam_out:
        subprocess.run(
            ["minimap2", *args_list, str(ref_fa), str(reads_path)],
            check=True,
            stdout=sam_out
        )

    # Scan CIGAR to get longest contiguous segment (tolerate small indels)
    sf = pysam.AlignmentFile(str(sam_path), 'r')
    passed_ids: List[str] = []
    longest_passed: Dict[str, int] = {}
    seen = set()
    for aln in sf.fetch(until_eof=True):
        qname = aln.query_name
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped or qname in seen:
            continue
        rec = record_dict.get(qname)
        if rec is None:
            continue

        seq_len = len(rec.seq)
        longest = 0
        current = 0
        # CIGAR ops: 0=M, 1=I, 2=D, 4=S, 5=H
        for op, length in (aln.cigartuples or []):
            if op == 0:
                current += length
                longest = max(longest, current)
            elif op in (1, 2) and length <= indel_tol:
                current += length
                longest = max(longest, current)
            elif op in (4, 5):
                # clipping does not break contiguity
                continue
            else:
                # big indel / ref skip / others reset contiguity
                current = 0

        # on-target decision
        if (seq_len and longest / seq_len >= min_self) or (longest / full_ref_len >= min_ref):
            passed_ids.append(qname)
            longest_passed[qname] = int(longest)
            seen.add(qname)
    sf.close()

    # Write out (FASTA)
    out_recs = [record_dict[i] for i in passed_ids if i in record_dict]
    _write_all_to_path(out_recs, out_path, "fasta")
    passed = len(out_recs)

    # Clean up intermediate SAM
    try:
        os.remove(str(sam_path))
    except OSError:
        pass

    return passed, total, longest_passed


def extract_contamination_ref(candidates: List[SeqRecord], ltr5: str, ltr3: str, fuzzy: int) -> str | None:
    """
    From candidate sequences carrying both motifs (5' Ψ, 3' WPRE), extract a seed
    reference spanning from the first 5' motif to the last 3' motif.
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
                best = seq[first5.start():last3.end()]
                best_len = inner_len
    return best


def contamination_stage_in_memory(
    records_in: List[SeqRecord],
    ltr5: str,
    ltr3: str,
    stage: int,
    work_dir: Path,
    suffix_tag: str,
    params: dict,
    indel_tol: int,
    minimap2_args: List[str],
) -> Tuple[Path | None, int, List[SeqRecord]]:
    """
    One iteration working purely in memory:
      - seed (dual-motif) → seed reference
      - align current pool to seed, write matched subset
      - return (seed_ref_path, matched_count, new_remaining_pool)
    """
    raw = records_in
    fuzzy = params.get('mismatch_max', 3)
    scan_len = params['scan_length']
    pat5 = regex.compile(f"({ltr5}){{e<={fuzzy}}}", flags=regex.IGNORECASE)
    pat3 = regex.compile(f"({ltr3}){{e<={fuzzy}}}", flags=regex.IGNORECASE)

    cand: List[SeqRecord] = []
    for rec in raw:
        s = str(rec.seq)
        if pat5.search(s[:scan_len]) and pat3.search(s[-scan_len:]):
            cand.append(rec)
        if len(cand) >= params['max_candidates']:
            break
    if not cand:
        return None, 0, raw

    ref_seq = extract_contamination_ref(cand, ltr5, ltr3, fuzzy)
    if not ref_seq:
        return None, 0, raw

    ref_fa = work_dir / f"contaminant_or_recombinant_ref_{stage}__{suffix_tag}.fasta"
    with open(ref_fa, 'w') as f:
        f.write(f">contaminant_or_recombinant_ref_{stage}__{suffix_tag}\n{ref_seq}\n")

    tmp_reads = _write_temp_fasta(raw)
    out_subset = work_dir / f"contaminant_or_recombinant_stage_{stage}__{suffix_tag}.fasta"
    try:
        passed, _, _ = align_and_filter_counts(
            tmp_reads, ref_fa, out_subset,
            params['min_self'], params['min_ref'], indel_tol,
            minimap2_args=minimap2_args
        )
    finally:
        try:
            tmp_reads.unlink()
        except OSError:
            pass

    passed_ids = {rec.id for rec in _parse_all_from_path(out_subset, "fasta")}
    remaining = [rec for rec in raw if rec.id not in passed_ids]

    return ref_fa, passed, remaining


def process_single_input(sample_dir: Path, ref_fa: Path, input_path: Path, config: dict):
    """Run the classification workflow for one consensus FASTA file inside a sample folder."""
    work_dir = sample_dir
    stem = input_path.stem  # e.g., "consensus_reads_ge_5"
    suffix_tag = stem

    # Minimap2 arguments for consensus
    mm_args = config.get("minimap2_args_consensus", "-a -x asm5").split()

    input_records = _parse_all_from_path(input_path, "fasta")
    original_total = len(input_records)

    # Primary on-target filtering
    filt = config['filter1']
    indel_tol = filt.get('indel_tol', 50)

    tmp_in = _write_temp_fasta(input_records)
    on_target_path = work_dir / f"fully_on_target_proviral_reads__{suffix_tag}.fasta"
    try:
        on_target_count, _, longest_passed = align_and_filter_counts(
            tmp_in, ref_fa, on_target_path,
            filt['min_self'], filt['min_ref'], indel_tol,
            minimap2_args=mm_args
        )
    finally:
        try:
            tmp_in.unlink()
        except OSError:
            pass

    # Complete-coverage fraction within on-target (use "longest contiguous aligned segment")
    ref_key = str(ref_fa)
    full_ref_len = default_ref_length_cache.get(ref_key)
    if full_ref_len is None:
        ref_record = SeqIO.read(str(ref_fa), "fasta")
        full_ref_len = len(ref_record.seq)
        default_ref_length_cache[ref_key] = full_ref_len

    complete_cov_ratio = filt.get('complete_cov_ratio', filt.get('hq_self_ratio', 0.95))
    complete_cov_threshold = complete_cov_ratio * full_ref_len
    complete_cov_count = sum(1 for L in longest_passed.values() if L >= complete_cov_threshold)
    complete_cov_fraction = (complete_cov_count / on_target_count) if on_target_count else 0.0

    # Remaining pool (in memory) — still need to know which reads were written
    on_target_recs = _parse_all_from_path(on_target_path, "fasta")
    on_ids = {rec.id for rec in on_target_recs}
    remaining_pool = [rec for rec in input_records if rec.id not in on_ids]

    # Iterative contaminant/recombinant harvesting
    ltr_cfg = config['ltr']
    ltr5 = ltr_cfg.get('5') or ltr_cfg.get(5)
    ltr3 = ltr_cfg.get('3') or ltr_cfg.get(3)
    max_stages = config.get('contamination', {}).get('max_stages', 50)

    contamination_counts = []
    stage = 1
    while stage <= max_stages and remaining_pool:
        ref_file, passed, new_pool = contamination_stage_in_memory(
            remaining_pool, ltr5, ltr3, stage, work_dir, suffix_tag,
            config['contamination'], indel_tol,
            minimap2_args=mm_args
        )
        if ref_file is None or passed == 0:
            break
        contamination_counts.append(passed)
        remaining_pool = new_pool
        stage += 1

    # Final residual
    residual_path = work_dir / f"non_specific_residual_products__{suffix_tag}.fasta"
    _write_all_to_path(remaining_pool, residual_path, "fasta")

    # Ratios CSV (unchanged keys; only the internal calc of "complete_coverage_within_on_target" changed)
    ratios = {
        'fully_on_target_proviral_reads': (on_target_count / original_total) if original_total else 0.0,
        'complete_coverage_within_on_target': complete_cov_fraction,
    }
    for idx, count in enumerate(contamination_counts, start=1):
        ratios[f'contaminant_or_recombinant_stage_{idx}'] = (count / original_total) if original_total else 0.0
    ratios['contaminant_or_recombinant_total'] = (sum(contamination_counts) / original_total) if original_total else 0.0
    ratios['non_specific_residual_products'] = (len(remaining_pool) / original_total) if original_total else 0.0

    csv_file = work_dir / f"{sample_dir.name}__{suffix_tag}_ratios.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name'] + list(ratios.keys()))
        writer.writerow([f"{sample_dir.name}__{suffix_tag}"] + [ratios[k] for k in ratios])
    logging.info("Processed %s / %s; results in %s", sample_dir.name, suffix_tag, csv_file)


def process_sample(sample_dir: Path, config: dict, config_dir: Path):
    sample = sample_dir.name
    # Resolve reference (prefer local <sample>.fasta; else copy from ref_folder under config_dir)
    ref_fa = resolve_reference(sample_dir, sample, config.get('ref_folder', 'ref'), config_dir)
    if not ref_fa:
        return

    # Find consensus inputs like consensus_reads_ge_*.fasta
    pattern = config.get("consensus_glob", "consensus_reads_ge_*.fasta")
    inputs = sorted(sample_dir.glob(pattern))
    if not inputs:
        logging.info("[%s] No inputs match '%s'.", sample, pattern)
        return

    for inp in inputs:
        process_single_input(sample_dir, ref_fa, inp, config)


def main():
    parser = argparse.ArgumentParser(description="Classify consensus FASTA inputs per sample (no on-disk remaining files).")
    parser.add_argument('input_folder', help="Parent folder containing sample subdirectories (e.g., ./sample)")
    parser.add_argument('-c', '--config', default='config.yaml', help="YAML config file path")
    args = parser.parse_args()

    setup_logging()
    config = load_config(args.config)
    config_dir = Path(args.config).resolve().parent  # resolve ref_folder relative to config.yaml

    base = Path(args.input_folder).resolve()
    for sample in base.iterdir():
        if sample.is_dir():
            process_sample(sample, config, config_dir)


if __name__ == '__main__':
    main()
