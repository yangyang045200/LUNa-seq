#!/usr/bin/env python3
import os
import time
import argparse
import yaml
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from fuzzysearch import find_near_matches

def load_config():
    parser = argparse.ArgumentParser(description='Process FASTQ for UMI extraction with YAML config.')
    parser.add_argument('-c', '--config', default='config.yaml', help='YAML config file path')
    args = parser.parse_args()
    with open(args.config, 'r', encoding='utf-8') as cf:
        return yaml.safe_load(cf)

def get_parameters():
    cfg = load_config()
    umi_cfg = cfg.get('umi', {}); out_cfg = cfg.get('output', {}); inp_cfg = cfg.get('input', {})
    return {
        'prefix': umi_cfg.get('prefix', 'GCCCAATCTG'),
        'suffix': umi_cfg.get('suffix', 'CAGTGGCGCC'),
        'max_errors': umi_cfg.get('max_errors', 1),
        'pattern': umi_cfg.get('pattern', 'HBNHVNBDNHVNBDNHBD'),
        'lengths': umi_cfg.get('lengths', [18, 17, 19, 16, 20, 15, 21]),
        'target_length': umi_cfg.get('target_length', 18),
        'prefix_window': umi_cfg.get('prefix_window', 48),
        'suffix_window': umi_cfg.get('suffix_window', 48),
        'input_dir': inp_cfg.get('dir', 'sample'),
        'input_patterns': inp_cfg.get('patterns', []),
        'write_umi_fastq': out_cfg.get('write_umi_fastq', True),
        'write_filtered_fastq': out_cfg.get('write_filtered_fastq', False),
        'umi_fastq_dir': out_cfg.get('umi_fastq_dir', 'umi_fastq_files'),
        'filtered_fastq': out_cfg.get('filtered_fastq', 'filtered_seqs.fastq')
    }

PARAMS = get_parameters()
PREFIX   = PARAMS['prefix'].upper()
SUFFIX   = PARAMS['suffix'].upper()
MAX_ERR  = int(PARAMS['max_errors'])
PATTERN  = PARAMS['pattern'].upper()
LENS     = list(PARAMS['lengths'])
TARGET   = int(PARAMS['target_length'])
PRE_WIN  = int(PARAMS['prefix_window'])
SUF_WIN  = int(PARAMS['suffix_window'])
IN_DIR   = PARAMS['input_dir']
IN_PATS  = PARAMS['input_patterns']
W_UMI    = bool(PARAMS['write_umi_fastq'])
W_MERGE  = bool(PARAMS['write_filtered_fastq'])
UMI_DIR  = PARAMS['umi_fastq_dir']
MERGE_FN = PARAMS['filtered_fastq']

ALLOW = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "R": set("AG"), "Y": set("CT"), "S": set("GC"), "W": set("AT"),
    "K": set("GT"), "M": set("AC"),
    "B": set("CGT"), "D": set("AGT"), "H": set("ACT"), "V": set("ACG"),
    "N": set("ACGT")
}

def filter_umi(umi: str, pattern: str) -> bool:
    umi = umi.upper(); pattern = pattern.upper()
    if len(umi) != len(pattern): return False
    # Rule: observed UMI must be A/C/G/T only; reject any ambiguity (e.g., N)
    if any(b not in "ACGT" for b in umi): return False
    for b, p in zip(umi, pattern):
        allowed = ALLOW.get(p, {p})
        if b not in allowed: return False
    return True

def find_input_fastqs(folder):
    files = []; name = os.path.basename(folder)
    for pat in IN_PATS:
        path = os.path.join(folder, pat.format(folder=name))
        if os.path.isfile(path): files.append(path)
    return files

def process_fastq(folder, fq):
    start = time.time()
    sample = os.path.basename(folder)
    suffix = os.path.splitext(os.path.basename(fq))[0]
    print(f"Processing {sample} / {suffix}")

    counts = defaultdict(int)         # original: all lengths that hit suffix check (diagnostic)
    counts_valid = defaultdict(int)   # CHANGED: only TARGET length & filter_umi passed
    seqs = defaultdict(list)
    total_reads = 0

    for rec in SeqIO.parse(fq, 'fastq'):
        total_reads += 1
        seq = str(rec.seq).upper()

        # find prefix (allowing edit distance MAX_ERR) on forward; else try RC
        m1 = find_near_matches(PREFIX, seq[:PRE_WIN], max_l_dist=MAX_ERR)
        if not m1:
            rc = rec.reverse_complement(id=True, name=True, description=True)
            seq_rc = str(rc.seq).upper()
            m1 = find_near_matches(PREFIX, seq_rc[:PRE_WIN], max_l_dist=MAX_ERR)
            if not m1:
                continue
            rec = rc
            seq = seq_rc

        pos = m1[0].end

        # First try TARGET length for valid UMI path
        hit = False
        L = TARGET
        start_idx = pos + L
        end_idx = start_idx + len(SUFFIX)
        if end_idx <= len(seq):
            anchor = seq[start_idx:end_idx]
            errors = sum(1 for a, b in zip(anchor, SUFFIX) if a != b)
            if errors <= MAX_ERR:
                umi = seq[pos:pos+L]
                if filter_umi(umi, PATTERN):
                    counts[umi] += 1
                    counts_valid[umi] += 1  # CHANGED: record valid UMI only here
                    if W_UMI:
                        seqs[umi].append(rec)
                    hit = True

        # If TARGET not hit, try other lengths only for diagnostic distribution
        if not hit:
            for L2 in LENS:
                if L2 == TARGET:
                    continue
                start_idx = pos + L2
                end_idx = start_idx + len(SUFFIX)
                if end_idx > len(seq):
                    continue
                anchor = seq[start_idx:end_idx]
                errors = sum(1 for a, b in zip(anchor, SUFFIX) if a != b)
                if errors <= MAX_ERR:
                    umi2 = seq[pos:pos+L2]
                    counts[umi2] += 1  # diagnostic-only
                    break

    out_dir = os.path.join(folder, f"{UMI_DIR}_{suffix}_files")
    if W_UMI and seqs:
        os.makedirs(out_dir, exist_ok=True)
        for umi, recs in seqs.items():
            with open(os.path.join(out_dir, f"{umi}.fastq"), "w") as fh:
                SeqIO.write(recs, fh, "fastq")

    if W_MERGE and seqs:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, MERGE_FN), "w") as agg:
            for recs in seqs.values():
                SeqIO.write(recs, agg, "fastq")

    # --------- Summary based on counts_valid (CHANGED) ---------
    if counts_valid:
        sg = pd.Series(counts_valid).reset_index()
        sg.columns = ['UMI', 'Count']
        static_summary = sg.groupby('Count').size().reset_index(name='Frequency')
    else:
        static_summary = pd.DataFrame(columns=['Count','Frequency'])

    unique_umi = int(static_summary['Frequency'].sum()) if not static_summary.empty else 0
    filtered_reads_valid = sum(counts_valid.values())  # CHANGED
    filtered_rate = (filtered_reads_valid / total_reads) if total_reads else 0.0  # CHANGED

    thr_stats = {
        t: (static_summary.loc[static_summary['Count'] >= t, 'Frequency'].sum()
            if not static_summary.empty else 0)
        for t in [3, 5, 10, 20]
    }

    singleton_freq = static_summary.loc[static_summary['Count'] == 1, 'Frequency']
    singleton_rate = (float(singleton_freq.iloc[0]) / unique_umi) if (not singleton_freq.empty and unique_umi) else 0.0

    # CHANGED: duplication_rate uses filtered_reads_valid (UMI-level duplicate extent)
    duplication_rate = 1 - (unique_umi / filtered_reads_valid) if filtered_reads_valid else 0.0

    # Redundant fraction also on valid path (kept consistent with new base)
    redundant_reads = (
        (static_summary.loc[static_summary['Count'] > 1, 'Count']
         * static_summary.loc[static_summary['Count'] > 1, 'Frequency']).sum()
        if not static_summary.empty else 0
    )
    redundant_frac = float(redundant_reads / filtered_reads_valid) if filtered_reads_valid else 0.0

    summary = {
        'Demultiplexed':            [total_reads],
        'Unique_UMI':               [unique_umi],
        'Filtered_rate(%)':         [round(filtered_rate * 100, 2)],
        'Singleton_rate(%)':        [round(singleton_rate * 100, 2)],
        'Duplication_rate(%)':      [round(duplication_rate * 100, 2)],   # CHANGED
        'Redundant_frac(%)':        [round(redundant_frac * 100, 2)],
    }
    for t, cnt in thr_stats.items():
        summary[f'UMI>={t}_count']   = [int(cnt)]
        summary[f'UMI>={t}_rate(%)'] = [round((cnt / unique_umi * 100), 2) if unique_umi else 0]

    df1 = pd.DataFrame(summary)

    # --------- LengthStats (diagnostic) still from counts (unchanged idea) ---------
    rows = []
    total_len_reads = sum(counts.values())
    for L in sorted(set(LENS + [TARGET])):
        reads = sum(c for u, c in counts.items() if len(u) == L)
        # kinds: number of distinct UMI species of this length
        kinds = len({u for u in counts if len(u) == L})
        ratio = f"{round(reads / total_len_reads * 100, 2)}%" if total_len_reads else '0%'
        rows.append([L, kinds, reads, ratio])
    df2 = pd.DataFrame(rows, columns=['Length', 'UMI_Count', 'Reads', 'Ratio'])

    # Outputs
    static_summary.to_csv(os.path.join(folder, f"static_{suffix}.csv"), index=False)
    excel_path = os.path.join(folder, f"result_record_{suffix}.xlsx")
    with pd.ExcelWriter(excel_path) as writer:
        df1.to_excel(writer, sheet_name='Summary', index=False)
        df2.to_excel(writer, sheet_name='LengthStats', index=False)

    print(f"{sample}/{suffix}: done in {time.time()-start:.2f}s")

def main():
    for d in os.listdir(IN_DIR):
        folder = os.path.join(IN_DIR, d)
        if os.path.isdir(folder):
            for fq in find_input_fastqs(folder):
                process_fastq(folder, fq)

if __name__ == '__main__':
    main()
