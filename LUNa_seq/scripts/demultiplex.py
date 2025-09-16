#!/usr/bin/env python3
import argparse
import logging
import os
import gzip
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import yaml
from Bio import SeqIO


# ---------- Config & Logging ----------

def load_config(path: str) -> dict:
    """
    Load YAML configuration. Supports two layouts:
      1) Top-level keys (fastq_path, demux_table, etc.)
      2) Nested under a 'demultiplex' key.
    Required keys: fastq_path, demux_table, output_dir, seqkit_cmd, search_range
    """
    with open(path, 'r') as f:
        cfg_full = yaml.safe_load(f)
    if not cfg_full:
        raise ValueError(f"Configuration file {path} is empty or invalid YAML.")

    cfg = cfg_full.get('demultiplex', cfg_full)

    required = ['fastq_path', 'demux_table', 'output_dir', 'seqkit_cmd', 'search_range']
    missing = [k for k in required if k not in cfg]
    if missing:
        raise KeyError(f"Missing required config keys: {', '.join(missing)}")

    return cfg


def setup_logging(level=logging.INFO):
    logging.basicConfig(level=level, format="%(asctime)s [%(levelname)s] %(message)s")


# ---------- Small IO helpers ----------

def open_fastq(path: Path, mode='rt'):
    """Open FASTQ file that may be gzipped."""
    if str(path).endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def count_total_reads_fastq(seqkit_cmd: str, filename: Path) -> int:
    """Count reads via `seqkit stats`."""
    try:
        result = subprocess.run(
            [seqkit_cmd, 'stats', str(filename)],
            capture_output=True, text=True, check=True
        )
        for line in result.stdout.splitlines():
            if not line.startswith('file'):
                cols = line.split()
                return int(cols[3].replace(',', ''))
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error running seqkit stats: {e}")
    return 0


def file_line_count(path: Path) -> int:
    """Count lines of a small text file (e.g., ID list)."""
    if not path.exists() or path.stat().st_size == 0:
        return 0
    out = subprocess.run(['wc', '-l', str(path)], capture_output=True, text=True, check=True)
    return int(out.stdout.strip().split()[0])


def safe_unlink(path: Path):
    """Silently remove a file if existing."""
    try:
        if path and path.exists():
            path.unlink()
    except Exception:
        pass


# ---------- Core per-barcode matching ----------

def process_one_barcode(
    bc: str,
    bid: str,
    pool_in: Path,
    out_dir: Path,
    seqkit_cmd: str,
    search_range: str,
    threads: int = 1,
    write_forward: bool = False,
    write_reverse: bool = False,
) -> tuple[int, int, Path, Path]:
    """
    Process one barcode on the current working pool.

    Steps (case-insensitive):
      1) Forward match on pool_in -> forward_tmp.fastq
         Extract forward IDs -> ids_fwd.txt
      2) Exclude forward IDs from pool_in -> pool_wo_fwd.fastq
      3) Reverse-complement pool_wo_fwd -> pool_wo_fwd.rc.fastq (DNA)
      4) Reverse match on RC pool -> reverse_tmp.fastq
         Extract reverse IDs -> ids_rev.txt
      5) Combine forward_tmp + reverse_tmp -> <ID>.fastq
      6) IDs_this_round = unique(ids_fwd ∪ ids_rev) -> ids_round.uniq.txt
         Next pool = pool_in minus ONLY this round of matched IDs -> pool_after_<ID>.fastq

    Returns:
      forward_count, reverse_count, ids_round_file, next_pool_path
    """
    bc = bc.upper()  # keep the pattern uppercase; grep uses -i anyway
    out_dir.mkdir(parents=True, exist_ok=True)

    # Temporary files for this barcode (all under per-ID folder)
    fwd_tmp = out_dir / f"forward_tmp_{bid}.fastq"
    fwd_ids = out_dir / f"ids_fwd_{bid}.txt"
    pool_wo_fwd = out_dir / f"pool_wo_fwd_{bid}.fastq"
    pool_wo_fwd_rc = out_dir / f"pool_wo_fwd_{bid}.rc.fastq"
    rev_tmp = out_dir / f"reverse_tmp_{bid}.fastq"
    rev_ids = out_dir / f"ids_rev_{bid}.txt"
    ids_round = out_dir / f"ids_round_{bid}.txt"
    ids_round_uniq = out_dir / f"ids_round_{bid}.uniq.txt"
    final_fastq = out_dir / f"{bid}.fastq"
    next_pool = out_dir / f"pool_after_{bid}.fastq"

    logging.info(f"[{bid}] Matching on working pool: {pool_in.name}")

    # 1) Forward match (case-insensitive)
    subprocess.run(
        [seqkit_cmd, 'grep', '-j', str(threads),
         '-R', str(search_range), '-s', '-P', '-i', '-p', bc,
         str(pool_in), '-o', str(fwd_tmp)],
        check=True
    )

    # Extract forward IDs (names only)
    with open(fwd_ids, 'w') as fh:
        subprocess.run([seqkit_cmd, 'seq', '-n', str(fwd_tmp)], check=True, stdout=fh)
    fwd_count = file_line_count(fwd_ids)

    # 2) Exclude forward IDs from pool_in -> pool_wo_fwd
    if fwd_count > 0:
        subprocess.run(
            [seqkit_cmd, 'grep', '-j', str(threads),
             '-n', '-f', str(fwd_ids), '-v',
             str(pool_in), '-o', str(pool_wo_fwd)],
            check=True
        )
    else:
        # No forward hits, simply copy pool_in to pool_wo_fwd
        shutil.copyfile(pool_in, pool_wo_fwd)

    # 3) Reverse-complement pool_wo_fwd (DNA) -> pool_wo_fwd.rc.fastq
    subprocess.run(
        [seqkit_cmd, 'seq', '-t', 'DNA', '-r', '-p',
         str(pool_wo_fwd), '-o', str(pool_wo_fwd_rc)],
        check=True
    )

    # 4) Reverse match on RC pool (case-insensitive)
    subprocess.run(
        [seqkit_cmd, 'grep', '-j', str(threads),
         '-R', str(search_range), '-s', '-P', '-i', '-p', bc,
         str(pool_wo_fwd_rc), '-o', str(rev_tmp)],
        check=True
    )

    # Extract reverse IDs
    with open(rev_ids, 'w') as rh:
        subprocess.run([seqkit_cmd, 'seq', '-n', str(rev_tmp)], check=True, stdout=rh)
    rev_count = file_line_count(rev_ids)

    # 5) Combine forward + reverse to final <ID>.fastq (under per-ID folder)
    with open(final_fastq, 'wb') as out_f:
        if fwd_tmp.exists() and fwd_tmp.stat().st_size > 0:
            with open(fwd_tmp, 'rb') as fh:
                shutil.copyfileobj(fh, out_f)
        if rev_tmp.exists() and rev_tmp.stat().st_size > 0:
            with open(rev_tmp, 'rb') as rh:
                shutil.copyfileobj(rh, out_f)

    # Optionally keep forward_/reverse_ files
    if not write_forward:
        safe_unlink(fwd_tmp)
    else:
        dst = out_dir / f"forward_{bid}.fastq"
        if fwd_tmp.exists():
            try:
                dst.unlink()
            except Exception:
                pass
            fwd_tmp.rename(dst)

    if not write_reverse:
        safe_unlink(rev_tmp)
    else:
        dst = out_dir / f"reverse_{bid}.fastq"
        if rev_tmp.exists():
            try:
                dst.unlink()
            except Exception:
                pass
            rev_tmp.rename(dst)

    # 6) Build per-round unique ID list = forward IDs ∪ reverse IDs
    with open(ids_round, 'wb') as w:
        if fwd_ids.exists():
            with open(fwd_ids, 'rb') as fh:
                shutil.copyfileobj(fh, w)
        if rev_ids.exists():
            with open(rev_ids, 'rb') as rh:
                shutil.copyfileobj(rh, w)

    if ids_round.exists() and ids_round.stat().st_size > 0:
        subprocess.run(['sort', '-u', str(ids_round), '-o', str(ids_round_uniq)], check=True)
        # Next pool = pool_in minus ONLY this round of matched IDs
        subprocess.run(
            [seqkit_cmd, 'grep', '-j', str(threads),
             '-n', '-f', str(ids_round_uniq), '-v',
             str(pool_in), '-o', str(next_pool)],
            check=True
        )
    else:
        # No matches this round -> carry pool forward
        shutil.copyfile(pool_in, next_pool)

    # Cleanup round-specific temps (keep final_fastq, next_pool)
    safe_unlink(fwd_ids)
    safe_unlink(rev_ids)
    safe_unlink(ids_round)
    safe_unlink(ids_round_uniq)
    safe_unlink(pool_wo_fwd)
    safe_unlink(pool_wo_fwd_rc)

    logging.info(f"[{bid}] Forward={fwd_count}, Reverse={rev_count}, Total={fwd_count + rev_count} -> {final_fastq}")
    return fwd_count, rev_count, final_fastq, next_pool


# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(description="Demultiplex FASTQ by barcodes with per-round mutual exclusion")
    parser.add_argument('-c', '--config', default='config.yaml', help='Path to YAML config file')
    args = parser.parse_args()

    cfg = load_config(args.config)
    setup_logging()

    fastq_path      = Path(cfg['fastq_path']).resolve()
    demux_table     = Path(cfg['demux_table']).resolve()
    output_dir      = Path(cfg['output_dir']).resolve()
    seqkit_cmd      = cfg['seqkit_cmd']
    search_range    = str(cfg['search_range'])
    threads         = int(cfg.get('threads', 1))
    write_forward   = bool(cfg.get('output_forward', False))
    write_reverse   = bool(cfg.get('output_reverse', False))
    write_unmatched = bool(cfg.get('output_unmatched', False))

    if not fastq_path.exists():
        raise FileNotFoundError(f"FASTQ not found: {fastq_path}")

    df = pd.read_csv(demux_table)
    if not {'ID', 'BC'}.issubset(df.columns):
        raise ValueError("Demux table must contain columns: ID, BC")

    output_dir.mkdir(parents=True, exist_ok=True)

    total_reads = count_total_reads_fastq(seqkit_cmd, fastq_path)
    logging.info(f"Total reads: {total_reads}")

    # Initialize working pool as the original FASTQ (flat under output_dir root)
    pool_curr = output_dir / "pool_initial.fastq"
    if pool_curr.exists():
        pool_curr.unlink()
    shutil.copyfile(fastq_path, pool_curr)

    forward_counts = {}
    reverse_counts = {}

    # Per-barcode processing with mutual exclusion between barcodes
    for _, row in df.iterrows():
        bid = str(row['ID'])
        bc  = str(row['BC']).upper()

        # >>> CHANGE: create per-sample subfolder under output_dir
        out_dir = output_dir / bid
        out_dir.mkdir(parents=True, exist_ok=True)

        fwd, rev, _final, pool_next = process_one_barcode(
            bc=bc,
            bid=bid,
            pool_in=pool_curr,
            out_dir=out_dir,                 # per-ID folder for this round
            seqkit_cmd=seqkit_cmd,
            search_range=search_range,
            threads=threads,
            write_forward=write_forward,
            write_reverse=write_reverse,
        )
        forward_counts[bid] = fwd
        reverse_counts[bid] = rev

        # advance working pool (may live in this per-ID folder; path is carried forward)
        safe_unlink(pool_curr)
        pool_curr = pool_next

    # At the end, pool_curr contains reads unmatched by all barcodes
    unmatched_fastq = output_dir / 'unmatched.fastq'
    if write_unmatched:
        if unmatched_fastq.exists():
            unmatched_fastq.unlink()
        shutil.move(str(pool_curr), str(unmatched_fastq))
        unmatched_count = count_total_reads_fastq(seqkit_cmd, unmatched_fastq)
        logging.info(f"Unmatched reads: {unmatched_count}")
    else:
        safe_unlink(pool_curr)

    # Write summary CSV at output_dir root
    stats = pd.DataFrame({
        'ID': df['ID'],
        'Forward': [int(forward_counts.get(i, 0)) for i in df['ID']],
        'Reverse': [int(reverse_counts.get(i, 0)) for i in df['ID']],
    })
    stats['Total'] = stats['Forward'] + stats['Reverse']
    stats.to_csv(output_dir / 'BC_Match_Statistics.csv', index=False)
    logging.info("BC match statistics saved.")
    logging.info("Done.")


if __name__ == '__main__':
    main()
