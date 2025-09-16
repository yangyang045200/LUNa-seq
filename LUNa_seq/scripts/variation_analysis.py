#!/usr/bin/env python3
# variation_analysis.py
#
# - Per-sample: align each input FASTA/FASTQ to its reference (minimap2+samtools)
# - Build per-base event matrix CSV and whole-length summaries + plot
# - Supports per-input minimap2 args via YAML `query_inputs`
# - Config can be top-level or nested under `variation_analysis:`
# - NEW: statistics & plots ignore the first/last EDGE_IGNORE positions.

import argparse
import yaml
import subprocess
import logging
from pathlib import Path
import pysam
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv
from typing import Iterable, List, Tuple, Dict, Any

# how many positions to ignore at each end in statistics/plots
EDGE_IGNORE = 5

def setup_logging(log_file: Path):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )

def run_command(cmd: List[str]):
    logging.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(result.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")

def _ensure_a(args_list: Iterable[str]) -> List[str]:
    args_list = list(args_list)
    return args_list if "-a" in args_list else ["-a", *args_list]

def align(ref: Path, query: Path, outdir: Path, minimap2_args: Iterable[str] = ("-a",)) -> Path:
    sam = outdir / f"{query.stem}.sam"
    bam = outdir / f"{query.stem}.bam"
    sorted_bam = outdir / f"{query.stem}.sorted.bam"

    mm_args = _ensure_a(minimap2_args)
    run_command(["minimap2", *mm_args, str(ref), str(query), "-o", str(sam)])
    run_command(["samtools", "view", "-b", "-o", str(bam), str(sam)])
    run_command(["samtools", "sort", "-o", str(sorted_bam), str(bam)])
    run_command(["samtools", "index", str(sorted_bam)])

    for f in (sam, bam):
        try: f.unlink()    
        except Exception: pass
    return sorted_bam

def generate_event_matrix(bam_path: Path, ref_fasta: Path, output_csv: Path):
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    fasta = pysam.FastaFile(str(ref_fasta))

    ref_name = fasta.references[0]
    ref_len = fasta.get_reference_length(ref_name)
    positions = list(range(1, ref_len + 1))

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name"] + [f"pos_{i}" for i in positions])
        # reference row must be uppercase "REFERENCE"
        w.writerow(["REFERENCE"] + [fasta.fetch(ref_name, i - 1, i).upper() for i in positions])

        for read in bam.fetch(ref_name, 0, ref_len):
            if read.is_unmapped:
                continue
            seq = read.query_sequence
            if not seq:
                logging.warning("Read %s has no sequence; skipping.", read.query_name)
                continue

            seqU = seq.upper()
            events = ["." for _ in positions]
            prev_pr = None
            ins_buf: List[str] = []

            for qr, pr in read.get_aligned_pairs(matches_only=False):
                # insertion relative to reference (pr None, qr not None)
                if pr is None and qr is not None:
                    if 0 <= qr < len(seqU):
                        ins_buf.append(seqU[qr])
                    continue

                # flush any buffered insertion to the anchor just after prev_pr
                if ins_buf:
                    anchor_1b = (prev_pr if prev_pr is not None else -1) + 1
                    # still record insertions at all positions; filtering happens in analyze_events
                    if 1 <= anchor_1b <= ref_len:
                        events[anchor_1b - 1] = f"INS:{''.join(ins_buf)}"
                    ins_buf = []

                # deletion relative to reference (qr None, pr not None)
                if qr is None and pr is not None:
                    p1 = pr + 1
                    events[p1 - 1] = "DEL"
                    prev_pr = pr

                # match/mismatch base (both not None)
                elif qr is not None and pr is not None:
                    p1 = pr + 1
                    refb = fasta.fetch(ref_name, pr, pr + 1).upper()
                    qb = seqU[qr]
                    # === CHANGE: treat 'N' as 'N' (unknown), not as mismatch ===
                    if qb == "N":
                        events[p1 - 1] = "N"
                    elif refb != qb:
                        events[p1 - 1] = qb
                    prev_pr = pr

            # flush tail insertion if exists
            if ins_buf:
                anchor_1b = (prev_pr if prev_pr is not None else -1) + 1
                if 1 <= anchor_1b <= ref_len:
                    events[anchor_1b - 1] = f"INS:{''.join(ins_buf)}"

            w.writerow([read.query_name] + events)

    bam.close()
    fasta.close()
    logging.info("Wrote events to %s", output_csv)

def _kept_position_columns(df_reads: pd.DataFrame) -> Tuple[List[str], int]:
    """Return list of position columns after removing first/last EDGE_IGNORE, and ref_len."""
    # columns like POS_1, POS_2, ...
    positions = df_reads.columns.tolist()
    # parse numeric indices
    nums = []
    for c in positions:
        try:
            nums.append(int(str(c).split("_", 1)[1]))
        except Exception:
            pass
    ref_len = max(nums) if nums else len(positions)

    left_cut = set(range(1, min(EDGE_IGNORE, ref_len) + 1))
    right_cut_start = max(ref_len - EDGE_IGNORE + 1, 1)
    right_cut = set(range(right_cut_start, ref_len + 1))
    ignore_set = left_cut | right_cut

    kept = []
    for c in positions:
        try:
            n = int(str(c).split("_", 1)[1])
            if n not in ignore_set:
                kept.append(c)
        except Exception:
            # if a weird column sneaks in, keep it out of stats
            pass
    return kept, ref_len

def analyze_events(evt_csv: Path, mutation_csv: Path, plot_file: Path):
    df = pd.read_csv(evt_csv)
    df.columns = [c.upper() for c in df.columns]
    if "NAME" not in df.columns:
        raise ValueError(f"{evt_csv} is missing 'NAME' column.")
    df["NAME"] = df["NAME"].astype(str).str.upper()

    for c in df.columns[1:]:
        df[c] = df[c].astype(str).str.upper()

    if (df["NAME"] == "REFERENCE").sum() == 0:
        raise ValueError(f"No 'REFERENCE' row found in {evt_csv}. NAME head={df['NAME'].head().tolist()}")

    ref_series = df[df["NAME"] == "REFERENCE"].iloc[0, 1:]
    df_reads = df[df["NAME"] != "REFERENCE"].iloc[:, 1:]

    # keep only inner positions (drop first/last EDGE_IGNORE)
    kept_cols, ref_len = _kept_position_columns(df_reads)

    n_reads = df_reads.shape[0]
    n_positions_eff = len(kept_cols)
    total_bases_eff = n_reads * n_positions_eff if n_positions_eff else 0

    # counts on kept columns only
    ins = sum(df_reads[c].str.startswith("INS:").sum() for c in kept_cols)
    dele = sum((df_reads[c] == "DEL").sum() for c in kept_cols)

    mm_count = 0
    for c in kept_cols:
        s = (~df_reads[c].isin([".", "DEL", "N"])) & (~df_reads[c].str.startswith("INS:"))
        mm_count += int(s.sum())

    nn = sum((df_reads[c] == "N").sum() for c in kept_cols)

    ins_rate = ins / total_bases_eff if total_bases_eff else 0.0
    del_rate = dele / total_bases_eff if total_bases_eff else 0.0
    mm_rate = mm_count / total_bases_eff if total_bases_eff else 0.0
    n_rate  = nn  / total_bases_eff if total_bases_eff else 0.0

    with open(mutation_csv, "w") as fh:
        fh.write(
            "sample,total_reads,total_bases_eff,insertion_count,deletion_count,"
            "mismatch_count,N_count,insertion_rate,deletion_rate,mismatch_rate,N_rate\n"
        )
        fh.write(
            f"{Path(evt_csv).stem},{n_reads},{total_bases_eff},{ins},{dele},"
            f"{mm_count},{nn},{ins_rate},{del_rate},{mm_rate},{n_rate}\n"
        )

        # mismatch spectrum on kept columns only
        mm_types: Dict[str, int] = {}
        for c in kept_cols:
            refb = str(ref_series.get(c, "N")).upper()
            col = df_reads[c]
            mask = (~col.isin([".", "DEL", "N"])) & (~col.str.startswith("INS:"))
            for val in col[mask]:
                key = f"{refb}_to_{val}"
                mm_types[key] = mm_types.get(key, 0) + 1
        total_mm = sum(mm_types.values())
        fh.write("\n# full_length mismatch_type,count,proportion\n")
        if total_mm > 0:
            for k, v in sorted(mm_types.items()):
                fh.write(f"{k},{v},{v/total_mm:.8f}\n")
        else:
            fh.write("NA,0,0.00000000\n")

    logging.info("Wrote mutation summary to %s", mutation_csv)

    # sliding-window plot on kept columns only
    all_pos_cols = df_reads.columns.tolist()
    # map position index → column name
    pos_idx_to_col = {}
    for c in all_pos_cols:
        try:
            n = int(str(c).split("_", 1)[1])
            pos_idx_to_col[n] = c
        except Exception:
            pass

    ref_len_num = max(pos_idx_to_col.keys()) if pos_idx_to_col else len(all_pos_cols)
    kept_set = set(kept_cols)

    window_size = 50
    idx, ins_rates, del_rates, mm_rates = [], [], [], []
    for start1 in range(1, ref_len_num + 1, window_size):
        end1 = min(start1 + window_size - 1, ref_len_num)
        window_cols = []
        for n in range(start1, end1 + 1):
            c = pos_idx_to_col.get(n)
            if c and c in kept_set:
                window_cols.append(c)

        denom = n_reads * len(window_cols) if window_cols else 0
        if denom == 0:
            idx.append(start1)
            ins_rates.append(0.0)
            del_rates.append(0.0)
            mm_rates.append(0.0)
            continue

        ins_w = sum(df_reads[w].str.startswith("INS:").sum() for w in window_cols)
        del_w = sum((df_reads[w] == "DEL").sum() for w in window_cols)
        mm_w = 0
        for w in window_cols:
            s = (~df_reads[w].isin([".", "DEL", "N"])) & (~df_reads[w].str.startswith("INS:"))
            mm_w += int(s.sum())

        idx.append(start1)
        ins_rates.append(ins_w / denom)
        del_rates.append(del_w / denom)
        mm_rates.append(mm_w / denom)

    plt.figure()
    plt.bar(idx, mm_rates, width=window_size, label="Mismatch")
    plt.bar(idx, del_rates, width=window_size, bottom=mm_rates, label="Deletion")
    plt.bar(idx, ins_rates, width=window_size,
            bottom=[x + y for x, y in zip(mm_rates, del_rates)],
            label="Insertion")
    plt.xlabel("Position")
    plt.ylabel("Event rate")
    plt.title(f"Stacked variant event rates (50 bp windows, edges ±{EDGE_IGNORE} ignored)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_file)
    logging.info("Plot saved to %s", plot_file)

def _collect_inputs_with_args(sample_dir: Path, specs: List[Tuple[str, List[str]]]) -> List[Tuple[Path, List[str]]]:
    out: List[Tuple[Path, List[str]]] = []
    seen = set()
    for pat, mm_args in specs:
        for p in sample_dir.glob(pat):
            if p.is_file():
                rp = p.resolve()
                if rp not in seen:
                    out.append((p, mm_args))
                    seen.add(rp)
    return out

def _load_cfg(cfg_raw: Dict[str, Any]) -> Dict[str, Any]:
    if "variation_analysis" in cfg_raw and isinstance(cfg_raw["variation_analysis"], dict):
        return cfg_raw["variation_analysis"]
    return cfg_raw

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="YAML config path")
    args = parser.parse_args()

    cfg_all = yaml.safe_load(open(args.config))
    cfg = _load_cfg(cfg_all)

    samples_dir = Path(cfg["samples_dir"])
    rp = cfg.get("ref_pattern", "{folder}.fasta")

    specs: List[Tuple[str, List[str]]] = []
    for ent in cfg.get("query_inputs", []) or []:
        if isinstance(ent, str):
            mm = cfg.get("minimap2_args", "-a")
            specs.append((ent, mm.split()))
        elif isinstance(ent, dict):
            pat = ent.get("pattern")
            if not pat:
                continue
            mm = ent.get("minimap2_args", cfg.get("minimap2_args", "-a"))
            specs.append((pat, mm.split()))
    for pat in cfg.get("query_files", []) or []:
        mm = cfg.get("minimap2_args", "-a")
        specs.append((pat, mm.split()))
    if not specs:
        raise ValueError("Please provide 'query_inputs' or 'query_files' in the YAML.")

    setup_logging(Path(samples_dir) / cfg.get("log", "variation.log"))
    logging.info("Ignoring first/last %d reference positions in statistics & plots.", EDGE_IGNORE)

    for folder in samples_dir.iterdir():
        if not folder.is_dir():
            continue
        sample = folder.name
        ref_file = folder / rp.format(folder=sample)
        if not ref_file.exists():
            logging.warning("Reference not found: %s", ref_file)
            continue

        inputs = _collect_inputs_with_args(folder, specs)
        if not inputs:
            logging.info("[%s] No inputs matched any of %s", sample, [s[0] for s in specs])
            continue

        for qf, mm_args in inputs:
            logging.info("[%s] Processing %s  (minimap2 args: %s)", sample, qf.name, " ".join(mm_args))
            bam = align(ref_file, qf, folder, minimap2_args=mm_args)

            stem_tag = f"{sample}_{qf.stem}"
            evt_csv = folder / f"events_{stem_tag}.csv"
            mut_csv = folder / f"mutation_{stem_tag}.csv"
            plot_file = folder / f"plot_{stem_tag}.png"

            generate_event_matrix(bam, ref_file, evt_csv)
            analyze_events(evt_csv, mut_csv, plot_file)

            for f in (bam, bam.with_suffix(bam.suffix + ".bai")):
                try: f.unlink()
                except Exception: pass

    logging.info("All samples processed.")

if __name__ == "__main__":
    main()
