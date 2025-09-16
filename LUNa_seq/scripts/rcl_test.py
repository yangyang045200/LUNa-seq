#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rcl_test.py

Recombination/contamination-like (RCL) integration check per sample, YAML-driven.

Purpose
-------
Given a query FASTA (~5.7 kb), detect whether it (partly or fully) integrates into:
  - contaminant/recombinant merged references (CR)
  - non-specific residual products (per consensus)

Positivity rule
---------------
Positive if ANY contiguous alignment block (M/= /X; small I/D tolerated) has
effective length >= threshold (bp). Small indels (I/D) with length <= indel_tol
do NOT break the block and DO contribute to the "effective block length".

Defaults aim for high recall:
  minimap2 args: "-a -x asm5 -N 100000 --secondary=yes -p 0"
  keep primary + secondary + supplementary; filter by mapq_min if > 0.

Outputs per sample (saved directly in sample/<ID>/):
  - RCL_summary.tsv
  - <hits_tsv> per target (configurable; defaults provided)
  - RCL_from_cr.png
  - RCL_from_non_specific.png

Everything is controlled by config.yaml under the `rcl_test:` section.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
from pathlib import Path
from typing import List, Tuple, Dict, Optional

import yaml
import pysam
from Bio import SeqIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --------------------------
# Logging & Config
# --------------------------
def setup_logging(level=logging.INFO):
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s: %(message)s")


def load_config(config_path: Path) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


# --------------------------
# Minimap2 & FASTA helpers
# --------------------------
def run_minimap2_sam(ref_fa: Path, query_fa: Path, sam_out: Path, mm2_args: List[str]):
    """Run minimap2 to produce SAM. Ensure '-a' is present."""
    if "-a" not in mm2_args:
        mm2_args = ["-a", *mm2_args]
    cmd = ["minimap2", *mm2_args, str(ref_fa), str(query_fa)]
    logging.info("Running: %s > %s", " ".join(cmd), sam_out)
    with open(sam_out, "w") as fout:
        res = subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        logging.error("minimap2 failed: %s", res.stderr)
        raise RuntimeError("minimap2 failed")


def read_fasta_lengths(fa: Path) -> Dict[str, int]:
    out = {}
    for rec in SeqIO.parse(str(fa), "fasta"):
        out[rec.id] = len(rec.seq)
    return out


def count_fasta_records(fa: Path) -> int:
    if not fa.exists():
        return 0
    return sum(1 for _ in SeqIO.parse(str(fa), "fasta"))


# --------------------------
# CIGAR parsing with tolerance
# --------------------------
def extract_blocks_with_tolerance(aln: pysam.AlignedSegment, indel_tol: int) -> List[Tuple[int, int, int]]:
    """
    Extract contiguous blocks on the reference with small indels tolerated.
    Returns [(ref_start_0b, ref_end_0b_excl, eff_len)].
    eff_len = M/= /X + tolerated small I + tolerated small D.
    """
    blocks = []
    ref_pos = aln.reference_start
    cigar = aln.cigartuples or []

    cur_start = None
    cur_ref_end = None   # exclusive end on ref
    eff_len = 0

    def close_block():
        nonlocal cur_start, cur_ref_end, eff_len
        if cur_start is not None and cur_ref_end is not None and eff_len > 0:
            blocks.append((cur_start, cur_ref_end, eff_len))
        cur_start = None
        cur_ref_end = None
        eff_len = 0

    for op, length in cigar:
        if op in (0, 7, 8):  # M, =, X
            if cur_start is None:
                cur_start = ref_pos
                cur_ref_end = ref_pos + length
                eff_len = length
            else:
                cur_ref_end += length
                eff_len += length
            ref_pos += length

        elif op == 2:  # D (advance ref)
            if length <= indel_tol:
                if cur_start is None:
                    cur_start = ref_pos
                    cur_ref_end = ref_pos + length
                    eff_len = length
                else:
                    cur_ref_end += length
                    eff_len += length
                ref_pos += length
            else:
                close_block()
                ref_pos += length

        elif op == 1:  # I (advance query)
            if length <= indel_tol:
                if cur_start is None:
                    cur_start = ref_pos
                    cur_ref_end = ref_pos
                    eff_len = 0
                eff_len += length
            else:
                close_block()

        elif op in (3, 4, 5, 6):  # N/S/H/P
            close_block()

        else:
            close_block()

    close_block()
    return blocks


def find_positive_events(ref_fa: Path,
                         query_fa: Path,
                         threshold: int,
                         indel_tol: int,
                         mm2_args: List[str],
                         mapq_min: int) -> Tuple[Dict[str, List[Tuple[int, int]]], Dict[str, int]]:
    """
    Align query to ref, parse SAM, and collect all contiguous aligned blocks
    with effective length >= threshold (small I/D tolerated).

    KEEP primary + secondary + supplementary (skip only unmapped).
    Apply MAPQ threshold if mapq_min > 0.

    Returns:
      hits_by_target: {target_id: [(start0, end0_excl), ...]}
      target_lengths: {target_id: length}
    """
    if not ref_fa.exists():
        return {}, {}

    target_lengths = read_fasta_lengths(ref_fa)
    sam_path = ref_fa.with_suffix(ref_fa.suffix + ".query.sam")
    run_minimap2_sam(ref_fa, query_fa, sam_path, mm2_args=mm2_args)

    hits_by_target: Dict[str, List[Tuple[int, int]]] = {}
    total_alns = 0
    kept_mapq = 0
    kept_blocks = 0

    sf = pysam.AlignmentFile(str(sam_path), "r")
    for aln in sf.fetch(until_eof=True):
        total_alns += 1
        if aln.is_unmapped:
            continue
        if aln.mapping_quality < mapq_min:
            continue
        kept_mapq += 1

        tname = sf.get_reference_name(aln.reference_id)
        blocks = extract_blocks_with_tolerance(aln, indel_tol=indel_tol)
        blocks = [(s, e) for (s, e, eff) in blocks if eff >= threshold]
        if not blocks:
            continue
        kept_blocks += len(blocks)
        arr = hits_by_target.get(tname, [])
        arr.extend(blocks)
        hits_by_target[tname] = arr
    sf.close()

    try:
        sam_path.unlink()
    except Exception:
        pass

    # De-duplicate intervals per target
    for t in list(hits_by_target.keys()):
        uniq = sorted(set(hits_by_target[t]))
        hits_by_target[t] = uniq

    logging.info("Align stats for %s: total=%d, pass_mapq=%d, blocks(>=thr)=%d, targets_with_hits=%d",
                 ref_fa.name, total_alns, kept_mapq, kept_blocks, len(hits_by_target))
    return hits_by_target, target_lengths


# --------------------------
# Plotting
# --------------------------
def plot_hits(out_png: Path, hits_by_target: Dict[str, List[Tuple[int, int]]],
              target_lengths: Dict[str, int], title: str):
    """Save summary plot of positive intervals per target. PNG is placed in sample/<ID>/."""
    if not hits_by_target:
        plt.figure(figsize=(6, 2))
        plt.title(f"{title} (No positive events)")
        plt.axis("off")
        out_png.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()
        return

    targets = sorted(hits_by_target.keys())
    n = len(targets)
    height = max(2.0, 0.4 + n * 0.35 + 0.4)
    fig, ax = plt.subplots(figsize=(12, height))

    for i, tname in enumerate(targets):
        y = n * 0.35 - i * 0.35
        tlen = target_lengths.get(tname, 0)
        ax.hlines(y, 1, max(tlen, 1), linewidth=1)
        for (s0, e0) in sorted(hits_by_target[tname]):
            s1, e1 = s0 + 1, e0
            blen = max(e0 - s0, 1)
            ax.add_patch(plt.Rectangle((s1, y - 0.06), blen, 0.12, alpha=0.7))
            # Optional label next to track:
            ax.text(0, y + 0.12, f"{tname}  {s1}-{e1} ({blen})",
                    ha="left", va="bottom", fontsize=7)

    xmax = max((target_lengths.get(t, 0) for t in targets), default=0) + 5
    ax.set_xlim(0, max(5, xmax))
    ax.set_ylim(-0.2, n * 0.35 + 0.4)
    ax.set_yticks([])
    ax.set_xlabel("Target coordinate (1-based)")
    ax.set_title(title)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


# --------------------------
# Core
# --------------------------
def _resolve_with_ext(p: Path) -> Path:
    """Return existing path; if missing and no suffix, try '.fasta'."""
    if p.exists():
        return p
    if p.suffix == "":
        alt = p.with_suffix(".fasta")
        if alt.exists():
            return alt
    return p


def process_sample(sample_dir: Path, query_fa: Path, cfg: dict):
    """
    For each sample folder:
      - Count consensus reads (for denominator)
      - For each configured target:
          align -> collect positive blocks -> write hits tsv -> plot -> append summary
    """
    sample = sample_dir.name
    rcl_cfg = cfg.get("rcl_test", {})

    threshold = int(rcl_cfg.get("threshold", 100))
    indel_tol = int(rcl_cfg.get("indel_tol", 5))
    mapq_min = int(rcl_cfg.get("mapq_min", 0))
    mm2_args = (rcl_cfg.get("minimap2_args") or "-a -x asm5 -N 100000 --secondary=yes -p 0").split()

    # Default inputs/outputs if not provided
    targets = rcl_cfg.get("targets") or [
        {
            "label": "cr",
            "path": "contaminant_or_recombination_judgement_raw_data.fasta",
            "hits_tsv": "RCL_from_cr_hits.tsv",
            "plot_png": "RCL_from_cr.png",
        },
        {
            "label": "non_specific",
            "path": "non_specific_residual_products__consensus_reads_ge_5.fasta",
            "hits_tsv": "RCL_from_non_specific_hits.tsv",
            "plot_png": "RCL_from_non_specific.png",
        },
    ]

    consensus_fa = sample_dir / "consensus_reads_ge_5.fasta"
    consensus_count = count_fasta_records(consensus_fa)

    summary_rows = []
    summary_tsv = sample_dir / "RCL_summary.tsv"

    for t in targets:
        label = t.get("label", "target")
        tgt_rel = t["path"]
        hits_name = t.get("hits_tsv", f"{label}_hits.tsv")
        plot_name = t.get("plot_png", f"{label}.png")

        tgt_path = _resolve_with_ext(sample_dir / tgt_rel)
        hits_tsv = sample_dir / hits_name
        plot_png = sample_dir / plot_name

        if not tgt_path.exists():
            logging.warning("[%s] Missing target file: %s", sample, tgt_path.name)
            hits_tsv.write_text("target_id\tstart_1b\tend_1b\tlength\n")
            summary_rows.append([label, "No", 0, 0, consensus_count, 0.0])
            plot_hits(plot_png, {}, {}, f"{sample} - {label} (No target file)")
            continue

        hits_by_target, target_lengths = find_positive_events(
            ref_fa=tgt_path,
            query_fa=query_fa,
            threshold=threshold,
            indel_tol=indel_tol,
            mm2_args=mm2_args,
            mapq_min=mapq_min
        )

        with open(hits_tsv, "w") as w:
            w.write("target_id\tstart_1b\tend_1b\tlength\n")
            for tid, blocks in sorted(hits_by_target.items()):
                for (s0, e0) in sorted(blocks):
                    s1, e1 = s0 + 1, e0
                    blen = e0 - s0
                    w.write(f"{tid}\t{s1}\t{e1}\t{blen}\n")

        positive_target_count = len(hits_by_target)             # targets with >=1 positive block
        positive_segment_count = sum(len(v) for v in hits_by_target.values())
        presence = "Yes" if positive_target_count > 0 else "No"
        ratio = (positive_target_count / consensus_count) if consensus_count else 0.0

        summary_rows.append([label, presence, positive_target_count, positive_segment_count, consensus_count, ratio])

        plot_hits(plot_png, hits_by_target, target_lengths,
                  title=f"{sample} - {label} (≥{threshold} bp, indel_tol={indel_tol})")

    with open(summary_tsv, "w") as w:
        w.write("file\tpositive\tpositive_target_count\tpositive_segment_count\tconsensus_reads_count\tratio_targets_over_consensus\n")
        for r in summary_rows:
            w.write("\t".join([str(x) for x in r]) + "\n")

    logging.info("[%s] RCL done. Summary -> %s", sample, summary_tsv)


# --------------------------
# CLI
# --------------------------
def main():
    ap = argparse.ArgumentParser(description="RCL integration check per sample (YAML-driven).")
    ap.add_argument("input_folder", help="Parent folder containing sample subdirectories (e.g., ./sample)")
    ap.add_argument("-c", "--config", default="config.yaml", help="YAML config file path")
    args = ap.parse_args()

    setup_logging()
    cfg = load_config(Path(args.config))
    config_dir = Path(args.config).resolve().parent

    base = Path(args.input_folder).resolve()
    if not base.exists():
        raise FileNotFoundError(f"Samples dir not found: {base}")

    rcl_cfg = cfg.get("rcl_test", {})
    query_fa_cfg = rcl_cfg.get("query_fasta")
    if not query_fa_cfg:
        raise ValueError("rcl_test.query_fasta must be set in config.yaml")
    # Resolve query path relative to config.yaml directory for portability
    query_fa = (config_dir / query_fa_cfg).resolve() if not Path(query_fa_cfg).is_absolute() else Path(query_fa_cfg)
    if not query_fa.exists():
        raise FileNotFoundError(f"Query FASTA not found: {query_fa}")

    for sample_dir in sorted(p for p in base.iterdir() if p.is_dir()):
        process_sample(sample_dir, query_fa, cfg)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logging.error("Interrupted by user.")
        raise SystemExit(130)
