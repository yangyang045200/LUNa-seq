#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
visualize_igv_snapshots.py
Batch IGV visualization for per-sample subfolders, driven by YAML config.

Changes from the standalone prototype:
- All tunables (IGV/minimap2/threads/downsampling/merge patterns/inputs) are in YAML.
- PNG snapshots are saved DIRECTLY into each sample subfolder (requested behavior).
- BAM/session/log artifacts go to a shared "artifact_dir" (default: ./visualization/<ID>).
- Raw merge EXCLUDES files containing '__consensus_reads_ge_5'.

Original core logic is preserved:
minimap2 -> samtools -> two IGV runs -> snapshot of full reference locus (first contig).
"""

from __future__ import annotations

import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple, Dict
import random
import yaml
import xml.etree.ElementTree as ET

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# -------------------------
# Logging
# -------------------------
def setup_logging(level=logging.INFO):
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s: %(message)s")


# -------------------------
# Config
# -------------------------
def load_config(config_path: Path) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


# -------------------------
# Reference resolution (same pattern as your other script)
# -------------------------
def resolve_reference(sample_dir: Path, sample_name: str, ref_folder: str, config_dir: Path) -> Optional[Path]:
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
        dest.write_bytes(candidate.read_bytes())
        logging.info("[%s] Copied reference from %s to %s", sample_name, candidate, dest)
        return dest

    logging.warning("[%s] Reference not found: %s OR %s", sample_name, local, candidate)
    return None


# -------------------------
# Utilities
# -------------------------
def run(cmd: List[str], **kw):
    logging.info("[RUN] %s", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True, **kw)


def ensure_dir(path: Path):
    if path.exists() and not path.is_dir():
        bak = path.with_suffix(path.suffix + ".bak")
        logging.warning("%s exists and is a file. Renaming to %s", path, bak)
        path.rename(bak)
    path.mkdir(parents=True, exist_ok=True)


def ensure_fai(ref_fa: Path) -> Path:
    fai = ref_fa.with_suffix(ref_fa.suffix + ".fai")
    if not fai.exists():
        run(["samtools", "faidx", str(ref_fa)])
    return fai


def first_seq(fai: Path) -> Tuple[str, int]:
    with open(fai, "r") as f:
        line = f.readline().strip()
    if not line:
        raise RuntimeError(f"{fai} is empty")
    name, length = line.split("\t")[:2]
    return name, int(length)


def write_session_xml(ref_fa: Path, sorted_bam: Path, session_xml: Path):
    root = ET.Element("Session", {
        "genome": str(ref_fa.resolve()),
        "hasGeneTrack": "false",
        "hasSequenceTrack": "true",
        "version": "8"
    })
    res = ET.SubElement(root, "Resources")
    ET.SubElement(res, "Resource", {"path": str(sorted_bam.resolve())})
    ET.ElementTree(root).write(session_xml, encoding="UTF-8", xml_declaration=True)


def print_log_tail(log_path: Path, n=80):
    try:
        lines = log_path.read_text(errors="ignore").splitlines()
        tail = "\n".join(lines[-n:])
        logging.warning("---- IGV LOG TAIL ----\n%s\n-----------------------", tail)
    except Exception:
        pass


def get_bam_track_name(log_path: Path) -> str:
    txt = log_path.read_text(errors="ignore")
    for line in txt.splitlines():
        if ".bam" in line:
            return line.strip()
    raise RuntimeError("Failed to find BAM track name in IGV log.")


def detect_format(path: Path) -> str:
    suf = path.suffix.lower()
    if suf in [".fa", ".fasta"]:
        return "fasta"
    if suf in [".fq", ".fastq"]:
        return "fastq"
    return "fasta"


def sample_reads(input_path: Path, out_dir: Path, max_n: int) -> Path:
    """Random downsampling of fasta/fastq up to max_n records; keep all if <= max_n."""
    fmt = detect_format(input_path)
    records = list(SeqIO.parse(str(input_path), fmt))
    n = len(records)
    if n == 0:
        raise RuntimeError(f"No sequences found in {input_path}")
    if max_n and n > max_n:
        logging.info("[%s] %d records; randomly keeping %d", input_path.name, n, max_n)
        records = random.sample(records, max_n)
        out_path = out_dir / f"{input_path.stem}_subsampled{input_path.suffix}"
        SeqIO.write(records, str(out_path), fmt)
        return out_path
    else:
        logging.info("[%s] %d records; keeping all", input_path.name, n)
        return input_path


def merge_fastas(sample_dir: Path, include_glob: str, out_name: str, exclude_substrings: Optional[List[str]] = None) -> Optional[Path]:
    """
    Concatenate FASTA files matching include_glob into out_name, skipping any whose name
    contains a forbidden substring in exclude_substrings.
    """
    in_paths = sorted(sample_dir.glob(include_glob))
    if exclude_substrings:
        in_paths = [p for p in in_paths if not any(sub in p.name for sub in exclude_substrings)]
    if not in_paths:
        logging.info("[%s] No files match '%s'; skip merge to %s", sample_dir.name, include_glob, out_name)
        return None

    out_path = sample_dir / out_name
    logging.info("[%s] Merging %d files -> %s", sample_dir.name, len(in_paths), out_name)
    with open(out_path, "w") as w:
        for p in in_paths:
            txt = p.read_text()
            if not txt.endswith("\n"):
                txt += "\n"
            w.write(txt)
    return out_path


def igv_cmdline(igv_cmd: str, use_xvfb: bool, batch_path: Path, ref_fa: Path) -> List[str]:
    """
    Build the IGV command line. When use_xvfb is True, wraps IGV with 'xvfb-run -a'.
    """
    if use_xvfb:
        return ["xvfb-run", "-a", igv_cmd, "-g", str(ref_fa.resolve()), "-b", str(batch_path)]
    else:
        return [igv_cmd, "-g", str(ref_fa.resolve()), "-b", str(batch_path)]


# -------------------------
# Core
# -------------------------
def process_sample(sample_dir: Path, cfg: dict, config_dir: Path):
    vcfg = cfg.get("visualization_igv", {})
    sample = sample_dir.name

    # Reference resolution
    ref_fa = resolve_reference(sample_dir, sample, vcfg.get("ref_folder", "ref"), config_dir)
    if not ref_fa:
        return

    # Artifact directory for BAM/session/log
    artifact_root = Path(vcfg.get("artifact_dir", "visualization"))
    out_dir = artifact_root / sample
    ensure_dir(out_dir)

    # Merge raw & consensus contaminant/recombinant references
    merge_cfg = vcfg.get("merge", {})
    raw_cfg = merge_cfg.get("raw", {})
    cons_cfg = merge_cfg.get("consensus", {})

    merge_fastas(
        sample_dir,
        include_glob=raw_cfg.get("include_glob", "contaminant_or_recombinant_ref_*.fasta"),
        out_name=raw_cfg.get("out", "contaminant_or_recombination_judgement_raw_data.fasta"),
        exclude_substrings=raw_cfg.get("exclude", ["__consensus_reads_ge_5"]),
    )
    merge_fastas(
        sample_dir,
        include_glob=cons_cfg.get("include_glob", "contaminant_or_recombinant_ref_*__consensus_reads_ge_5.fasta"),
        out_name=cons_cfg.get("out", "contaminant_or_recombination_judgement_consensus.fasta"),
        exclude_substrings=cons_cfg.get("exclude", None),
    )

    # Inputs -> snapshots mapping (PNG saved into sample_dir)
    inputs_cfg = vcfg.get("inputs") or [
        {"path": "{ID}.fastq",                                                "snapshot": "{ID}_raw_data.png"},
        {"path": "fully_on_target_proviral_reads.fastq",                      "snapshot": "fully_on_target_raw_data.png"},
        {"path": "contaminant_or_recombination_judgement_raw_data.fasta",     "snapshot": "contaminant_or_recombination_judgement_raw_data.png"},
        {"path": "consensus_reads_ge_5.fasta",                                "snapshot": "{ID}_consensus.png"},
        {"path": "fully_on_target_proviral_reads__consensus_reads_ge_5.fasta","snapshot": "{ID}_fully_on_target_consensus.png"},
        {"path": "contaminant_or_recombination_judgement_consensus.fasta",    "snapshot": "contaminant_or_recombination_judgement_consensus.png"},
    ]

    # Runtime knobs
    mm_args = (vcfg.get("minimap2_args") or "-a -x map-ont").split()
    threads = vcfg.get("threads")
    if threads:
        # minimap2 supports -t; samtools sort below will also use -@ if provided
        pass
    sleep_ms = int(vcfg.get("igv", {}).get("sleep_ms", 15000))
    timeout_sec = int(vcfg.get("igv", {}).get("timeout_sec", 900))
    igv_cmd = vcfg.get("igv", {}).get("cmd") or os.environ.get("IGV_CMD") or os.environ.get("IGV_SH") or "igv.sh"
    use_xvfb = bool(vcfg.get("igv", {}).get("use_xvfb", True))
    max_reads = int(vcfg.get("max_reads", 2000))

    # Prepare reference locus
    fai = ensure_fai(ref_fa)
    contig, length = first_seq(fai)
    locus = f"{contig}:1-{length}"

    # Iterate inputs
    for item in inputs_cfg:
        in_path_tmpl = item["path"]
        snap_tmpl = item["snapshot"]

        in_path = in_path_tmpl.replace("{ID}", sample)
        snapshot_name = snap_tmpl.replace("{ID}", sample)

        reads_path = sample_dir / in_path
        png_path = sample_dir / snapshot_name

        if not reads_path.exists():
            logging.info("[%s] Missing input: %s (skip)", sample, reads_path.name)
            continue

        # Downsample if needed
        try:
            reads_used = sample_reads(reads_path, out_dir, max_reads=max_reads)
        except Exception as e:
            logging.warning("[%s] Sampling failed for %s (%s); using full file", sample, reads_path.name, e)
            reads_used = reads_path

        # Filenames for artifacts
        prefix = reads_used.stem
        sam_path     = out_dir / f"{sample}_{prefix}.sam"
        bam_path     = out_dir / f"{sample}_{prefix}.bam"
        sorted_bam   = out_dir / f"{sample}_{prefix}_sorted.bam"
        batch_path   = out_dir / f"{prefix}_igv.batch"
        session_path = out_dir / f"{prefix}_session.xml"
        igv_log      = out_dir / f"{prefix}_igv.log"

        logging.info("=== [%s] Input: %s → Snapshot: %s ===", sample, reads_path.name, png_path.name)

        # 1) minimap2 alignment (SAM)
        mm2 = ["minimap2", *mm_args]
        if threads:
            mm2 += ["-t", str(threads)]
        with open(sam_path, "w") as sam_out:
            run(mm2 + [str(ref_fa), str(reads_used)], stdout=sam_out)

        # 2) SAM -> BAM -> sort -> index
        run(["samtools", "view", "-b", str(sam_path), "-o", str(bam_path)])
        sort_cmd = ["samtools", "sort", "-o", str(sorted_bam), str(bam_path)]
        if threads:
            sort_cmd = ["samtools", "sort", "-@", str(threads), "-m", "256M", "-o", str(sorted_bam), str(bam_path)]
        run(sort_cmd)
        run(["samtools", "index", str(sorted_bam)])

        # 3) IGV session
        write_session_xml(ref_fa, sorted_bam, session_path)

        # 4) IGV pass 1: list tracks
        with open(batch_path, "w") as w:
            w.write(f'genome "{ref_fa.resolve()}"\n')
            w.write(f'load "{session_path.resolve()}"\n')
            w.write(f"sleep {sleep_ms}\n")
            w.write("list\n")
            w.write("exit\n")

        cmd1 = igv_cmdline(igv_cmd, use_xvfb, batch_path, ref_fa)
        with open(igv_log, "w") as lf:
            subprocess.run(cmd1, check=True, stdout=lf, stderr=lf, timeout=timeout_sec)

        bam_track = get_bam_track_name(igv_log)

        # 5) IGV pass 2: configure + snapshot
        with open(batch_path, "w") as w:
            w.write(f'genome "{ref_fa.resolve()}"\n')
            w.write(f'load "{session_path.resolve()}"\n')
            w.write(f"sleep {sleep_ms}\n")
            w.write("colorBy READ_STRAND\n")
            w.write("squish\n")
            w.write(f"squish {bam_track}\n")
            w.write("maxPanelHeight 800\n")
            w.write(f"goto {locus}\n")
            w.write("sort base\n")
            w.write(f'snapshotDirectory "{sample_dir.resolve()}"\n')  # IMPORTANT: snapshot in sample folder
            w.write(f'snapshot "{png_path.name}"\n')
            w.write("exit\n")

        cmd2 = igv_cmdline(igv_cmd, use_xvfb, batch_path, ref_fa)
        with open(igv_log, "a") as lf:
            subprocess.run(cmd2, check=True, stdout=lf, stderr=lf, timeout=timeout_sec)

        if png_path.exists():
            logging.info("[DONE] %s -> %s", reads_path.name, png_path)
        else:
            logging.warning("[WARN] Snapshot not created: %s", png_path)
            print_log_tail(igv_log)


# -------------------------
# CLI
# -------------------------
def main():
    parser = argparse.ArgumentParser(description="IGV snapshot visualizer for per-sample subfolders (YAML-driven).")
    parser.add_argument("input_folder", help="Parent folder containing sample subdirectories (e.g., ./sample)")
    parser.add_argument("-c", "--config", default="config.yaml", help="YAML config file path")
    args = parser.parse_args()

    setup_logging()
    cfg = load_config(Path(args.config))
    config_dir = Path(args.config).resolve().parent

    base = Path(args.input_folder).resolve()
    if not base.exists():
        logging.error("Input folder not found: %s", base)
        sys.exit(1)

    for sample in sorted(p for p in base.iterdir() if p.is_dir()):
        process_sample(sample, cfg, config_dir)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logging.error("Interrupted by user.")
        sys.exit(130)
