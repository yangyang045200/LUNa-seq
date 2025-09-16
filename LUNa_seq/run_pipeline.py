#!/usr/bin/env python3
"""
Master pipeline orchestrator
- raw_data_classification.py is OPTIONAL (enable via YAML flag)
- consensus_classification.py is REQUIRED
- IGV visualization is OPTIONAL (visualization_igv.enabled)
- RCL test is OPTIONAL (rcl_test.enabled)
- Checkpoint/resume via .pipeline_state/*.done
- Per-step wall time written to .pipeline_state/timings.csv
"""

import argparse
import csv
import logging
import subprocess
import sys
import time
from pathlib import Path
import yaml

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

def load_config(path: Path):
    if not path.exists():
        logging.error("Config file %s not found", path)
        sys.exit(1)
    try:
        with open(path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    except Exception as e:
        logging.error("Failed to parse config: %s", e)
        sys.exit(1)

def resolve_script(name: str, base_dir: Path) -> Path:
    """
    Resolve a script by name. Search order:
      1) Same directory as this orchestrator
      2) ./scripts/<name>
    """
    here = base_dir
    p1 = here / name
    if p1.exists():
        return p1
    p2 = here / "scripts" / name
    if p2.exists():
        return p2
    return Path()  # not found

def run_step(marker_dir: Path, step_id: str, step_name: str, script_path: Path,
             args_list: list[str], resume: bool, timings_csv: Path):
    done_file = marker_dir / f"{step_id}_{step_name}.done"
    if resume and done_file.exists():
        logging.info("Skipping %s (already completed)", step_name)
        return

    if not script_path or not script_path.exists():
        logging.error("Script not found for step %s: %s", step_name, script_path)
        sys.exit(1)

    cmd = [sys.executable, str(script_path)] + args_list
    logging.info(">>> %s", " ".join(cmd))

    t0 = time.time()
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error("Step %s failed with exit code %d", step_name, e.returncode)
        sys.exit(e.returncode)
    t1 = time.time()

    done_file.parent.mkdir(parents=True, exist_ok=True)
    done_file.write_text("")

    timings_csv.parent.mkdir(parents=True, exist_ok=True)
    write_header = not timings_csv.exists()
    with open(timings_csv, "a", newline="") as fh:
        w = csv.writer(fh)
        if write_header:
            w.writerow(["step_id", "step_name", "script", "start_ts", "end_ts", "elapsed_sec"])
        w.writerow([step_id, step_name, script_path.name, int(t0), int(t1), round(t1 - t0, 3)])

def main():
    parser = argparse.ArgumentParser(description="Run full pipeline via YAML config")
    parser.add_argument("-c", "--config", required=True, type=Path, help="Path to config.yaml")
    parser.add_argument("--resume", action="store_true", help="Resume from last checkpoint")
    args = parser.parse_args()

    setup_logging()
    cfg = load_config(args.config)

    base_dir = Path(__file__).resolve().parent
    state_dir = base_dir / ".pipeline_state"
    timings_csv = state_dir / "timings.csv"

    if not args.resume and state_dir.exists():
        for f in state_dir.iterdir():
            try:
                f.unlink()
            except Exception:
                pass
    state_dir.mkdir(exist_ok=True)

    demux_cfg = cfg.get("demultiplex") or cfg
    sample_root = Path(demux_cfg.get("output_dir", "sample"))

    # Resolve scripts
    demultiplex_py              = resolve_script("demultiplex.py", base_dir)
    raw_data_classification_py  = resolve_script("raw_data_classification.py", base_dir)
    umi_grouped_py              = resolve_script("UMI_grouped.py", base_dir)
    consensus_abPOA_py          = resolve_script("consensus_abPOA.py", base_dir)
    consensus_classification_py = resolve_script("consensus_classification.py", base_dir)
    variation_analysis_py       = resolve_script("variation_analysis.py", base_dir)
    visualize_igv_py            = resolve_script("visualize_igv_snapshots.py", base_dir)
    rcl_test_py                 = resolve_script("rcl_test.py", base_dir)

    if not consensus_classification_py or not consensus_classification_py.exists():
        logging.error("consensus_classification.py is REQUIRED but was not found.")
        sys.exit(1)

    rdc_cfg = cfg.get("raw_data_classification", {})
    run_raw_classification = bool(rdc_cfg.get("enabled", False)) and raw_data_classification_py.exists()

    viz_cfg = cfg.get("visualization_igv", {})
    run_visualization = bool(viz_cfg.get("enabled", False)) and visualize_igv_py.exists()

    rcl_cfg = cfg.get("rcl_test", {})
    run_rcl = bool(rcl_cfg.get("enabled", False)) and rcl_test_py.exists()

    # Build steps
    steps: list[tuple[str, str, Path, list[str]]] = []

    # 01 Demultiplex
    steps.append(("01", "demultiplex", demultiplex_py, ["-c", str(args.config)]))

    # 02 (optional) raw data classification
    if run_raw_classification:
        steps.append(("02", "raw_data_classification", raw_data_classification_py, [str(sample_root), "-c", str(args.config)]))
        next_id = "03"
    else:
        logging.info("raw_data_classification: disabled (set raw_data_classification.enabled: true to enable).")
        next_id = "02"

    # Next steps
    steps.append((next_id, "UMI_grouped", umi_grouped_py, ["-c", str(args.config)]))
    next_id = f"{int(next_id) + 1:02d}"

    steps.append((next_id, "consensus_abPOA", consensus_abPOA_py, ["-c", str(args.config)]))
    next_id = f"{int(next_id) + 1:02d}"

    steps.append((next_id, "consensus_classification", consensus_classification_py, [str(sample_root), "-c", str(args.config)]))
    next_id = f"{int(next_id) + 1:02d}"

    steps.append((next_id, "variation_analysis", variation_analysis_py, ["-c", str(args.config)]))
    next_id = f"{int(next_id) + 1:02d}"

    # OPTIONAL IGV visualization
    if run_visualization:
        steps.append((next_id, "igv_visualization", visualize_igv_py, [str(sample_root), "-c", str(args.config)]))
        next_id = f"{int(next_id) + 1:02d}"
    else:
        logging.info("igv_visualization: disabled (set visualization_igv.enabled: true to enable).")

    # OPTIONAL RCL test (placed after visualization)
    if run_rcl:
        steps.append((next_id, "rcl_test", rcl_test_py, [str(sample_root), "-c", str(args.config)]))
        next_id = f"{int(next_id) + 1:02d}"
    else:
        logging.info("rcl_test: disabled (set rcl_test.enabled: true to enable).")

    # Execute
    for step_id, name, script, argv in steps:
        run_step(state_dir, step_id, name, script, argv, args.resume, timings_csv)

    logging.info("🎉 All steps completed. Timings → %s", timings_csv)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logging.error("Interrupted by user.")
        sys.exit(130)
