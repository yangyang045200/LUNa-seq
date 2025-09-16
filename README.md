# LUNa-seq
A reproducible end-to-end pipeline to demultiplex, UMI-tag, build per-UMI consensuses, classify on-/off-target reads, visualize in IGV, RCL-test, and summarize results.

1) Folder layout
project_root/
├── config.yaml                    # Global pipeline config (edit this)
├── run_pipeline.py                # Orchestrator (supports --resume)
├── scripts/                       # All sub-steps
│   ├── demultiplex.py
│   ├── raw_data_classification.py         (optional)
│   ├── UMI_grouped.py
│   ├── consensus_abPOA.py
│   ├── consensus_classification.py  (required)
│   ├── variation_analysis.py
│   └── rcl_test.py
│   └── visualize_igv_snapshots.py
├── sample/                        # Per-sample outputs (created by pipeline)
│   └── <SAMPLE_ID>/ …
├── visualization/                 # Output visualization file (bam; sam...) by sample
├── ref/                           # Put per-sample references here as <SAMPLE_ID>.fasta
├── Demultiplexed.csv              # Sample sheet: sample ID ↔ barcode
├── raw.fastq.gz                   # Your raw input (can be .fastq or .fastq.gz)
├── vendor/                        # Put IGV here (manual install)
├── process.png                    # Flowchart for the interactive page
└── show.py                       # Interactive summary UI (optional)
The pipeline will copy <SAMPLE_ID>.fasta from ref/ into each sample/<SAMPLE_ID>/ as needed.

2) Create the conda environment (recommended)
Use the classic solver and the conda-forge/bioconda channels.
# 2.1 Create environment
conda create -y -n luna-seq -c conda-forge -c bioconda \
  python=3.11 minimap2=2.28 samtools=1.20 abpoa=1.5.4 openjdk=17 \
  biopython=1.83 pysam='0.23.*' regex pyyaml=6.0.2 matplotlib=3.10.1 \
  pandas=2.2.3 tqdm=4.* openpyxl=3.1.* seaborn=0.13.*
# 2.2 Activate
conda activate luna-seq

Linux (Ubuntu/WSL): IGV headless snapshots require Xvfb and fonts:
sudo apt-get update
sudo apt-get install -y xvfb fonts-dejavu-core

3) Install IGV (manual)
IGV is not available as a conda package at the pinned version; download a 2.19.x build and set its path.

# Example: place IGV under vendor/
mkdir -p vendor
cd vendor
# Download IGV_2.19.6 and unpack so you have vendor/igv-2.19.6/igv.sh
# (on Linux, a .zip can be unzipped; on macOS, use the app bundle or igv.sh in the distribution)

# Back to project root
cd ..

# Export IGV location in your active environment
export IGV_CMD="$PWD/vendor/igv-2.19.6/igv.sh"
chmod +x "$IGV_CMD"

# Persist it for future activations of luna-seq
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
printf 'export IGV_CMD=%q\n' "$IGV_CMD" > "$CONDA_PREFIX/etc/conda/activate.d/igv_cmd.sh"

4) Prepare your inputs

Raw reads
Put your input as raw.fastq.gz (or raw.fastq) at the project root.

Sample sheet (Demultiplexed.csv)
Minimal example (no header or with header; your demultiplex script dictates the exact format—keep the same you already use):
ID	BC
LV_P180_1	aggtgtct
LV_P2009_1_H	ccgcatgt
LV_P1946_1	ccttcgga

References (ref/)
Provide one FASTA per sample named exactly <SAMPLE_ID>.fasta.
Example: ref/LV_P180_2.fasta

5) Configure the pipeline
Open config.yaml and adjust:
Input file name (raw.fastq.gz)
Demultiplexing options (barcodes, output dir sample/)
Threads (e.g., THREADS: 8)
Any paths relevant to your environment
Only change documented keys; the internal logic and file naming conventions are preserved.

6) Run
conda activate luna-seq
python run_pipeline.py -c config.yaml
Resume from the last completed step:
python run_pipeline.py -c config.yaml --resume
Progress/state files are written under .pipeline_state/.

7) Outputs

Per sample (sample/<SAMPLE_ID>/):
Consensus FASTA/TSV/CSV and ratios
Contaminant/recombinant stage outputs
RCL/RCR integration scans and plots (e.g., RCL_from_cr.png, RCL_from_non_specific.png)
IGV snapshots:
visualization/<SAMPLE_ID>/*.png (e.g., <ID>_raw_data.png, contaminant_or_recombination_judgement_consensus.png, etc.)
Logs & timings:
.pipeline_state/timings.csv

8) Interactive summary (optional)
If you maintain an interactive report:
↳python inter.py
Prepare process.png and any assets the UI expects.

9) Package not found / ModuleNotFoundError
Install the missing package via conda and add it to your environment instructions:
conda install -y -c conda-forge <package>

10) Citation
If you use LUNa-seq in your work, please cite the accompanying manuscript (add when available).
