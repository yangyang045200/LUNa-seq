#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BC1/LUNa-seq 单文件静态 HTML 报告生成器（微调版）
- 左侧固定导航 + 右侧内容区
- 图片 Base64 内嵌，离线可读
- 自动发现 sample/ 下名称含 P 的样本目录（自然序）
- “I. Background” 与 “II. BC1 Results” 总览 + 每个样本七个小节
- 支持可选 PDF 导出（weasyprint）
依赖：pandas、numpy、matplotlib、openpyxl
"""
import argparse
import base64
import io
import os
import re
import sys
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------- 通用工具 ---------------------------

def natural_key(s: str):
    """自然序排序：P1 < P2 < P10"""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(s))]

def ensure_dir(p: Path):
    p.parent.mkdir(parents=True, exist_ok=True)

def guess_sep_lines(lines: List[str]) -> str:
    candidates = [",", "\t", ";", "|", " "]
    best = ","
    best_score = -1
    for sep in candidates:
        scores = [len([c for c in L.split(sep) if c != ""]) for L in lines if L.strip()]
        if not scores:
            continue
        score = sum(scores)
        if score > best_score:
            best_score = score
            best = sep
    return best

def read_text(p: Path, encoding="utf-8"):
    try:
        return p.read_text(encoding=encoding)
    except UnicodeDecodeError:
        return p.read_text(encoding="utf-8-sig", errors="ignore")

def img_to_base64(p: Path) -> Optional[str]:
    try:
        data = p.read_bytes()
        b64 = base64.b64encode(data).decode("ascii")
        ext = p.suffix.lower().lstrip(".")
        mime = {"png":"image/png","jpg":"image/jpeg","jpeg":"image/jpeg","svg":"image/svg+xml"}.get(ext, "image/png")
        return f"data:{mime};base64,{b64}"
    except Exception:
        return None

def render_table(df: pd.DataFrame, table_class="table", transpose_if_wide=True, index=True, format_mode: str = 'auto') -> str:
    """渲染 DataFrame 为 HTML 表格；列明显多于行且>3时转置并重命名列为 sample_1..
       format_mode='auto' 正常数值格式化；'raw' 保留文件中小数显示（不改动字符串）。
    """
    if transpose_if_wide and df.shape[1] > max(df.shape[0], 3):
        df = df.T.copy()
        df.columns = [f"sample_{i+1}" for i in range(df.shape[1])]
    if format_mode == 'raw':
        df_fmt = df.copy()
        return df_fmt.to_html(classes=table_class, escape=False, index=index, border=0)
    # 整数列不显示小数；浮点保留2位
    def fmt(v):
        try:
            if pd.isna(v):
                return ""
            fv = float(v)
            if fv.is_integer():
                return str(int(fv))
            return f"{fv:.2f}"
        except Exception:
            return str(v)
    df_fmt = df.copy()
    for c in df_fmt.columns:
        df_fmt[c] = df_fmt[c].map(fmt)
    return df_fmt.to_html(classes=table_class, escape=False, index=index, border=0)

def make_anchor_id(s: str) -> str:
    s = re.sub(r'\s+', '-', s.strip())
    s = re.sub(r'[^a-zA-Z0-9\-_]+', '', s)
    return s or "section"

def find_dirs_with_P(sample_root: Path) -> List[Path]:
    if not sample_root.exists():
        return []
    dirs = [d for d in sample_root.iterdir() if d.is_dir() and re.search(r'p', d.name, re.IGNORECASE)]
    return sorted(dirs, key=lambda d: natural_key(d.name))

def find_files_contains(root: Path, substr: str, exts: Optional[List[str]] = None, exact_name: bool=False) -> List[Path]:
    results = []
    if not root.exists():
        return results
    substr_lc = substr.lower()
    for p in sorted(root.rglob("*"), key=lambda x: natural_key(x.name)):
        if not p.is_file():
            continue
        if exts and p.suffix.lower() not in exts:
            continue
        name = p.name.lower()
        if exact_name:
            if p.name.lower() == substr_lc:
                results.append(p)
        else:
            if substr_lc in name:
                results.append(p)
    return results

def wrap_star_markup(text: str) -> str:
    """标记：*加粗*；**标红（不加粗）**"""
    def esc(s):
        return (s.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;"))
    text = esc(text)
    text = re.sub(r'\*\*([^*]+)\*\*', r'<span class="red">\1</span>', text)
    text = re.sub(r'\*([^*]+)\*', r'<b>\1</b>', text)
    return text

def plt_to_base64(plt_fig) -> str:
    buf = io.BytesIO()
    plt_fig.savefig(buf, format="png", bbox_inches="tight", dpi=160)
    plt.close(plt_fig)
    buf.seek(0)
    return f"data:image/png;base64,{base64.b64encode(buf.read()).decode('ascii')}"

# --------------------------- 图表绘制 ---------------------------

def plot_lines_from_csv(df: pd.DataFrame, title: Optional[str]=None) -> str:
    """第一列为 X，其余列为系列；跳过全空行；返回 base64 PNG"""
    if df is None or df.empty or df.shape[1] < 2:
        return ""
    df = df.dropna(how="all")
    x = df.columns[0]
    melted = df.melt(id_vars=[x], var_name="Series", value_name="Value")
    melted = melted.dropna(subset=["Value"])
    fig, ax = plt.subplots(figsize=(6, 3.2))
    for key, sub in melted.groupby("Series"):
        ax.plot(sub[x], sub["Value"], marker="o", label=str(key))
    ax.set_xlabel(str(x))
    ax.set_ylabel("Value")
    if title:
        ax.set_title(title)
    if melted["Series"].nunique() <= 10:
        ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    return plt_to_base64(fig)

def plot_mutation_donut(df_mtype: pd.DataFrame) -> str:
    """df_mtype 需包含: mismatch_type, proportion (0~1 或百分号), 可含 count"""
    if df_mtype is None or df_mtype.empty:
        return ""
    df = df_mtype.copy()
    if "proportion" in df.columns:
        def to_float(x):
            try:
                s = str(x).strip().rstrip("%")
                v = float(s)
                return v/100.0 if v > 1.0 else v
            except Exception:
                return np.nan
        df["prop"] = df["proportion"].map(to_float)
    else:
        total = float(df["count"].sum())
        df["prop"] = df["count"] / total if total > 0 else 0.0
    major = df[df["prop"] >= 0.04].copy()
    minor = df[df["prop"] < 0.04].copy()
    if not minor.empty:
        others = pd.DataFrame({
            "mismatch_type": ["Others"],
            "count": [minor["count"].sum() if "count" in minor.columns else 0],
            "prop": [minor["prop"].sum()]
        })
        major = pd.concat([major, others], ignore_index=True)
    major = major.sort_values("prop", ascending=False)
    labels = [str(x) for x in major["mismatch_type"].tolist()]
    sizes = major["prop"].tolist()
    cmap = matplotlib.cm.get_cmap("Blues")
    colors = [cmap(0.3 + 0.6*(i/max(1, len(sizes)-1))) for i in range(len(sizes))]
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    ax.pie(
        sizes,
        labels=[f"{l}" for l in labels],
        autopct=lambda pct: f"{pct:.1f}%",
        startangle=90,
        counterclock=False,
        colors=colors,
        wedgeprops=dict(width=0.35, edgecolor="white")
    )
    ax.axis('equal')
    fig.tight_layout()
    return plt_to_base64(fig)

# --------------------------- HTML 模板 ---------------------------

CSS = r"""
*{box-sizing:border-box}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,"Noto Sans","PingFang SC","Hiragino Sans GB","Microsoft YaHei","WenQuanYi Micro Hei",sans-serif;margin:0;color:#111;background:#fafafa}
a{color:#0b57d0;text-decoration:none}
a:hover{text-decoration:underline}
.sidebar{position:fixed;left:0;top:0;bottom:0;width:260px;background:#fff;border-right:1px solid #eee;overflow:auto;padding:18px}
.sidebar h2{font-size:16px;margin:8px 0 12px 0}
.sidebar .nav a{display:block;padding:6px 8px;border-radius:8px;margin:3px 0;color:#333}
.sidebar .nav a:hover{background:#f0f4ff}
.content{margin-left:260px;padding:24px 32px;max-width:1080px}
.center{max-width:960px;margin:0 auto}
.title{font-size:30px;font-weight:700;text-align:center;margin:10px 0 18px 0}
.subtitle{font-size:13px;color:#666;text-align:center;margin-bottom:18px}
.section-title{font-size:22px;margin:28px 0 12px 0;border-bottom:1px solid #eee;padding-bottom:6px}
.subsection-title{font-size:18px;margin:20px 0 10px 0;color:#333}
.table{border-collapse:collapse;width:100%;margin:6px 0 14px 0}
.table th,.table td{border:1px solid #e5e7eb;padding:6px 8px;text-align:left;vertical-align:top;font-size:13px}
.table th{background:#f9fafb}
.small-table .table th,.small-table .table td{font-size:12px}
.kv{display:grid;grid-template-columns:200px 1fr;gap:6px 12px}
.caption{font-size:12px;color:#555;margin:6px 0 18px 0;text-align:justify}
.img{max-width:70%;display:block;margin:4px auto 10px auto;border:1px solid #eee;border-radius:8px;background:white;box-shadow:0 1px 3px rgba(0,0,0,0.04)}
.row{display:flex;gap:14px;align-items:flex-start;flex-wrap:wrap}
.row .half{flex:1 1 420px}
.chips{display:flex;flex-wrap:wrap;gap:8px;margin:10px 0 6px 0}
.chips a{background:#eef2ff;padding:8px 12px;border-radius:999px;color:#1e293b;font-size:13px;border:1px solid #e5e7eb}
.narrow{max-width:70%;}
.code{background:#f6f8fa;border:1px solid #eee;border-radius:8px;padding:10px 12px;font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,"Liberation Mono","Courier New",monospace;font-size:12px;overflow:auto}
.note{background:#fff8e1;border:1px solid #ffecb3;padding:10px;border-radius:8px;color:#7a5900}
.bold{font-weight:700}
.red{color:#d93025}
img{image-rendering:auto}
hr{border:none;height:1px;background:#eee;margin:16px 0}
.small-table{display:inline-block;text-align:left;margin:0;}
.justified{text-align:justify;text-justify:inter-word;}
"""

HTML_HEAD = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
<style>{css}</style>
</head>
<body>
<div class="sidebar">
  <h2>Navigation</h2>
  <div class="nav">
    <a href="#home">Home</a>
    <a href="#background">I. Background</a>
    <a href="#bc1-results">II. BC1 Results</a>
    <div style="margin-top:10px;font-weight:600">Samples</div>
    {sample_nav}
  </div>
</div>
<div class="content">
<div class="center">
<a id="home"></a>
<div class="title">{title}</div>
<div class="subtitle">Single-file offline report · Base64-embedded images · Print-friendly</div>
"""

HTML_TAIL = r"""
</div></div>
</body>
</html>
"""

# --------------------------- 业务渲染 ---------------------------

BACKGROUND_TEXT = """I. Background
Lentiviral vectors (LVs) are widely deployed in gene and cell therapy manufacturing to stably integrate therapeutic payloads into mammalian genomes. Ensuring the genomic fidelity and biosafety of LVs is critical for both research and clinical applications. Conventional quality control (QC) strategies—such as infectivity assays, qPCR titers, and RCL testing—provide important functional metrics but may overlook low-frequency contaminants, structural variants, or stress-induced genomic lesions. To address these gaps, we developed *LUNa-seq*, a UMI-guided long-read sequencing workflow that consolidates multiple QC endpoints—including full-length integrity, point mutation burden, structural deletions, contamination, and vector–packaging recombination—into a single assay with molecule-level resolution.
"""

def render_background(process_image_path: Optional[Path]) -> str:
    parts = []
    parts.append('<a id="background"></a>')
    parts.append('<div class="section-title">I. Background</div>')
    parts.append('<div class="subsection-title">1.1 Background Content</div>')
    parts.append(f"<div class=\"justified\">{wrap_star_markup(BACKGROUND_TEXT)}</div>")
    parts.append('<div class="subsection-title">1.2 Workflow</div>')
    if process_image_path and process_image_path.exists():
        b64 = img_to_base64(process_image_path)
        if b64:
            parts.append(f'<img class="img" src="{b64}" alt="Workflow">')
        else:
            parts.append('<div class="note">Failed to load workflow image.</div>')
    else:
        parts.append('<div class="note">Workflow image not provided or not found.</div>')
    return "\n".join(parts)

def render_bc1_overview(sample_root: Path) -> Tuple[str, List[str]]:
    parts = []
    parts.append('<a id="bc1-results"></a>')
    parts.append('<div class="section-title">II. BC1 Results</div>')
    stats_csv = sample_root / "BC_Match_Statistics.csv"
    if stats_csv.exists():
        try:
            df = pd.read_csv(stats_csv)
            parts.append('<div class="subsection-title">BC_Match_Statistics</div>')
            parts.append(render_table(df, "table"))
        except Exception as e:
            parts.append(f'<div class="note">Failed to load BC_Match_Statistics.csv: {e}</div>')
    else:
        parts.append('<div class="note">未找到 sample/BC_Match_Statistics.csv</div>')
    sample_dirs = find_dirs_with_P(sample_root)
    parts.append('<div class="subsection-title">Samples</div>')
    parts.append('<div class="chips">')
    sample_names = []
    for d in sample_dirs:
        name = d.name
        sample_names.append(name)
        parts.append(f'<a href="#sample-{make_anchor_id(name)}">Sample: {name}</a>')
    parts.append('</div>')
    return "\n".join(parts), sample_names

def render_images_with_caption(img_paths: List[Path], caption: str) -> str:
    parts = []
    if not img_paths:
        parts.append('<div class="note">未找到相关图片。</div>')
        return "\n".join(parts)
    for p in img_paths:
        b64 = img_to_base64(p)
        if b64:
            parts.append(f'<img class="img" src="{b64}" alt="{p.name}">')
        else:
            parts.append(f'<div class="note">Failed to load image: {p.name}</div>')
    parts.append(f'<div class="caption">{wrap_star_markup(caption)}</div>')
    return "\n".join(parts)

def render_small_table(df: pd.DataFrame) -> str:
    return f'<div class="small-table narrow">{render_table(df, "table small", transpose_if_wide=False, format_mode="raw")}</div>'

def render_sample_block(sample_root: Path, sample_name: str) -> str:
    sdir = sample_root / sample_name
    sid = f"sample-{make_anchor_id(sample_name)}"
    out = []
    out.append(f'<a id="{sid}"></a>')
    out.append(f'<div class="section-title">Sample: {sample_name}</div>')
    # 1 Raw Data
    out.append('<div class="subsection-title">1. Raw Data</div>')
    imgs = []
    imgs += find_files_contains(sdir, f"{sample_name}_raw_data", exts=[".png"])
    if not imgs:
        imgs += find_files_contains(sdir, "raw_data", exts=[".png"])
    out.append(render_images_with_caption(
        imgs,
        "*The IGV overview displays all reads demultiplexed to the given sample. *The gray coverage bar at the top indicates sequencing depth. The pink track represents aligned matches. Base mismatches are color-coded by nucleotide: green = A, blue = C, orange = G, and yellow = T. Insertions are shown in purple, and horizontal dashes represent deletions."
    ))
    # 2 Raw Data Classification
    out.append('<div class="subsection-title">2. Raw Data Classification</div>')
    imgs = []
    imgs += find_files_contains(sdir, "fully_on_target_raw_data.png", exts=[".png"])
    if imgs:
        out.append(render_images_with_caption(imgs, "*The IGV view of fully on-target raw reads.*"))
    else:
        out.append('<div class="note">未找到 fully_on_target_raw_data 图。</div>')
    imgs = find_files_contains(sdir, "contaminant_or_recombination_judgement_raw_data", exts=[".png"])
    out.append(render_images_with_caption(imgs, "*The IGV view of raw data classification judgment.*"))
    ratios_csvs = find_files_contains(sdir, f"{sample_name}__ratios", exts=[".csv"]) or find_files_contains(sdir, f"{sample_name}_ratios", exts=[".csv"])
    if ratios_csvs:
        try:
            # 重要：以字符串读入，保留小数点后全部数字；渲染 raw 模式
            df = pd.read_csv(ratios_csvs[0], dtype=str)
            out.append(render_table(df, "table", format_mode='raw'))
            out.append('<div class="caption">{}</div>'.format(wrap_star_markup(
                "*Categories reported in the classification table include: (i) fully_on_target_proviral_reads* – Reads mapping contiguously to the lentiviral reference, with ≥90% of their length aligned or covering ≥95% of the reference sequence. These represent bona fide proviral molecules captured by the assay. *(ii) complete_coverage_within_on_target* – The subset of **fully on-target proviral reads **that span ≥95% of the entire reference length. This metric reflects the proportion of molecules reconstructed at near full-length resolution. Used as the full-length integrity measure for QC. *(iii) contaminant_or_recombinant_stage_n* – Residual reads classified iteratively by detecting dual motifs (5′ HIV-1 Ψ packaging signal and 3′ WPRE, ≤3 nt edit distance) in correct orientation. Each stage represents a new seed-derived contaminant/recombinant reference, with remaining reads remapped and assigned until no further dual-motif seeds remained.* (iv) non_specific_residual_products* – Reads not assignable to either the on-target proviral reference or to any contaminant/recombinant seed reference. *Except for the values marked in red, all other denominators are based on the total number of reads.*"
            )))
        except Exception as e:
            out.append(f'<div class="note">Failed to load ratios csv: {e}</div>')
    else:
        out.append('<div class="note">未找到 ratios CSV。</div>')
    # 3 UMI Group
    out.append('<div class="subsection-title">3. UMI Group</div>')
    record_xlsx = find_files_contains(sdir, "result_record", exts=[".xlsx"])
    if record_xlsx:
        try:
            df1 = pd.read_excel(record_xlsx[0], sheet_name=0, engine="openpyxl")
            out.append(render_small_table(df1.T))
            out.append('<div class="caption">{}</div>'.format(wrap_star_markup(
                "The table summarizes UMI diversity, filtering, and consensus support metrics. Unique_UMI: total number of distinct UMIs observed. Filtered_rate: proportion of all UMIs that are 18 nt in length and pass the expected motif check. Singleton_rate: fraction of UMIs supported by a single read only. Duplication_rate: proportion of reads identified as duplicates relative to unique UMIs. Redundant_frac: percentage of total reads collapsing into redundant UMI families. UMI≥3/5/10/20_count: number of UMI groups supported by at least 3, 5, 10, or 20 reads, respectively. UMI≥3/5/10/20_rate: proportion of UMI groups above each read-depth threshold. These metrics jointly reflect sequencing depth distribution across UMI groups and the robustness of consensus generation."
            )))
        except Exception as e:
            out.append(f'<div class="note">Failed to load result_record Sheet1: {e}</div>')
    else:
        out.append('<div class="note">未找到 result_record*.xlsx</div>')
    line_left = ""
    line_right = ""
    if record_xlsx:
        try:
            df2 = pd.read_excel(record_xlsx[0], sheet_name=1, engine="openpyxl")
            line_left = plot_lines_from_csv(df2, title="UMI length distribution")
        except Exception:
            pass
    statics = [p for p in sdir.glob("**/*") if p.is_file() and ("static" in p.name.lower()) and p.suffix.lower()==".csv"]
    if statics:
        try:
            df_static = pd.read_csv(statics[0])
            if df_static.shape[0] >= 2:
                df_static = df_static.iloc[1:].reset_index(drop=True)
            line_right = plot_lines_from_csv(df_static, title="Distribution of UMI family sizes")
        except Exception:
            pass
    if line_left or line_right:
        out.append('<div class="row">')
        if line_left:
            out.append(f'<div class="half"><img class="img" src="{line_left}"><div class="caption">{wrap_star_markup("*Figure A. UMI length distribution. *The plot shows the observed distribution of UMI lengths in the dataset. The majority of reads correspond to the intended target length (18 nt), with progressively fewer reads at shorter or longer lengths. This central enrichment with tapering at both ends follows the expected probability distribution and indicates that most UMIs were generated at the designed length. If a pronounced reduction is observed at the expected length, it may indicate that the defined UMI window is too narrow and does not sufficiently cover the ±10 nt flanking region, warranting re-evaluation of the UMI extraction parameters.")}</div></div>')
        if line_right:
            out.append(f'<div class="half"><img class="img" src="{line_right}"><div class="caption">{wrap_star_markup("*Figure B. Distribution of UMI family sizes.* The plot shows the frequency of UMI groups stratified by the number of supporting reads per UMI. Each point corresponds to the number of UMI families (y-axis) observed at a given read support threshold (x-axis). Singleton UMIs (x = 1) were excluded to improve visualization, as they dominated the distribution. The remaining data illustrate a monotonically decreasing trend, with progressively fewer UMI families supported by larger read counts.")}</div></div>')
        out.append('</div>')
    else:
        out.append('<div class="note">未生成 UMI 分布/族群规模曲线图。</div>')
    # 4 Consensus Sequence Extraction
    out.append('<div class="subsection-title">4. Consensus Sequence Extraction</div>')
    imgs = find_files_contains(sdir, f"{sample_name}_consensus", exts=[".png"])
    out.append(render_images_with_caption(imgs, "*The IGV view of all consensus reads.*"))
    # 5 Consensus Sequence Classification
    out.append('<div class="subsection-title">5. Consensus Sequence Classification</div>')
    imgs = find_files_contains(sdir, f"{sample_name}_fully_on_target_consensus", exts=[".png"])
    out.append(render_images_with_caption(imgs, "*The IGV view of fully on-target consensus reads.*"))
    imgs = find_files_contains(sdir, "contaminant_or_recombination_judgement_consensus", exts=[".png"])
    out.append(render_images_with_caption(imgs, "*The IGV view of Consensus classification judgment.* Reads retaining the conserved lentiviral sequences at both ends but containing large internal insertions are most often attributable to contamination events. Such patterns can be used for preliminary screening, and the putative contaminant references can then be aligned against a suspected lentiviral reference database for confirmation."))
    csvs = find_files_contains(sdir, f"{sample_name}__consensus_reads_ge_5_ratios", exts=[".csv"]) or find_files_contains(sdir, f"{sample_name}_consensus_reads_ge_5_ratios", exts=[".csv"])
    if csvs:
        try:
            # 重要：以字符串读入，保留小数点后全部数字；渲染 raw 模式
            df = pd.read_csv(csvs[0], dtype=str)
            out.append(render_table(df, "table", format_mode='raw'))
            out.append('<div class="caption">{}</div>'.format(wrap_star_markup(
                "Categories reported in the consensus classification table include:* (i) fully_on_target_proviral_reads – Reads mapping contiguously to the lentiviral reference, with ≥90% of their length aligned or covering ≥95% of the reference sequence. These represent bona fide proviral molecules captured by the assay. (ii) complete_coverage_within_on_target – The subset of **fully on-target proviral reads** that span ≥95% of the entire reference length. This metric reflects the proportion of molecules reconstructed at near full-length resolution, used as the full-length integrity measure for QC. (iii) contaminant_or_recombinant_stage_n – Residual reads classified iteratively by detecting dual motifs (5′ HIV-1 Ψ packaging signal and 3′ WPRE, ≤3 nt edit distance) in correct orientation. Each stage represents a new seed-derived contaminant/recombinant reference, with remaining reads remapped and assigned until no further dual-motif seeds remained. (iv) non_specific_residual_products – Reads not assignable to either the on-target proviral reference or to any contaminant/recombinant seed reference. *Except for the values marked in red, all other denominators are based on the total number of reads.*"
            )))
        except Exception as e:
            out.append(f'<div class="note">Failed to load consensus ratios csv: {e}</div>')
    else:
        out.append('<div class="note">未找到 consensus ratios CSV。</div>')
    # 6 Mutation Analysis
    out.append('<div class="subsection-title">6. Mutation Analysis</div>')
    mut_csvs = find_files_contains(sdir, "mutation", exts=[".csv"])
    if mut_csvs:
        mut_path = mut_csvs[0]
        try:
            raw = read_text(mut_path)
            lines = [L for L in raw.splitlines() if L is not None]
            # 先展示前两行数据作为小表
            header_guess_lines = [ln for ln in lines[:4] if ln.strip()]
            sep = guess_sep_lines(header_guess_lines)
            df_head = pd.DataFrame([
                [c.strip() for c in header_guess_lines[0].split(sep)] if len(header_guess_lines)>=1 else [],
                [c.strip() for c in header_guess_lines[1].split(sep)] if len(header_guess_lines)>=2 else []
            ])
            df_head.columns = [f"col_{i+1}" for i in range(df_head.shape[1])]
            out.append(render_small_table(df_head))
            out.append('<div class="caption">{}</div>'.format(wrap_star_markup("*Base-level mutation profile of fully on-target consensus reads (≥5 UMI support)*")))
            # 识别以 "# full_length mismatch_type count proportion" 为表头的段（把 full_length mismatch_type 视为一个整体）
            start_idx = None
            for i, ln in enumerate(lines):
                norm = ln.strip()
                if re.search(r'^\s*#\s*full_length(?:\s+|_)mismatch_type\b', norm, re.IGNORECASE):
                    start_idx = i
                    break
            if start_idx is not None:
                hdr_line = lines[start_idx].lstrip("#").strip()
                hdr_line = hdr_line.replace("、", " ")
                hdr_line = re.sub(r'\s+', ' ', hdr_line)
                hdr_protected = re.sub(r'full_length\s+mismatch_type', 'full_length_mismatch_type', hdr_line, flags=re.IGNORECASE)
                header = re.split(r'[\s,]+', hdr_protected)
                header = [h.replace('full_length_mismatch_type','mismatch_type') for h in header]
                # 收集数据行
                rows = []
                for j in range(start_idx+1, len(lines)):
                    if lines[j].strip().startswith("#") or not lines[j].strip():
                        break
                    cells = re.split(r'[\s,]+', lines[j].strip().replace('、',' '))
                    if len(cells) >= len(header):
                        rows.append(cells[:len(header)])
                if rows:
                    df_m = pd.DataFrame(rows, columns=header)
                    # 规范列名
                    lower = [c.lower() for c in df_m.columns]
                    # 将第一列（full_length mismatch_type）规范为 mismatch_type
                    if "mismatch_type" not in lower:
                        df_m = df_m.rename(columns={df_m.columns[0]: "mismatch_type"})
                    # count / proportion 容错
                    if "count" not in lower:
                        if "counts" in lower:
                            df_m = df_m.rename(columns={df_m.columns[lower.index("counts")]: "count"})
                    if "proportion" not in lower and "ratio" in lower:
                        df_m = df_m.rename(columns={df_m.columns[lower.index("ratio")]: "proportion"})
                    # donut 绘制
                    donut_b64 = plot_mutation_donut(df_m)
                    # 左表右图
                    out.append('<div class="subsection-title">*Distribution of nucleotide substitution types (%)*</div>')
                    out.append('<div class="row">')
                    out.append(f'<div class="half">{render_table(df_m, "table", transpose_if_wide=False, format_mode="raw")}</div>')
                    if donut_b64:
                        out.append(f'<div class="half"><img class="img" src="{donut_b64}" alt="mutation donut"></div>')
                    out.append('</div>')
                else:
                    out.append('<div class="note">未解析到 mismatch_type 表格数据。</div>')
            else:
                # 未匹配到标记行，整体按普通 CSV 渲染
                try:
                    df_mut = pd.read_csv(mut_path)
                except Exception:
                    df_mut = pd.read_csv(mut_path, sep=None, engine="python")
                out.append(render_table(df_mut, "table"))
            # 展示含 plot 的 PNG
            plot_pngs = [p for p in sdir.glob("**/*.png") if "plot" in p.name.lower()]
            if plot_pngs:
                out.append(render_images_with_caption(plot_pngs, "*Genomic distribution of variant events across the lentiviral vector*. Variant event rates were calculated in 50 bp sliding windows (edges ±5 bp ignored) and plotted across the vector genome. Mismatches (blue), deletions (orange), and insertions (green) are shown as stacked bars."))
        except Exception as e:
            out.append(f'<div class="note">解析 Mutation CSV 失败：{e}</div>')
    else:
        out.append('<div class="note">未找到 mutation CSV。</div>')
    # 7 RCL Check
    out.append('<div class="subsection-title">7. RCL Check</div>')
    imgs = find_files_contains(sdir, "RCL_from_cr", exts=[".png"])
    out.append(render_images_with_caption(imgs, "*Recombinant fragments containing RCL, derived from contaminant_or_recombinant_consensus*"))
    imgs = find_files_contains(sdir, "RCL_from_non_specific", exts=[".png"])
    out.append(render_images_with_caption(imgs, "*Recombinant fragments containing RCL, derived from non_specific_residual_products_consensus*"))
    tsvs = find_files_contains(sdir, "RCL_summary", exts=[".tsv"])
    if tsvs:
        try:
            df = pd.read_csv(tsvs[0], sep="\t")
            out.append(render_table(df, "table"))
            out.append('<div class="caption">{}</div>'.format(wrap_star_markup("*Summary of RCL fragment detection*")))
        except Exception as e:
            out.append(f'<div class="note">Failed to load RCL_summary.tsv: {e}</div>')
    else:
        out.append('<div class="note">未找到 RCL_summary.tsv。</div>')
    return "\n".join(out)

# --------------------------- 主流程 ---------------------------

def build_html(data_root: Path, title: str, process_image: Optional[Path]) -> str:
    sample_root = data_root / "sample"
    sample_dirs = find_dirs_with_P(sample_root)
    sample_nav = []
    for d in sample_dirs:
        name = d.name
        sample_nav.append(f'<a href="#sample-{make_anchor_id(name)}">{name}</a>')
    html_parts = [HTML_HEAD.format(title=title, css=CSS, sample_nav="\n".join(sample_nav))]
    html_parts.append(render_background(process_image))
    html_parts.append("<hr>")
    overview_html, sample_names = render_bc1_overview(sample_root)
    html_parts.append(overview_html)
    html_parts.append("<hr>")
    for name in sample_names:
        html_parts.append(render_sample_block(sample_root, name))
        html_parts.append("<hr>")
    html_parts.append(HTML_TAIL)
    return "\n".join(html_parts)


def main():
    ap = argparse.ArgumentParser(description="BC1/LUNa-seq 单文件静态 HTML 报告生成器（微调版）")
    ap.add_argument("--data-root", required=True, help="项目根目录（包含 sample/）")
    ap.add_argument("--output", default="BC1_report.html", help="输出 HTML 路径（默认 BC1_report.html）")
    ap.add_argument("--title", default="慢病毒", help="报告主标题")
    ap.add_argument("--process-image", default=None, help="流程图图片路径（可选）")
    ap.add_argument("--pdf", default=None, help="可选，输出 PDF 路径（需安装 weasyprint）")
    args = ap.parse_args()

    data_root = Path(args.data_root).resolve()
    if not data_root.exists():
        print(f"[ERROR] data-root 不存在：{data_root}", file=sys.stderr)
        sys.exit(2)

    process_image = Path(args.process_image).resolve() if args.process_image else None
    html_str = build_html(data_root, args.title, process_image)

    out_html = Path(args.output).resolve()
    ensure_dir(out_html)
    out_html.write_text(html_str, encoding="utf-8")
    print(f"[OK] Wrote: {out_html}")



if __name__ == "__main__":
    main()
