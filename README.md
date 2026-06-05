# ST-CCC Benchmark
<p align="left">

  <img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R">

  <img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">

  <img src="https://img.shields.io/badge/COMMOT-Spatial%20Communication-9cf?style=for-the-badge" alt="COMMOT">

  <img src="https://img.shields.io/badge/spaCI-Spatial%20CCI-ff69b4?style=for-the-badge" alt="spaCI">

  <img src="https://img.shields.io/badge/Spacia-Spatial%20Analysis-ffb347?style=for-the-badge" alt="Spacia">

  <img src="https://img.shields.io/badge/VGAE--CCI-Graph%20Neural%20Network-lightgrey?style=for-the-badge" alt="VGAE-CCI">

</p>

**Spatial transcriptomics (ST)** technologies have transformed our understanding of gene expression in tissues by preserving spatial context. However, understanding how cells communicate with each other—especially through ligand-receptor interactions—in these spatial environments remains a computational challenge.

While many tools infer **cell–cell communication (CCC)** from single-cell RNA-seq (scRNA-seq), these approaches often ignore spatial proximity, leading to biologically incomplete results. To address this, a new wave of tools has emerged that work natively on spatial transcriptomics data.

This benchmark systematically evaluates those native **Spatial Transcriptomics** tools that perform CCC inference directly from spatial transcriptomic datasets—without relying on external single-cell references.


## 📌 Objective

The benchmark evaluates:

 - Cell-type communication patterns.
 - Number of inferred ligand-receptor interactions.
 - Overlap of predicted ligand-receptor pairs between methods.
 - Scalability across datasets of different sizes.

## 🧰 Tools Included

| Tool           | ST-Native | Methodology | Language | GitHub |
|----------------|-----------|------------|----------|--------|
| COMMOT           | ✅     | Optimal Transport  | Python   | [COMMOT](https://github.com/zcang/COMMOT) |
| spaCI            | ✅     | Graph Neural Network | Python   | [spaCI](https://github.com/QSong-github/spaCI) |
| Spacia           | ✅     | Bayesian Inference (MCMC) | Python   | [Spacia](https://github.com/yunguan-wang/Spacia/blob/main/README.md) |
| VGAE-CCI         | ✅     | Variational Graph Autoencoder | Python   | [VGAE-CCI](https://github.com/zhangxiangz/VGAECCI) |


## 🧬 Dataset 
### MERFISH Datasets
 - Conoral mouse brain
 - S1R1 mouse brain
   
 Source:
 - https://www.nature.com/articles/s41586-024-08334-8#data-availability
 - https://info.vizgen.com/mouse-brain-map

### MOSTA Dataset
 - Mouse brain
   
 Source:
 - https://db.cngb.org/stomics/mosta/download/

## Installation
Each tool was executed in an independent Conda environment to ensure reproducibility and avoid dependency conflicts.

### Usage Examples

```bash
conda env create -f environments/ST_CCC_COMMOT.yml
conda activate ST_CCC_COMMOT
```
Environment files are available in the repository.

## Spatial Transcriptomics Format Converter

A Python utility for converting spatial transcriptomics datasets between common formats (CSV, H5AD, TXT, and SpaCI-compatible CSV) with optional CPM + log1p normalization for raw counts. This tool ensures proper handling of spatial coordinates and metadata, making datasets ready for downstream analysis.


### Arguments
 ```--from ```: Input format (csv or h5ad)

 ```--to ```: Output format (txt, h5ad, csv, csv_spaCI)

 ```--expr ```: Path to gene expression CSV (required for CSV input)

``` --meta ```: Path to metadata CSV (required for CSV input)

 ```--h5ad ```: Path to H5AD file (required for H5AD input)

 ```--out ```: Output directory (default: ./output)

 ```--name``` : Sample name prefix (default: sample)

 ```--normalized ```: Apply CPM + log1p normalization (only for TXT output; if omitted, assumes data is already normalized)
 
### Usage Examples

### Convert a CSV dataset to TXT with normalization:

```bash
python Preprocessing_formats.py --from csv --to txt --expr counts.csv --meta meta.csv --out output_dir --name sample_name --normalized 
```

### Convert an H5AD file to SpaCI CSV without normalization:

```bash
python Preprocessing_formats.py --from h5ad --to csv_spaCI --h5ad sample.h5ad --out output_dir --name sample_name 
```

## Citation 

If you use this repository, please cite:
Akram Abushmais.

Comparative Analysis of Spatial Transcriptomics Approaches for Cell–Cell Communication.


