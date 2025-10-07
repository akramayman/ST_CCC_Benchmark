# ST-CCC Benchmark

**Spatial transcriptomics (ST)** technologies have transformed our understanding of gene expression in tissues by preserving spatial context. However, understanding how cells communicate with each other—especially through ligand-receptor interactions—in these spatial environments remains a computational challenge.

While many tools infer **cell–cell communication (CCC)** from single-cell RNA-seq (scRNA-seq), these approaches often ignore spatial proximity, leading to biologically incomplete results. To address this, a new wave of tools has emerged that work natively on spatial transcriptomics data.

This benchmark systematically evaluates those native **Spatial Transcriptomis** tools that perform CCC inference directly from spatial transcriptomic datasets—without relying on external single-cell references.


## 📌 Objective

To compare and evaluate tools that perform CCC inference directly on native simulation data and two kind of spatial transcriptomics (ST) Imageing (Vizgen) and Sequenceing data (MOSTA).

## 🧰 Tools Included

| Tool           | ST-Native | Language | GitHub |
|----------------|-----------|----------|--------|
| COMMOT           | ✅        | Python   | [COMMOT](https://github.com/zcang/COMMOT) |
| spaCI            | ✅        | Python   | [spaCI](https://github.com/QSong-github/spaCI) |
| Spacia           | ✅        | Python   | [Spacia](https://github.com/yunguan-wang/Spacia/blob/main/README.md) |
| VGAE-CCI         | ✅        | Python   | [VGAE-CCI](https://github.com/zhangxiangz/VGAECCI) |


## Spatial Transcriptomics Format Converter

A Python utility for converting spatial transcriptomics datasets between common formats (CSV, H5AD, TXT, and SpaCI-compatible CSV) with optional CPM + log1p normalization for raw counts. This tool ensures proper handling of spatial coordinates and metadata, making datasets ready for downstream analysis.


### Arguments
 --from : Input format (csv or h5ad)

 --to : Output format (txt, h5ad, csv, csv_spaCI)

 --expr : Path to gene expression CSV (required for CSV input)

 --meta : Path to metadata CSV (required for CSV input)

 --h5ad : Path to H5AD file (required for H5AD input)

 --out : Output directory (default: ./output)

 --name : Sample name prefix (default: sample)

 --normalized : Apply CPM + log1p normalization (only for TXT output; if omitted, assumes data is already normalized)
 
### Usage Examples

### Convert a CSV dataset to TXT with normalization:

```bash
python Preprocessing_formats.py --from csv --to txt --expr counts.csv --meta meta.csv --out output_dir --name sample_name --normalized 
```

### Convert an H5AD file to SpaCI CSV without normalization:

```bash
python Preprocessing_formats.py --from h5ad --to csv_spaCI --h5ad sample.h5ad --out output_dir --name sample_name 
```



### Copyright (c) 2025 Akram Abushmais

