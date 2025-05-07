# ST-CCC Benchmark

**Spatial transcriptomics (ST)** technologies have transformed our understanding of gene expression in tissues by preserving spatial context. However, understanding how cells communicate with each other—especially through ligand-receptor interactions—in these spatial environments remains a computational challenge.

While many tools infer **cell–cell communication (CCC)** from single-cell RNA-seq (scRNA-seq), these approaches often ignore spatial proximity, leading to biologically incomplete results. To address this, a new wave of tools has emerged that work natively on spatial transcriptomics data.

This benchmark systematically evaluates those native **Spatial Transcriptomis** tools that perform CCC inference directly from spatial transcriptomic datasets—without relying on external single-cell references.


## 📌 Objective

To compare and evaluate tools that perform CCC inference directly on native spatial transcriptomics (ST) data.

## 🧰 Tools Included

| Tool           | ST-Native | Language | GitHub |
|----------------|-----------|----------|--------|
| COMMOT         | ✅        | Python   | [COMMOT](https://github.com/zcang/COMMOT) |
| SpatialScope   | ✅        | Python   | [SpatialScope](https://github.com/YangLabHKUST/SpatialScope) |
| spaCI          | ✅        | Python   | [spaCI](https://github.com/QSong-github/spaCI) |
| Spacia         | ✅        | Python   | [Spacia](https://github.com/yunguan-wang/Spacia/blob/main/README.md) |
| NCEM           | ✅        | Python   | [NCEM](https://github.com/theislab/ncem) |
| CellCharter    | ✅        | Python   | [CellCharter](https://github.com/CSOgroup/cellcharter) |
| stLearn        | ✅        | Python   | [stLearn](https://github.com/BiomedicalMachineLearning/stLearn) |
| MISTy          | ✅        | R        | [MISTy](https://github.com/saezlab/MISTy) |


## 📦 Installation

Each tool has its own Conda environment.

- Detailed installation instructions for each tool are in Tools.
  
    -Environment (.yml) file
  
    -Instructions INSTALL.md
  
Copyright (c) 2025 Akram Abushmais

