# 🐄 RNA Stability Dynamics in Bovine Mastitis

[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs12864--025--12396--x-blue)](https://link.springer.com/article/10.1186/s12864-025-12396-x)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.1.2+-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.9+-green.svg)](https://www.python.org/)

> **A novel computational framework integrating RNA stability dynamics with gene co-expression networks to identify key regulatory modules in bovine mastitis**

## 📖 Overview

This repository contains the complete analysis pipeline for our published study in **BMC Genomics** (2025):

**"RNA stability: a novel perspective on gene regulatory networks in bovine mastitis"**

Mastitis is the most economically important infectious disease in dairy cattle. While gene expression profiling has been extensively studied, **RNA stability dynamics** - a critical post-transcriptional regulatory mechanism - remained unexplored in this disease context.

### Key Contributions

- 🧬 **Novel approach**: First study to investigate RNA stability patterns in bovine mastitis
- 📊 **WGCNA integration**: Constructed weighted gene co-expression networks based on stability profiles
- 🔗 **Network propagation**: Developed a novel algorithm to validate module functional connectivity
- 🎯 **Candidate genes**: Identified both known (IL6, TLR2, NFKB1) and novel (RELB, CACTIN, DHX9) mastitis-associated genes

## 🔬 Methodology

```
RNA-Seq Data (GSE51856)
         ↓
┌─────────────────────────────────────┐
│   REMBRANDTS: RNA Stability         │
│   (Exon/Intron ratio analysis)      │
└─────────────────────────────────────┘
         ↓
┌─────────────────────────────────────┐
│   WGCNA: Module Detection           │
│   (Stability-based clustering)      │
└─────────────────────────────────────┘
         ↓
┌─────────────────────────────────────┐
│   Network Propagation Validation    │
│   (NetColoc on PCNet interactome)   │
└─────────────────────────────────────┘
         ↓
┌─────────────────────────────────────┐
│   Functional Enrichment Analysis    │
│   (GO, KEGG pathways)               │
└─────────────────────────────────────┘
```

## 📁 Repository Structure

```
bovine-mastitis-rna-stability/
├── scripts/
│   ├── R/
│   │   ├── 01_rna_stability_calculation.R    # REMBRANDTS preprocessing
│   │   ├── 02_wgcna_analysis.R               # WGCNA module detection
│   │   └── 03_functional_enrichment.R        # GO/KEGG analysis
│   ├── python/
│   │   └── network_propagation.py            # NetColoc implementation
│   └── preprocessing/
│       └── download_and_align.sh             # Raw data processing
├── notebooks/
│   └── Network_Propagation_WGCNA_Mastitis.ipynb
├── data/
│   ├── raw/                                  # Input stability matrices
│   ├── processed/                            # Normalized data
│   └── results/                              # Module assignments, enrichment
├── figures/                                  # Publication-ready visualizations
├── docs/                                     # Additional documentation
├── environment.yml                           # Conda environment
├── requirements.txt                          # Python dependencies
└── README.md
```

## 🚀 Quick Start

### Prerequisites

- R ≥ 4.1.2
- Python ≥ 3.9
- Conda (recommended)

### Installation

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/bovine-mastitis-rna-stability.git
cd bovine-mastitis-rna-stability

# Create conda environment
conda env create -f environment.yml
conda activate mastitis-stability

# Install R packages
Rscript scripts/R/install_packages.R
```

### Running the Analysis

```bash
# Step 1: Calculate RNA stability (if starting from raw counts)
Rscript scripts/R/01_rna_stability_calculation.R

# Step 2: WGCNA module detection
Rscript scripts/R/02_wgcna_analysis.R

# Step 3: Network propagation validation
python scripts/python/network_propagation.py

# Step 4: Functional enrichment
Rscript scripts/R/03_functional_enrichment.R
```

## 📊 Key Results

### Module Detection
| Module | Genes | Stability Pattern | Key Function |
|--------|-------|-------------------|--------------|
| 🔴 Red | 127 | ↑ Increased post-infection | Innate immunity |
| 🟡 Yellow | 549 | ↑ Increased post-infection | Cytokine signaling |
| 🔵 Blue | 892 | ↓ Decreased post-infection | Homeostasis |

### Top Hub Genes (Yellow Module)
- **IL6** - Interleukin 6, master inflammatory regulator
- **TNFAIP3** - TNF-induced protein 3 (A20), NF-κB inhibitor
- **NFKB1** - Nuclear factor kappa B subunit 1
- **VEGFA** - Vascular endothelial growth factor A

### Novel Candidates
Genes not previously associated with mastitis but identified through our stability-based approach:
- `RELB`, `ARHGEF2`, `TNIP2`, `CACTIN`, `DHX9`, `IFNGR1`, `ATF4`, `RIPK2`, `IRAK2`

## 📈 Visualizations

<details>
<summary>Click to expand sample figures</summary>

### Module-Trait Relationships
![Module-trait heatmap](figures/module_trait_heatmap.png)

### Functional Enrichment
![Enrichment analysis](figures/functional_enrichment.png)

### Network Propagation
![Network propagation results](figures/network_propagation.png)

</details>

## 📚 Data Sources

| Dataset | Accession | Description |
|---------|-----------|-------------|
| RNA-Seq | [GSE51856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51856) | Bovine milk somatic cells, healthy vs. mastitis |
| Interactome | [PCNet](https://www.ndexbio.org/) | Parsimonious Composite Network |
| Reference | ARS-UCD1.2 | Bovine genome assembly |

## 🛠️ Dependencies

### R Packages
```r
WGCNA          # Weighted correlation network analysis
DESeq2         # Differential expression (preprocessing)
clusterProfiler # Functional enrichment
org.Bt.eg.db   # Bovine annotations
ggplot2        # Visualization
ComplexHeatmap # Heatmaps
```

### Python Packages
```python
networkx       # Graph analysis
ndex2          # NDEx network access
netcoloc       # Network colocalization
pandas         # Data manipulation
matplotlib     # Visualization
```

## 📝 Citation

If you use this code or methodology, please cite:

```bibtex
@article{shirazi2025rna,
  title={RNA stability: a novel perspective on gene regulatory networks in bovine mastitis},
  author={Shirazi, M.A. and Garaghani, [Supervisor]},
  journal={BMC Genomics},
  volume={26},
  pages={XX},
  year={2025},
  publisher={BioMed Central},
  doi={10.1186/s12864-025-12396-x}
}
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📬 Contact

- **Author**: [Your Name]
- **Email**: [your.email@example.com]
- **Institution**: Iran University of Medical Sciences
- **Supervisor**: Dr. Garaghani

---

<p align="center">
  <i>This work was conducted as part of an M.Sc. thesis in Animal Genetics at Iran University of Medical Sciences</i>
</p>
