[![R Version](https://img.shields.io/badge/R-v4.5.2-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-v3.18-green.svg)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Transcriptomic Analysis of HIV Controllers vs cART-treated Patients

Transcriptomic analysis to identify differences in gene expression regulation between **HIV Controllers** and individuals with **chronic HIV infection receiving combination antiretroviral therapy (cART)** using public gene expression data.

---

## Project Overview

Human Immunodeficiency Virus (HIV) infection causes changes in host cellular gene regulation and immune responses. Although combination antiretroviral therapy (cART) effectively suppresses viral replication, it does not completely eliminate latent reservoirs or fully restore immune function.

A small group of individuals known as **HIV Controllers** can naturally maintain very low viral loads without antiretroviral treatment.

This project aims to explore transcriptomic differences between these groups to identify genes and biological pathways potentially involved in natural HIV suppression.

---

## Objective

- Identify **Differentially Expressed Genes (DEGs)**
- Compare transcriptomic profiles between:
  - HIV Controllers
  - HIV patients receiving cART
- Perform biological interpretation using:
  - Gene Ontology (GO)
  - KEGG Pathway enrichment

---

## Dataset

| Item | Description |
|---|---|
| Dataset | GSE108296 |
| Source | Gene Expression Omnibus (GEO) |
| Platform | Illumina HumanHT-12 V4.0 |
| Organism | *Homo sapiens* |
| Data Type | Microarray Expression Data |

---

# Workflow

<p align="center">
<img src="workflow/workflow.png" width="650">
</p>

---

## Methods

### 1. Data Acquisition
Gene expression data were downloaded from GEO using:

- GEOquery

### 2. Preprocessing
- Expression matrix extraction
- Log2 transformation
- Metadata processing
- Sample grouping:
  - Controller
  - cART

### 3. Differential Expression Analysis
Performed using:

- limma
- Linear modeling
- Empirical Bayes

Filtering criteria:

```text
Adjusted p-value < 0.05
|log2 Fold Change| ≥ 1
```

### 4. Visualization
Generated outputs:

- Boxplot
- Density Plot
- UMAP
- Volcano Plot
- Heatmap (Top 50 DEGs)

### 5. Functional Enrichment
Performed using:

- clusterProfiler
- Gene Ontology (GO)
- KEGG Pathway

---

## Tools & Packages

### Programming
- R
- RStudio

### Packages
- GEOquery
- limma
- ggplot2
- pheatmap
- dplyr
- clusterProfiler
- org.Hs.eg.db
- illuminaHumanv4.db
- umap

---

# Key Findings and Analysis

## Sample Distribution Validation

<p align="center">
<img src="results/boxplot.png" width="600">
</p>
Initial quality control via sample-wide boxplots confirmed highly uniform data distributions. The median expression lines across all profiles align perfectly at a log2 baseline scale of approximately 5.5. This structural consistency demonstrates that systemic batch effects are absent, confirming the data is ready for downstream linear modeling.

## Volcano Plot

Shows significantly upregulated and downregulated genes.

<p align="center">
<img src="results/volcano plot.png" width="800">
</p>
Using an expression cutoff of |log2 Fold Change| > 1 and an FDR-adjusted p-value < 0.05, significant transcripts were isolated and visualized using a volcano plot.
  - UP-regulated features (Red) represent biomarkers significantly elevated in target patient groups.
  - DOWN-regulated features (Blue) represent strongly suppressed genetic pathways.
  - Background stable features are colored gray, separating real biological signal from experimental noise.
  
---

## Heatmap (Top 50 Differentially Expressed Genes)

Visualization of expression patterns across samples.

<p align="center">
<img src="results/heatmap.png" width="800">
</p>
Hierarchical clustering of the Top 50 highly significant DEGs yielded perfect classification boundaries.
-The Controller cohort displays a uniform blocks of expression, showing strong down-regulated (blue) and up-regulated (red) clusters.
-The cART cohort exhibits a precise reciprocal pattern, demonstrating that this 50-gene panel serves as an exceptionally accurate diagnostic biomarker panel.

---

## Gene Ontology Enrichment

Top enriched biological processes.

<p align="center">
<img src="results/dotplot go.png" width="800">
</p>
Highly enriched terms visualized via a dotplot are heavily localized around functional immune networks, including mucosal immune responses, leukocyte migration, and acute inflammatory responses.

---

## KEGG Pathway Enrichment

Biological pathway enrichment analysis.

<p align="center">
<img src="results/dotplot kegg.png" width="800">
</p>
Exactly 5 highly specific pathways survived the strict false-discovery rate thresholds, forming a connected biological narrative:
-Phagocytosis & Neutrophil extracellular trap (NET) formation (Active innate immune cell defense mechanisms)
-Staphylococcus aureus infection & Leishmaniasis (Host-pathogen interaction interface)
-Hematopoietic cell lineage (Upstream immune cell differentiation)

---

## Output

Generated files:

- Differentially Expressed Genes (CSV)
- GO Enrichment Results
- KEGG Enrichment Results
- Figures and Visualizations

---

## Author

**Raysha Tryfhatya Nurhaidha**  
Biology Graduate  
Interest: Molecular Biology • Bioinformatics • Microbiology

Linkedln: `linkedin.com/in/rayshatn`
GitHub: `github.com/rayshatn`
Email: `rayshatryfhatya@gmail.com`

---
