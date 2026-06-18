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
<img src="workflow/workflow.png" width="850">
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

# Analysis Outputs

## Volcano Plot

Shows significantly upregulated and downregulated genes.

<p align="center">
<img src="results/volcano plot.png" width="800">
</p>

---

## Heatmap (Top 50 Differentially Expressed Genes)

Visualization of expression patterns across samples.

<p align="center">
<img src="results/heatmap.png" width="800">
</p>

---

## Gene Ontology Enrichment

Top enriched biological processes.

<p align="center">
<img src="results/dotplot go.png" width="800">
</p>

---

## KEGG Pathway Enrichment

Biological pathway enrichment analysis.

<p align="center">
<img src="results/dotplot kegg.png" width="800">
</p>

---

## Key Findings

Significant transcriptomic differences were observed between HIV Controllers and cART-treated patients.

Enrichment analysis highlighted pathways associated with:

- Immune response
- Chemotaxis
- Leukocyte migration
- Antibody-dependent cellular cytotoxicity (ADCC)
- Phagosome
- Neutrophil Extracellular Trap (NET) formation

These findings suggest that natural HIV control may involve more coordinated innate and adaptive immune regulation.

---

## Output

Generated files:

- Differentially Expressed Genes (CSV)
- GO Enrichment Results
- KEGG Enrichment Results
- Figures and Visualizations

---

## Author

**Raysha**  
Biology Graduate  
Interest: Molecular Biology • Bioinformatics • Microbiology

Linkedln: `linkedin.com/in/rayshatn`
GitHub: `github.com/rayshatn`
Email: `rayshatryfhatya@gmail.com`

---
