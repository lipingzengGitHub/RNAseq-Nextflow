# RNA-seq Analysis Pipeline (Nextflow + Docker)

This repository provides a full RNA-seq pipeline using Nextflow with Docker support, including differential expression analysis and visualization.

# Project structure:
# .
# ├── .github/workflows/test_pipeline.yml
# ├── Dockerfile
# ├── nextflow.config
# ├── RNAseq-Nextflow.nf
# ├── scripts/R_script_DEG.R
# └── README.md



## 📦 Requirements
- Docker
- Nextflow

## 📁 Directory Structure
- `RawData/` — Paired-end FASTQ files
- `Ref/` — Reference genome and annotation (FASTA, GTF, HISAT2 index)
- `scripts/R_script_DEG.R` — R script for DE analysis

## 🚀 Running the Pipeline
```bash
nextflow run RNAseq-Nextflow.nf -profile docker
```

## 🧪 GitHub Actions CI
Automatically runs the pipeline on every push and pull request.

## 📊 Output
- `DEG_results.csv`
- `PCA_plot.png`
- `Volcano_plot.png`
- `Heatmap.png`
- `GO_plot.png`
