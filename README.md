# RNA-seq Analysis Pipeline (Nextflow + Docker)

This repository provides a full RNA-seq pipeline using Nextflow with Docker support, including differential expression analysis and visualization.

# Project Structure
.
├── .github/workflows/test_pipeline.yml   # GitHub Actions CI workflow to test the pipeline
├── Dockerfile                            # Docker environment definition for reproducible execution
├── nextflow.config                       # Nextflow configuration file (resources, params, paths)
├── RNAseq-Nextflow.nf                    # Main Nextflow pipeline script for RNA-seq processing
├── scripts/
│   └── R_script_DEG.R                    # R script for DEG analysis using DESeq2 (used downstream)
└── README.md                             # Project overview, usage instructions, dependencies


## 📦 Requirements
- Docker
- Nextflow

## 📁 Directory Structure
- `RawData/` — Paired-end FASTQ files
- `Ref/` — Reference genome and annotation (FASTA, GTF, HISAT2 index)
- `scripts/R_script_DEG.R` — R script for DEG analysis

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
