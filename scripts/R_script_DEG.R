#This R script performs:

#DESeq2-based differential expression
#PCA plot
#Volcano plot for DEG
#Heatmap for DEG
#GO term enrichment with clusterProfiler


#!/usr/bin/env Rscript

# Required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.At.tair.db)  # Replace with appropriate OrgDb if not Arabidopsis

# Input: path to folder containing sample counts
args <- commandArgs(trailingOnly = TRUE)
count_dir <- args[1]

# Gather count files
count_files <- list.files(count_dir, pattern = "*.txt", full.names = TRUE)
sample_ids <- gsub(".*/(.*)/counts.txt", "\\1", count_files)

# Build count matrix
counts_list <- lapply(count_files, function(f) read.table(f, header = FALSE, row.names = 1))
count_matrix <- do.call(cbind, lapply(counts_list, `[`, , 1))
colnames(count_matrix) <- sample_ids

# Experimental design — adjust for your case!
condition <- factor(c(rep("A", 3), rep("B", 3)))  # <-- Update this!
coldata <- data.frame(row.names = colnames(count_matrix), condition = condition)

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered), file = "DEG_results.csv")

# Apply cutoff: padj < 0.05 and |log2FC| >= 1
res_df <- as.data.frame(resOrdered)
res_df$gene <- rownames(res_df)
res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, "DEG", "Not Significant")

# PCA Plot
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  ggtitle("PCA Plot") +
  theme_minimal()
ggsave("PCA_plot.png")

# Volcano Plot
ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.7) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  ggtitle("Volcano Plot") +
  xlab("log2 Fold Change") +
  ylab("-log10(FDR)")
ggsave("Volcano_plot.png")

# Heatmap of DEGs only
deg_genes <- rownames(subset(res_df, padj < 0.05 & abs(log2FoldChange) >= 1))
heatmap_mat <- assay(vsd)[deg_genes, ]
pheatmap(heatmap_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         scale = "row",
         filename = "Heatmap.png")

# GO Enrichment
ego <- enrichGO(gene = deg_genes,
                OrgDb = org.At.tair.db,  # Change for your organism
                keyType = "TAIR",        # e.g., "ENSEMBL", "SYMBOL"
                ont = "BP",
                readable = TRUE)
dotplot(ego) + ggtitle("GO Enrichment")
ggsave("GO_plot.png")


