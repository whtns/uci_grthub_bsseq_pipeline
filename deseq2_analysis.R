#!/usr/bin/env Rscript 
# DESeq2 analysis script
# Usage: Rscript deseq2_analysis.R counts.txt metadata.csv output_dir

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 3) stop("Usage: Rscript deseq2_analysis.R counts.txt metadata.csv output_dir")

counts_file <- args[1]
meta_file <- args[2]
out_dir <- args[3]

dir.create(out_dir, showWarnings=FALSE)

counts <- read.table(counts_file, header=TRUE, row.names=1)

# Extract the count matrix (assuming first 5 columns are metadata)
count_matrix <- counts[,6:ncol(counts)]

# format column name to match sample name
format_sample_id <- function(colname){
    colname |>
    str_remove("results.hisat2_alignment\\.") %>%
  # Remove the suffix
  str_remove("_align_sorted.bam") %>%
  # Replace all periods with hyphens
  str_replace_all("\\.", "-")
}

colnames(count_matrix) <- sapply(colnames(count_matrix), format_sample_id)

meta <- read.csv(meta_file, row.names=1)

print(rownames(meta))
print(colnames(count_matrix))
# Ensure sample names match
count_matrix <- count_matrix[, rownames(meta)]

dds <- DESeqDataSetFromMatrix(countData=count_matrix, 
  colData=meta, 
  design=~condition)

dds <- DESeq(dds)

saveRDS(dds, file=file.path(out_dir, "dds.rds"))

res <- results(dds)
write.csv(as.data.frame(res), file=file.path(out_dir, "deseq2_results.csv"))
