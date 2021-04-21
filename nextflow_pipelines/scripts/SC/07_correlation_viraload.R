#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
library(gtools)
library(parallel)

options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

# Read in Seurat Object 
file = args[1]
myeloids <- readRDS(file)
expression_matrix <- myeloids@assays$RNA@data

# Output directory 
#robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/03_viralload"
robjectsdir <- args[2]
dir.create(robjectsdir, recursive = T,  showWarnings = F)
genes <- unique(rownames(myeloids))


calc_correlation_viralload <- function(gene1_vector,viral_load,gene1, type = "spearman"){
  mask <- gene1_vector!=0 
  #mask <- rep(TRUE, length(gene1_vector))
  if(sum(mask)>30){ 
  pc <- cor.test(gene1_vector[mask], viral_load[mask], method = c(type))
  df <- data.frame( pval=pc$p.value, rho=pc$estimate)
  df$gene <- gene1
  return(df)
  }
}

viral_load <- myeloids$viraload 

correlations <- mclapply(unique(genes), function(gene) calc_correlation_viralload(expression_matrix[rownames(expression_matrix) == gene, ], viral_load, gene),  mc.cores = 48)

correlations_df <- do.call(rbind, correlations)

saveRDS(correlations_df, file.path(robjectsdir, paste("ViralLoad_spearman-correlations_standard.rds")))


