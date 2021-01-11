#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
library(gtools)
library(parallel)

options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

file = args[1]
myeloids <- readRDS(file)
expression_matrix <- myeloids@assays$RNA@data


robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/03_viralload"
robjectsdir <- args[2]
lnc <- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/lnc.rds")
pc <- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/pc.rds")

# Only calculate for the DE genes 



calc_correlation_viralload <- function(gene1_vector,viral_load,gene1, type = "spearman"){
  mask <- gene1_vector!=0 
  #mask <- rep(TRUE, length(gene1_vector))
  if(sum(mask)>30){ 
  pc <- cor.test(gene1_vector[mask], viral_load[mask], method = c(type))
  df <- data.frame( pval=pc$p.value, rho=pc$estimate)
  df$gene <- gene1
  return(df)
  }
  else{
    return(NA)
  }
  
}

viral_load <- myeloids$viral_load 


correlations <- mclapply(unique(lnc), function(gene) calc_correlation_viralload(expression_matrix[rownames(expression_matrix) == gene, ], viral_load, gene),  mc.cores = 48)
correlations_pc <- mclapply(unique(pc), function(gene) calc_correlation_viralload(expression_matrix[rownames(expression_matrix) == gene, ], viral_load, gene),  mc.cores = 48)


correlations_df <- do.call(rbind, correlations)
correlations_pc_df <- do.call(rbind, correlations_pc)

saveRDS(correlations_df, file.path(robjectsdir, paste("spearman-correlations_lnc.rds")))
saveRDS(correlations_pc_df, file.path(robjectsdir, paste("spearman-correlations_pc.rds")))










