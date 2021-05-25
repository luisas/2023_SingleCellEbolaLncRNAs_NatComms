#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
library(gtools)
library(parallel)

options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

file = args[1]
infected_myeloids_24 <- readRDS(file)


robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/03_viralload"
robjectsdir <- args[2]
lnc <- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/lnc.rds")
pc <- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/pc.rds")

# Only calculate for the DE genes 
expression_matrix <- infected_myeloids_24@assays$RNA@data
viral_load <- infected_myeloids_24$viral_load


cal_corr_viral_load_average_window <- function(expression_matrix, viral_load,gene,  window = 20){
  df <- data.frame(exp = expression_matrix[gene, ], viral_load)
  df <- df[df$exp !=0, ]
  # Order cells by viral load
  df_ordered <- df[order(df$viral_load),]
  if(nrow(df_ordered) > 50 ){
    # Average by 100  window gene expression
    n <- window
    df_ordered$bin <- cut(1:length(df_ordered$viral_load), breaks =c( seq(0, length(df_ordered$viral_load), n), ceiling(length(df_ordered$viral_load)/n)*n) )
    print(c( seq(0, length(df_ordered$viral_load), n), ceiling(length(df_ordered$viral_load)/n)*n))
    avg_bins <- df_ordered %>% dplyr::group_by(bin) %>% dplyr::summarise(avg = mean(exp))
    avg_bins_vl <- df_ordered %>% dplyr::group_by(bin) %>% dplyr::summarise(avg = mean(viral_load))
    df_ordered$avg <- NULL
    df_ordered$avg <- unlist(lapply(df_ordered$bin, function(x) avg_bins[x== avg_bins$bin,]$avg))
    df_ordered$avg_vl <- unlist(lapply(df_ordered$bin, function(x) avg_bins_vl[x== avg_bins_vl$bin,]$avg))
    cor <- cor.test(df_ordered$avg,df_ordered$viral_load, method = c("spearman"))
    return(data.frame(gene = gene, rho = cor$estimate, pval = cor$p.value,  n_cell = nrow(df_ordered)))
  }
}


correlations <- mclapply(unique(lnc), function(gene) cal_corr_viral_load_average_window(expression_matrix, viral_load, gene, 10),  mc.cores = 48)
correlations_pc <- mclapply(unique(pc), function(gene) cal_corr_viral_load_average_window(expression_matrix, viral_load, gene, 100),  mc.cores = 48)


correlations_df <- do.call(rbind, correlations)
correlations_pc_df <- do.call(rbind, correlations_pc)

saveRDS(correlations_df, file.path(robjectsdir, paste("spearman-correlations_lnc.rds")))
saveRDS(correlations_pc_df, file.path(robjectsdir, paste("spearman-correlations_pc.rds")))
