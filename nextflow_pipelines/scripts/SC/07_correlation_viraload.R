#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
library(gtools)
library(parallel)

options(future.globals.maxSize = 10000 * 1024^2)
args = commandArgs(trailingOnly=TRUE)


# 0. Read command line arguments
file <-  args[1]
type <- args[2]
out <- args[3]
out2 <- args[4]
immune.combined <- readRDS(file)


#immune.combined <- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready_infectionstatus.rds"))


# 1. Extract Monocytes 
cellsubset <- subset(immune.combined, ident = "Monocyte")
# 2. Extract Infected cells 
cellsubset <- cellsubset[,cellsubset$viral_load>0]
# 3. Retain  late only 
if(type == "invivo"){
  cellsubset <- cellsubset[,cellsubset$group_dpi %in% c("DPI005","DPI006","DPI007","DPI008")]
}else if(type == "exvivo"){
  cellsubset <- cellsubset[,cellsubset$cond == "live"]
  cellsubset <- cellsubset[,cellsubset$dpi == "H024"]
}
saveRDS(cellsubset, out2)

# Extract expression matrix
expression_matrix <- cellsubset@assays$RNA@data


# METHOD 1 --- correlation on all values 
calc_correlation_viralload <- function(expression_matrix,viral_load,gene1, type = "spearman", dropzeros = F){
  gene1_vector <- expression_matrix[gene1, ]
  if(dropzeros == F){
    mask <- rep(TRUE, length(gene1_vector))
  }else{
    mask <- gene1_vector!=0 
  }
  if(sum(mask)>30){ 
    pc <- cor.test(gene1_vector[mask], viral_load[mask], method = c(type))
    df <- data.frame( pval=pc$p.value, rho=pc$estimate)
    df$gene <- gene1
    return(df)
  }
}

# METHOD 2 ----- correlation on averaged windows ´
calc_correlation_viralload_average <- function(expression_matrix, viral_load,gene,  window = 10, dropzeros = F){
  # 1. Create a DF with gene expression and viral load 
  df <- data.frame(exp = expression_matrix[gene, ], viral_load)
  # 2. Remove zeros
  if(dropzeros == T){
    df <- df[df$exp != 0,]
  }
  # 2. Order cells by viral load
  df_ordered <- df[order(df$viral_load),]
  if(nrow(df_ordered) > 30 ){
    # Average gene expression by n-cells sliding window
    n <- window
    df_ordered$bin <- cut(1:length(df_ordered$viral_load), breaks =c( seq(0, length(df_ordered$viral_load), n), ceiling(length(df_ordered$viral_load)/n)*n) )
    #print(c( seq(0, length(df_ordered$viral_load), n), ceiling(length(df_ordered$viral_load)/n)*n))
    avg_bins <- df_ordered %>% dplyr::group_by(bin) %>% dplyr::summarise(avg = mean(exp))
    avg_bins_vl <- df_ordered %>% dplyr::group_by(bin) %>% dplyr::summarise(avg = mean(viral_load))
    df_ordered$avg <- NULL
    df_ordered$avg <- unlist(lapply(df_ordered$bin, function(x) avg_bins[x== avg_bins$bin,]$avg))
    df_ordered$avg_vl <- unlist(lapply(df_ordered$bin, function(x) avg_bins_vl[x== avg_bins_vl$bin,]$avg))
    print(df_ordered)
    # Calculate correlation value
    cor <- cor.test(df_ordered$avg,df_ordered$viral_load, method = c("spearman"))
    return(data.frame(gene = gene, rho = cor$estimate, pval = cor$p.value,  n_cell = nrow(df_ordered)))
  }
}
# log10 viral load 
viral_load <- log10(cellsubset$viral_load)
genes <- unique(rownames(cellsubset))

# Get directory 
outputdir <- dirname(out)

# Spearman on  all cells (A_approach) 
correlations <- mclapply(unique(genes), function(gene) calc_correlation_viralload(expression_matrix, viral_load, gene, dropzeros = F),  mc.cores = 48)
correlations_df <- do.call(rbind, correlations)
saveRDS(correlations_df, file.path(outputdir, "Approach_A.rds"))

# Spearman only on cells in which the gene is expressed (B_approach) 
correlations <- mclapply(unique(genes), function(gene) calc_correlation_viralload(expression_matrix, viral_load, gene, dropzeros = T),  mc.cores = 48)
correlations_df <- do.call(rbind, correlations)
saveRDS(correlations_df, file.path(outputdir, "Approach_B.rds"))


# Spearman on 10-cells sliding window on all cells  (C_approach)
correlations <- mclapply(unique(genes), function(gene) calc_correlation_viralload_average(expression_matrix, viral_load, gene, dropzeros = F),  mc.cores = 48)
correlations_df <- do.call(rbind, correlations)
saveRDS(correlations_df, file.path(outputdir, "Approach_C.rds"))

# Spearman on 10-cells sliding window on cells in which the gene is expressed  (D_approach)
correlations <- mclapply(unique(genes), function(gene) calc_correlation_viralload_average(expression_matrix, viral_load, gene, dropzeros = T),  mc.cores = 48)
correlations_df <- do.call(rbind, correlations)
saveRDS(correlations_df, file.path(outputdir, "Approach_D.rds"))
















