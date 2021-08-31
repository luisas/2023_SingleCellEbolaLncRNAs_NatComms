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



# 1. Extract Monocytes 
cellsubset <- subset(immune.combined, ident = "Monocyte")
# 2. Extract Infected cells 
cellsubset <- cellsubset[,cellsubset$infection == "Infected"]
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
# log10 viral load 
viral_load <- log10(cellsubset$viral_load)
# Genes 
genes <- unique(rownames(cellsubset))
# Get directory 
outputdir <- dirname(out)


# METHOD 1 --- correlation on all values 
calc_correlation_viralload <- function(expression_matrix,viral_load,gene1, type = "spearman", dropzeros = T){
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


# Spearman only on cells in which the gene is expressed (B_approach) 
correlations <- mclapply(unique(genes), function(gene) calc_correlation_viralload(expression_matrix, viral_load, gene, dropzeros = T),  mc.cores = 48)
correlations_df <- do.call(rbind, correlations)
correlations_df$qval <- p.adjust(correlations_df$pval, method = "BH", n = length(correlations_df$pval))
saveRDS(correlations_df, file.path(out))



















