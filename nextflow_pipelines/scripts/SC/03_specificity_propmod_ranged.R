
library(Seurat)
library(gtools)
library(parallel)
shhh <- suppressPackageStartupMessages
#shhh(library(SingleCellExperiment))
shhh(library(Seurat))
shhh(library(stringr))
shhh(library(dplyr))
shhh(library(MAST))

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
#      This script calculates the cell-type specificity score of a gene
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


# Read in argumants
args = commandArgs(trailingOnly=TRUE)
# Seurat object (pre-processed and normalized)
df_celltype <- readRDS(args[1])
# Where the result file should be saved
outfile <- args[2]
df_celltype_permutation <-readRDS(args[3])

# Prepare dir for storing output file
dir.create(dirname(file.path(outfile)), showWarnings = FALSE)
df_celltype$observed <- df_celltype$n_cells

get_score_prop_mod <- function(df_celltype, gene){
  df_celltype_gene <- df_celltype[df_celltype$gene_id == gene, ]
  df_celltype_gene$prop <- df_celltype_gene$tot_cells/sum(df_celltype_gene$tot_cells)
  # Calculate the expected 
  df_celltype_gene$expected <- sum(df_celltype_gene$n_cells)*df_celltype_gene$prop
  df_celltype_gene$obs_prop <- df_celltype_gene$observed/sum(df_celltype_gene$observed)
  
  calc_score_1  <- function(df_celltype_gene, ident){
    df_celltype_gene_celltype <- df_celltype_gene[as.character(df_celltype_gene$ident) == ident,];  
    #s <- (df_celltype_gene_celltype$obs_prop -df_celltype_gene_celltype$prop)
    s <- max((df_celltype_gene_celltype$obs_prop -df_celltype_gene_celltype$prop)*100/(100-100*df_celltype_gene_celltype$prop))[1]
    return(data.frame(score = s, identity = ident, stringsAsFactors = F))
  }
  score <- Reduce("rbind", (lapply(unique(df_celltype$ident), function(ident){calc_score_1(df_celltype_gene, ident)} ))) 
  score$gene <- gene
  return(score[score$score == max(score$score),][1,])
}

identities <- as.character(unique(df_celltype$ident))

# Calculate the real ones 
#specificity_scores <- Reduce("rbind", lapply(unique(df_celltype$gene_id), function(gene) get_score_prop_mod(df_celltype, gene)))
#saveRDS(specificity_scores, outfile)


# Calculate it for the permutation 
df_celltype_permutation$observed <- df_celltype_permutation$n_cells

specificity_scores <- Reduce("rbind", lapply(unique(df_celltype_permutation$gene_id), function(gene) get_score_prop_mod(df_celltype_permutation, gene)))
saveRDS(specificity_scores, file.path(dirname(outfile), paste("Permutations_propmodranged_1.rds", sep = "")))


