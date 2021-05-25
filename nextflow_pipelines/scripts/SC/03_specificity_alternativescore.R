
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

# Prepare dir for storing output file
dir.create(dirname(file.path(outfile)), showWarnings = FALSE)
df_celltype$observed <- df_celltype$n_cells


# df_celltype contains stats about the number of cells and summaries of expression of each gene

get_score_prop_mod <- function(df_celltype, gene){
  
  # Extract gene entry from summary DF 
  df_celltype_gene <- df_celltype[df_celltype$gene_id == gene, ]
  
  # What are the real proportion of cell-types in dataset 
  df_celltype_gene$prop <- df_celltype_gene$tot_cells/sum(df_celltype_gene$tot_cells)
  
  # Calculate the expected expression proportion per cell-type
  df_celltype_gene$expected <- sum(df_celltype_gene$n_cells)*df_celltype_gene$prop
  # Calculate the observed expression proportion per cell-type
  df_celltype_gene$obs_prop <- df_celltype_gene$observed/sum(df_celltype_gene$observed)
  
  calc_score_1  <- function(df_celltype_gene, ident){
    # Extract one cell-type
    df_celltype_gene_celltype <- df_celltype_gene[as.character(df_celltype_gene$ident) == ident,];  
    # Get the maximum 
    # Out of the difference of observed and expected and the maximum prop that can be obtained by that cell-type
    s <- max((df_celltype_gene_celltype$obs_prop -df_celltype_gene_celltype$prop)*100/(100-100*df_celltype_gene_celltype$prop))[1]
    return(data.frame(score = s, identity = ident, stringsAsFactors = F))
  }
  score <- Reduce("rbind", (lapply(unique(df_celltype$ident), function(ident){calc_score_1(df_celltype_gene, ident)} ))) 
  score$gene <- gene
  return(score[score$score == max(score$score),][1,])
}

# Cell-types
identities <- as.character(unique(df_celltype$ident))

# Calculate for the complete dataset 
specificity_scores <- Reduce("rbind", lapply(unique(df_celltype$gene_id), function(gene) get_score_prop_mod(df_celltype, gene)))
saveRDS(specificity_scores, outfile)




