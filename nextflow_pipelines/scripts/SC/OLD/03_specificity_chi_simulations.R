
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

dir.create(dirname(file.path(outfile)), showWarnings = FALSE)

calc_score_gene <- function(df_celltype, gene){
  df_gene <- df_celltype[df_celltype$gene_id == gene,]
  if(any(df_gene$n_cells) == 0 ){
    df <- data.frame(pval= NA, chi = NA, gene = gene)
    rownames(df) <- df$gene
    return(df)
  }
  test <- prop.test( df_gene$n_cells, n=rep(sum(df_gene$n_cells), 3), p=(df_gene$tot_cells/sum(df_gene$tot_cells)))
  df <- data.frame(pval= test$p.value, chi = test$statistic, gene = gene)
  print(gene)
  rownames(df) <- df$gene
  return(df)
}

cl <- makeCluster(4, type="FORK")
clusterExport(cl,list("df_celltype","calc_score_gene"),envir=globalenv())
spec_scores <- parLapply(cl, as.character(unique((df_celltype$gene_id))), function(gene) calc_score_gene(df_celltype, gene))
final <- Reduce("rbind", spec_scores)
saveRDS(final, outfile)
stopCluster(cl)


