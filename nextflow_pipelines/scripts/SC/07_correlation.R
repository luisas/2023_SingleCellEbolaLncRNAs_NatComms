#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
library(gtools)
library(parallel)
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

file = args[1]
robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/06_correlation/01_pearson"
immune.combined <- readRDS(file)

# Only calculate for the DE genes 
de_all_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_genes.rds")
immune.combined <- subset(immune.combined, features=de_all_genes)


calc_correlation_vectors <- function(gene1_vector,gene2_vector, type = "pearson"){
  mask <- gene1_vector !=0 & gene2_vector !=0
  # calculate the correlation coefficient
  if(sum(mask)>50){ 
    #pc <- bayes.cor.test(c1, c2) 
    pc <- cor.test(gene1_vector[mask], gene2_vector[mask], method = c(type))
    df <- data.frame( pval=pc$p.value, rho=pc$estimate)
    return(df)
  }else{ 
    return(NA)
  }
}

get_correlation_df_genes <- function(gene1_vector, gene2_vector, gene1, gene2, identities, identity){
   mask <- identities == identity 
   df <- calc_correlation_vectors(gene1_vector[mask], gene2_vector[mask])
   if(!is.na(df)){
     df$celltype <- identity
     df$g1 <- gene1
     df$g2 <- gene2
     return(df)
   }
}


get_correlations <- function(gene1_vector, gene2_vector, gene1, gene2, identities){
  lapply(identities)
}

expand.grid.unique <- function(x, y, include.equals=FALSE){
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}


cl <- makeCluster(4, type="FORK")
# Prepare count matrix 
expressionmatrix <- immune.combined@assays$RNA@counts
colnames(expressionmatrix) <- Idents(immune.combined)
# Create all the pairs
pairs <- (expand.grid.unique(rownames(expressionmatrix), rownames(expressionmatrix)))

clusterExport(cl,list("expressionmatrix","pairs"),envir=globalenv())

correlations <-parApply(cl, pairs,1, function(genes) do.call(rbind,lapply(as.character(unique(colnames(expressionmatrix))), function(identity) get_correlation_df_genes(expressionmatrix[genes[1],], expressionmatrix[genes[2],], genes[1], genes[2], colnames(expressionmatrix), identity)
)))
final <- Reduce("rbind", correlations)
saveRDS(final, file.path(dirname(robjectsdir), "pearson_correlations.rds"))

stopCluster(cl)









