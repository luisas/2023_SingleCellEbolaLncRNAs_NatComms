#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
library(gtools)
library(parallel)
library(reshape2)
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

file = args[1]
robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/06_correlation/01_pearson"
immune.combined <- readRDS(file)

# Only calculate for the DE genes 
de_all_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_genes.rds")
de_lnc_all<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_lnc.rds")
de_lnc <- unique(unlist(de_lnc_all))




BaCo <- function(X){
  alpha0 <- rep(1/nrow(X),ncol(X))
  beta0=1-alpha0
  nrowsX <- nrow(X)
  k <- ncol(X)
  cs <- colSums(X)
  alphas <- alpha0 + X
  betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
  alphasPLUSbetas <- alphas + betas
  Psi <- alphas/alphasPLUSbetas - matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE) 
  var_vec <- as.matrix( ( rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + rowSums(Psi^2) )/k )
  cov_mtrx <- (Psi %*% t(Psi))/k
  Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
  diag(Bcorrvals) <- 1
  Bcorrvals
}




# Compute the Bayesian correlation matrix
get_baco <- function(expressionmatrix, ident, genes){
  X <- expressionmatrix[rownames(expressionmatrix) %in% genes,colnames(expressionmatrix) == ident ]
  B <- BaCo(X)
  rownames(B) <- colnames(B) <- genes
  df <- melt(as.matrix(B))
  df$ident <- ident
  return(df)
}


#cl <- makeCluster(4, type="FORK")
# Prepare count matrix 
subset <- subset(immune.combined, features=unique(de_all_genes))
expressionmatrix <- as.matrix(subset@assays$RNA@counts)
colnames(expressionmatrix) <- Idents(subset)
identities <- unique(unlist(colnames(expressionmatrix)))
# Create all the pairs

correlations <- lapply(identities, function(ident) get_baco(expressionmatrix, ident, de_all_genes))

#clusterExport(cl,list("expressionmatrix", "de_all_genes", "identities"),envir=globalenv())
#correlations <-parLapply(cl, identities, function(ident) get_baco(expressionmatrix, ident, de_all_genes)) 
final <- Reduce("rbind", correlations)
saveRDS(final, file.path(dirname(robjectsdir), "bayesan.rds"))
#stopCluster(cl)









