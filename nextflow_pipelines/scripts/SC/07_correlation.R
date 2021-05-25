#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
library(gtools)
library(parallel)
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

file = args[1]
immune.combined <- readRDS(file)
identity <- args[2]

robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/06_correlation/02_spearman_all"
# Only calculate for the DE genes
de_all_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_genes.rds")
de_lnc_all<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_lnc.rds")
de_lnc <- unique(unlist(de_lnc_all))


calc_correlation_genes <- function(gene1_vector,gene2_vector,gene1, gene2, type = "spearman"){
  mask <- gene1_vector !=0 & gene2_vector !=0
  # calculate the correlation coefficient
  if(sum(mask)>50){
    #pc <- bayes.cor.test(c1, c2)
    pc <- cor.test(gene1_vector[mask], gene2_vector[mask], method = c(type))
    df <- data.frame( pval=pc$p.value, rho=pc$estimate)
    df$g1 <- gene1
    df$g2 <- gene2
    return(df)
  }else{
    return(NA)
  }
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



if(identity == "all"){
  subset <-  subset(immune.combined, features=unique(de_all_genes))
}else{
  subset <- subset(immune.combined, features=unique(de_all_genes), ident=identity)
}

expressionmatrix <- subset@assays$RNA@counts
colnames(expressionmatrix) <- Idents(subset)

# Create all the pairs
pairs <- (expand.grid.unique(rownames(expressionmatrix)[rownames(expressionmatrix) %in% de_lnc], rownames(expressionmatrix)))
pairs <- as.list(data.frame(t(pairs), stringsAsFactors = F))

print(identity)
length(unique(pairs))
print("Beginning Calculation .. ")


unique(unlist(lapply(pairs, function(x) length(x))))



correlations <- mclapply(pairs, function(genes) {calc_correlation_genes(expressionmatrix[genes[1], ], expressionmatrix[genes[2],], genes[1], genes[2])}, mc.cores = 48)
print("DONE!")
saveRDS(correlations, file.path(robjectsdir, paste("pearson_correlations_",identity,".rds")))
