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

robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/06_correlation/01_spearman"
# Only calculate for the DE genes 
robjectsdir_stats <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/"
all_lncrnas <- readRDS(file.path(robjectsdir_stats, "all_lncrnas.rds"))
annotated_mrnas <- readRDS(file.path(robjectsdir_stats,"annotated_mrnas.rds"))

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


mamu_dr <- c("MAMU-DRA", "MAMU-DRB1")
genes <- c(all_lncrnas, mamu_dr)
if(identity == "all"){
  subset <-  subset(immune.combined, features=unique(genes))
}else{
  subset <- subset(immune.combined, features=unique(genes), ident=identity)
}

expressionmatrix <- subset@assays$RNA@counts
colnames(expressionmatrix) <- Idents(subset)

# Create all the pairs

pairs <- (expand.grid.unique(rownames(expressionmatrix)[rownames(expressionmatrix) %in% all_lncrnas], rownames(expressionmatrix)[gsub("-unknown","",rownames(expressionmatrix)) %in% unique(mamu_dr)]))
pairs <- as.list(data.frame(t(pairs), stringsAsFactors = F))


print(identity)
length(unique(pairs))
print("Beginning Calculation .. ")


correlations <- mclapply(pairs, function(genes) {calc_correlation_genes(expressionmatrix[genes[1], ], expressionmatrix[genes[2],], genes[1], genes[2])}, mc.cores = 48)
saveRDS(correlations, file.path(robjectsdir, paste("Spearman_correlations_MAMU_alllnc_",identity,".rds")))
print("DONE lnc!")


genes <- c(annotated_mrnas, mamu_dr)
if(identity == "all"){
  subset <-  subset(immune.combined, features=unique(genes))
}else{
  subset <- subset(immune.combined, features=unique(genes), ident=identity)
}
expressionmatrix <- subset@assays$RNA@counts
colnames(expressionmatrix) <- Idents(subset)

pairs <- (expand.grid.unique(rownames(expressionmatrix)[rownames(expressionmatrix) %in% annotated_mrnas], rownames(expressionmatrix)[gsub("-unknown","",rownames(expressionmatrix)) %in% unique(mamu_dr)]))
pairs <- as.list(data.frame(t(pairs), stringsAsFactors = F))
correlations <- mclapply(pairs, function(genes) {calc_correlation_genes(expressionmatrix[genes[1], ], expressionmatrix[genes[2],], genes[1], genes[2])}, mc.cores = 48)


saveRDS(correlations, file.path(robjectsdir, paste("Spearman_correlations_MAMU_allpc_",identity,".rds")))










