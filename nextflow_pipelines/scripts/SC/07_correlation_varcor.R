#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
shhh(library(scran))
shhh(library(SingleCellExperiment))
library(gtools)
library(parallel)
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

file = args[1]
robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/06_correlation/01_pearson"
immune.combined <- readRDS(file)

# Only calculate for the DE genes 
de_all_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_genes.rds")
subset_ic <- subset(immune.combined, features=unique(de_all_genes)[1:100])
de_lnc_all<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_lnc.rds")
de_lnc <- unique(unlist(de_lnc_all))
identities <- unique(Idents(subset_ic))





#cl <- makeCluster(4, type="FORK")
#clusterExport(cl,list("subset", "identities"),envir=globalenv())
print(identities)
identity <- "Monocyte"
#correlations <-parLapply(cl, identities,  function(ident) get_cor(ident))
#correlations <-lapply( identities, function(ident) )
imm.comb <- subset(subset_ic, idents = identity)
imm.comb_sc <- as.SingleCellExperiment(imm.comb)
var.cor <- scran::correlatePairs(imm.comb_sc);
var.cor$celltype <- ident;

saveRDS(var.cor, file.path(dirname(robjectsdir), "pearson_correlations_varcor.rds"))

#stopCluster(cl)









