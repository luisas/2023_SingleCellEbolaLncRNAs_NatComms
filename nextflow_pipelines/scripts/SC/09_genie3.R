#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(GENIE3))
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

immune.combined<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds")
de_all_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_genes.rds")
all_lncrnas <- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/all_lncrnas.rds")
annotated_mrnas <- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/annotated_mrnas.rds")
de_lnc <- unique(de_all_genes[de_all_genes %in% all_lncrnas])
de_pc <- unique(de_all_genes[de_all_genes %in% annotated_mrnas])

table(de_pc%in% rownames(immune.combined))
table(de_lnc%in% rownames(immune.combined))


immune.combined.subset <- subset(immune.combined, features=unique(de_all_genes), idents = "Monocyte")

exp_mat <- immune.combined.subset@assays$integrated@data

colnames(exp_mat) <- Idents(immune.combined.subset)
exprMatr <- as.matrix(exp_mat)

colnames(exprMatr)
rownames(exprMatr)

set.seed(123) # For reproducibility of results
#weightMat <- GENIE3(exprMatr)
weightMat <- GENIE3(exprMatr, regulators=de_lnc,  nCores=16, verbose=TRUE)

saveRDS(weightMat, "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/weigth_matrix_degenesONLY_regulatorsLNC_allcells.rds")
