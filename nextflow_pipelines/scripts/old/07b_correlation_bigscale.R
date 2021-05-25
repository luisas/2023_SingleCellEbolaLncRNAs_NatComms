#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(bigSCale))
shhh(library(Seurat))

args = commandArgs(trailingOnly=TRUE)
file = args[1]
immune.combined <- readRDS(file)
de_all_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_genes.rds")
immune.combined <- subset(immune.combined, ident = "Monocyte")
expr.ctl <- immune.combined@assays$RNA@counts
gene.names <- rownames(immune.combined)[rownames(immune.combined) %in% de_all_genes]

results.ctl=compute.network(expr.data = expr.ctl,gene.names = gene.names[1:100], clustering = "direct")

saveRDS(results.ctl,"/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/06_correlation/02_bigscale/in_vivo_correlation_test.rds" )



