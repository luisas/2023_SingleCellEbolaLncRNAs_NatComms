#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
library(gtools)
library(parallel)
options(future.globals.maxSize = 10000 * 1024^2)
args = commandArgs(trailingOnly=TRUE)

file <-  args[1]
out <- args[2]

immune.combined <- readRDS(file)


immune.combined <- NormalizeData(immune.combined)
immune.combined <- FindVariableFeatures(immune.combined)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
immune.combined <- ScaleData(immune.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(immune.combined))
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca",dims = 1:20 )
immune.combined <- FindClusters(immune.combined, resolution = 0.2)

saveRDS(immune.combined, out) 