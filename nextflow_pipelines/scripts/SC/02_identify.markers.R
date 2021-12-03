#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))

immune.combined <- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10_NOVELFILTERTEST_counts/05_RObjects/03_prep/03_immune.combined_preprocessed_noHBB.rds")
pbmc.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(pbmc.markers,"/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10_NOVELFILTERTEST_counts/05_RObjects/03_prep/markers_exvivo_normal_noHBB.rds")