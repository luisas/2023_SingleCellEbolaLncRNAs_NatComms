#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(SingleCellExperiment))
shhh(library(Seurat))
shhh(library(stringr))
shhh(library(dplyr))

options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

file = file.path("/home/luisas/Desktop/cluster/proj/code/ebola/src/scripts/Analysis/results/seurat_pbmc_rhemac10_merged_aftercellandgeneqc_afterScrublet.rds")
file = file.path(args[1])
integration_variable = file.path(args[2]) # e.g. "cond"

# Extract dir and file names from file path
results_folder = dirname(file)
file_name_no_extension = gsub(pattern = "\\.rds$", "", basename(file))

# Read in Seurat Object: should already have all metadata in order
pbmc <- readRDS(file)


pbmc <- SplitObject(pbmc, split.by = integration_variable)

pbmc <- lapply(X = pbmc, FUN = function(x) {
  x <- NormalizeData(x)
})

print("LogNorm completed!")
new_name = paste(file_name_no_extension, "_00_normalized.rds", sep  ="_")
saveRDS(pbmc, file.path(results_folder,new_name))


immune.anchors <- FindIntegrationAnchors(object.list = pbmc, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)


new_name = paste(file_name_no_extension, "_01_integrated.rds", sep  ="_")
saveRDS(immune.combined, file.path(results_folder,new_name))

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

new_name = paste(file_name_no_extension, "_02_integrated_dimred.rds", sep  ="_")
saveRDS(immune.combined, file.path(results_folder,new_name))
print("Completed!")

