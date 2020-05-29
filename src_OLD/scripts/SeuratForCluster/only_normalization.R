#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(SingleCellExperiment))
shhh(library(Seurat))
shhh(library(stringr))
shhh(library(dplyr))

options(future.globals.maxSize = 10000 * 1024^2)
args = commandArgs(trailingOnly=TRUE)

#file = file.path("/home/luisas/Desktop/cluster/proj/code/ebola/src/scripts/Analysis/results/seurat_pbmc_rhemac10_merged_aftercellandgeneqc_afterScrublet.rds")
file = file.path(args[1])

# Extract all names from file path
results_folder = dirname(file)
file_name = basename(file)
file_name_no_extension = gsub(pattern = "\\.rds$", "", file_name)

# # Compute number of counts per row ( number of cells in which the gene is expressed)

#
# # Add metadata
#pbmc$cond <- unlist(lapply(pbmc$orig.ident, function(x) unlist(str_split(x, "-")[[1]][2])))
# table(pbmc$cond)
 # pbmc$hour <- unlist(lapply(pbmc$orig.ident, function(x) rev(str_split(x, "-")[[1]])[1]))
# table(pbmc$hour)
 # pbmc$sample <- unlist(lapply(pbmc$orig.ident, function(x) paste(unlist(str_split(x, "-")[[1]][2]),rev(str_split(x, "-")[[1]])[1])))
#
 #pbmc <- NormalizeData(pbmc)
 #pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
 #all.genes <- rownames(pbmc)
 #pbmc <- ScaleData(pbmc, features = all.genes)
 
 #print("Normalization completed!")
#new_name = paste(file_name_no_extension, "normalization_only_standard.rds", sep  ="_")
#saveRDS(pbmc, file.path(results_folder,new_name))
immune.combined<- readRDS(file.path(results_folder, paste(file_name_no_extension, "normalization_only_standard.rds", sep  ="_")))
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

new_name = paste(file_name_no_extension, "dim_red.rds", sep  ="_")
saveRDS(immune.combined, file.path(results_folder,new_name))
print("Completed!")
