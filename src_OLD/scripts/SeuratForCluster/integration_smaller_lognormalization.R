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
 pbmc <- readRDS(file.path(results_folder, file_name))
#
# # Add metadata
 pbmc$cond <- unlist(lapply(pbmc$orig.ident, function(x) unlist(str_split(x, "-")[[1]][2])))
# table(pbmc$cond)
 pbmc$hour <- unlist(lapply(pbmc$orig.ident, function(x) rev(str_split(x, "-")[[1]])[1]))
# table(pbmc$hour)
 pbmc$sample <- unlist(lapply(pbmc$orig.ident, function(x) paste(unlist(str_split(x, "-")[[1]][2]),rev(str_split(x, "-")[[1]])[1])))
#
# print("Added meatadata!")
# # Select all gene names
# all_genes <- rownames(pbmc)

 # --------------------------------------
 #           Integration
 # --------------------------------------

 pbmc <- SplitObject(pbmc, split.by = "cond")

 pbmc <- lapply(X = pbmc, FUN = function(x) {
   #x <- SCTransform(x, verbose = FALSE)
   x <- NormalizeData(x)
 })


print("LogNorm completed!")
new_name = paste(file_name_no_extension, "afterLogNorm.rds", sep  ="_")
saveRDS(pbmc, file.path(results_folder,new_name))


immune.anchors <- FindIntegrationAnchors(object.list = pbmc, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)


new_name = paste(file_name_no_extension, "after_integration_5k.rds", sep  ="_")
saveRDS(immune.combined, file.path(results_folder,new_name))

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

new_name = paste(file_name_no_extension, "after_integration_5k_dimred.rds", sep  ="_")
saveRDS(immune.combined, file.path(results_folder,new_name))
print("Completed!")
