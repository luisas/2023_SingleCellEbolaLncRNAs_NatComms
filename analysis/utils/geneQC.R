
print("hola")
shhh <- suppressPackageStartupMessages
shhh(library(SingleCellExperiment))
shhh(library(Seurat))

args = commandArgs(trailingOnly=TRUE)

#file = file.path("/home/luisas/Desktop/cluster/proj/code/ebola/src/scripts/Analysis/results/seurat_pbmc_rhemac10_merged_aftercellqc.rds")
file = file.path(args[1])

print(file)
# Extract all names from file path
results_folder = dirname(file)
file_name = basename(file)
file_name_no_extension = gsub(pattern = "\\.rds$", "", file_name)


# Compute number of counts per row ( number of cells in which the gene is expressed)
pbmc <- readRDS(file)
print("Here1")
pbmc


n_cells <- rowSums(counts(pbmc) > 0) 

#saveRDS(n_cells, file.path(results_folder,"number_of_cellss.rds"))
print("Here2")
# filter out the genes showing clusters with less than 10 cells
pbmc <- pbmc[n_cells > 10, ]
new_name = paste(file_name_no_extension, "aftercellandgeneqc.rds", sep  ="_")
saveRDS(pbmc, file.path(results_folder,new_name))

print("Completed!")

