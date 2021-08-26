#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
library(gtools)
library(parallel)

options(future.globals.maxSize = 10000 * 1024^2)
args = commandArgs(trailingOnly=TRUE)



file <- args[1]
genes <- args[2]
out <- args[3]
ebola_genome_percentage_df <- readRDS(args[4])
out2 <- args[5]

# 0. read in file 
immune.combined <- readRDS(file)
DefaultAssay(immune.combined) <- "RNA"
ebola_genes <- readRDS(genes)

# 1. remove ebola genes 
immune.combined_ebolaremoved <- immune.combined[setdiff(rownames(immune.combined), ebola_genes),]
immune.combined_ebolaremoved_norm <- NormalizeData(immune.combined_ebolaremoved)
print("Removed ebola reads")
saveRDS(immune.combined_ebolaremoved_norm, file.path(out))

# 2. Add infection information
print("Infection information matches information in seurat object")
print(table(colnames(immune.combined) == rownames(ebola_genome_percentage_df)))
immune.combined_ebolaremoved_norm$infection <- ebola_genome_percentage_df$classification 
immune.combined_ebolaremoved_norm$viral_load <- ebola_genome_percentage_df$percentage_viral_reads 
saveRDS(immune.combined_ebolaremoved_norm, file.path(out2))








