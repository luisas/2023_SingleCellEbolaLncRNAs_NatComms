#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
 


options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

# Read in Seurat Object 
#datadir <- "/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/"
datadir <- args[1]
file <- args[2]
ident <- args[3]
immune.combined <- readRDS(file.path(datadir,file))
#Filename where to store seurat object 
output <- file.path(datadir,"03_prep/02_immune.combined.infectionstatus_noebola.rds")
output_dir_subsets <- dirname(dirname(output))
#ident <- "Monocyte"
viralload_summary <-readRDS(file.path(datadir,"03_prep/df_viralpercentage.rds"))
  
# Identify Ebola Genes
ebola_genes <- readRDS(file.path(datadir,"05_stats/ebola_genes.rds"))


# I need to add infection status 
immune.combined$infection <- viralload_summary[colnames(immune.combined),]$classification
immune.combined$viraload <- viralload_summary[colnames(immune.combined),]$percentage_viral_reads
saveRDS(immune.combined, file.path(datadir,"03_prep/03_immune.combined.infectionstatus_neut.rds"))


# Remove ebola genes 
immune.combined_noebola <- immune.combined[setdiff(rownames(immune.combined), c(ebola_genes)),]

# Re-normalize
DefaultAssay(immune.combined_noebola) <- "RNA"
immune.combined_noebola <- NormalizeData(immune.combined_noebola)
immune.combined_noebola <- FindVariableFeatures(immune.combined_noebola)
immune.combined_noebola <- ScaleData(immune.combined_noebola)
#print(Idents(immune.combined_noebola))
#print(Idents(immune.combined))
#Idents(immune.combined_noebola) <- Idents(immune.combined)

# Save
saveRDS(immune.combined_noebola, output)

immune.combined_noebola <- readRDS(output)


# Prepare other subsets 
myeloids <- subset(immune.combined_noebola, ident = ident)
infected_myeloids <- myeloids[,myeloids$infection == "Infected" ]

if("group" %in% colnames(infected_myeloids@meta.data)){
  infected_myeloids_late <- infected_myeloids[,infected_myeloids$group == "late" ]
}else{
  infected_myeloids_late <- infected_myeloids[,infected_myeloids$dpi == "H024" ]
}

saveRDS(infected_myeloids_late, file.path(datadir,"03_prep/03_monocytes_infected_late.rds"))







