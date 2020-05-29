#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(SingleCellExperiment))
shhh(library(Seurat))
shhh(library(MAST))
shhh(library(data.table))
shhh(library(stringr))
shhh(library(zinbwave))
shhh(library(DESeq2))
options(future.globals.maxSize = 10000 * 1024^2)

args = commandArgs(trailingOnly=TRUE)

file = as.character(args[1])

# Read in input files
robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/RObjects/"
immune.combined  <- readRDS(paste0(robjectsdir, file ))
ebola_genes <- readRDS(paste0(robjectsdir, "ebola_genes.rds"))
all_lncrnas <- readRDS(paste0(robjectsdir, "all_lncrnas.rds"))
annotated_mrnas <- readRDS(paste0(robjectsdir, "annotated_mrnas.rds"))


de_condition_celltype_weighted <- function(subset){
    subset$group <- Idents(subset)
    sc <- as.SingleCellExperiment(subset)
    filter <- rowSums(assay(sc)>0)>0
    sc <- sc[filter, ]
    filter <- colSums(assay(sc)>0)>0
    sc <- sc[,filter]
    counts(sc) <- as.matrix(counts(sc))
    zinb <- zinbwave(sc,K= 0, epsilon = 1000, observationalWeights = TRUE)
    print("Got zinb distribution")
    weights <- assay(zinb, "weights")
    print("calcualted the weigths")
    # Compute p-calues and log fold changes 
    dds_zinb <- DESeqDataSet(zinb, design= ~ nFeature_RNA + individual + group)
    dds_zinb <- DESeq(dds_zinb, sfType="poscounts", useT=TRUE, minmu=1e-6)
    res_zinb <- lfcShrink(dds_zinb, contrast=c("group", levels(Idents(subset))),type = "normal")
    res_zinb$comparison <- paste0(unique(subset$cond), collapse ="-")
    return(res_zinb)
}

de_deseq <- function(myeloids_live_24){
    counts <- myeloids_live_24@assays$RNA@counts +1
    myeloids_live_24$group <- Idents(myeloids_live_24)
    coldata <-myeloids_live_24@meta.data
    # Set up design 
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = coldata,
                                  design= ~ nFeature_RNA + individual + group)
    # Run DE
    dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
    
    return(dds)
}

de_condition_celltype <- function(subset){
    subset$group <- Idents(subset)
    dds <- de_deseq(subset)
    ress <- lfcShrink(dds, contrast=c("group", levels(Idents(subset))), type = "normal")
    ress$comparison <- paste0(unique(subset$cond), collapse ="-")
    return(ress)
}

calc_destats_per_celltype <- function(immune.combined, features, cell_type, weighted = TRUE){
    live <- paste(cell_type,"live", sep ="_")
    media <- paste(cell_type,"media", sep="_")
    irrad <- paste(cell_type, "irrad", sep="_")
    
    # media irrad 
    subset <- subset(immune.combined, features = features, idents = c(media,irrad))
    if(weighted){
        media_irrad <- de_condition_celltype_weighted(subset)
    }else{
        media_irrad <- de_condition_celltype(subset)
    }
    print("Done media irrad")
    
    # irrad - live
    subset <- subset(immune.combined, features = features, idents = c(irrad,live))
    if(weighted){
        irrad_live <- de_condition_celltype_weighted(subset)
    }else{
        irrad_live <- de_condition_celltype(subset)
    }
    print("Done  irrad live")
    # media - live
    subset <- subset(immune.combined, features = features, idents = c(media,live))
   
    if(weighted){
        media_live <- de_condition_celltype_weighted(subset)
    }else{
        media_live <- de_condition_celltype(subset)
    }
    print("Done  media live")
    df <- rbind(irrad_live,media_irrad, media_live )
    df$celltype <- cell_type
    return(df)
}

immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$cond, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
table(Idents(immune.combined))

all_genes <-setdiff(rownames(immune.combined), c(ebola_genes, "ribodepl-MSTRG.72510-unkown"))
immune.combined <- NormalizeData(immune.combined)

print("Running weighted DESeq")
# Run the weighted one
deseq_celltypes_lnc <-rbindlist(lapply(unique(immune.combined$celltype), function(celltype) calc_destats_per_celltype(immune.combined, all_genes, celltype, weighted = TRUE)))
saveRDS(deseq_celltypes_lnc, paste0(robjectsdir,"DESEQWeighted-celltypes-all_genes_noebola.rds"))
print("Done with weighted DESeq")
# Run the unweighted one
deseq_celltypes_lnc <-rbindlist(lapply(unique(immune.combined$celltype), function(celltype) calc_destats_per_celltype(immune.combined, all_genes, celltype, weighted = FALSE)))
saveRDS(deseq_celltypes_lnc, paste0(robjectsdir,"DESEQNOTweighted-celltypes-all_genes_noebola.rds"))

