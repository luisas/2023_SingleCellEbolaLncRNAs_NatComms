#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(SingleCellExperiment))
shhh(library(Seurat))
shhh(library(MAST))
shhh(library(data.table))
shhh(library(stringr))
options(future.globals.maxSize = 10000 * 1024^2)

args = commandArgs(trailingOnly=TRUE)

file = as.character(args[1])

# Read in input files
robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/RObjects/"
immune.combined  <- readRDS(paste0(robjectsdir, file ))
ebola_genes <- readRDS(paste0(robjectsdir, "ebola_genes.rds"))
all_lncrnas <- readRDS(paste0(robjectsdir, "all_lncrnas.rds"))
annotated_mrnas <- readRDS(paste0(robjectsdir, "annotated_mrnas.rds"))


get_mast_model  <- function(subset, relevel = "bystander", contarst = "conditioninfected"){
    subset <- as.SingleCellExperiment(subset)
    cond<-factor(colData(subset)$ident)
    cond<-relevel(cond,relevel)
    colData(subset)$cond<-cond
    colData(subset)$wellKey<-colnames(subset)
    subset$individual <- unlist(lapply(subset$orig.ident , function(x) paste(unlist(str_split(x, "-")[[1]][1]))))
    scam <- as(subset, 'SingleCellAssay')
    rowData(scam)$primerid <- rownames(rowData(scam))
    cdr2 <-colSums(assay(scam)>0)
    colData(scam)$cngeneson <- scale(cdr2)
    #We’ll fit a hurdle model, modeling the condition and (centered) ngeneson factor, thus adjusting for the cellular detection rate.
    #In order to have more interpretable coefficients, we’ll set the reference level of the factor to be the “unstimulated” cells.
    cond<-factor(colData(scam)$cond)
    cond<-relevel(cond,relevel)
    colData(scam)$condition<-cond
    zlmCond <- zlm(~condition + cngeneson+individual, scam)
    # Summarize data
    summaryCond <- summary(zlmCond, doLRT=contarst)
    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(summaryDt[contrast==contarst & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                          summaryDt[contrast==contarst & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    return(fcHurdle)
}

calc_destats_per_celltype <- function(immune.combined, features, cell_type){
    live <- paste(cell_type,"live", sep ="_")
    media <- paste(cell_type,"media", sep="_")
    irrad <- paste(cell_type, "irrad", sep="_")

    # media live
    subset <- subset(immune.combined, features = features, idents = c(media,live))
    hurdle_media_live <- get_mast_model(subset, relevel = paste0(cell_type, "_media"), contarst = paste0("condition",cell_type,"_live"))
    scam <- as(as.SingleCellExperiment(subset), 'SingleCellAssay')
    rowData(scam)$primerid <- rownames(rowData(scam))
    hurdle_media_live_sig <- merge(hurdle_media_live, as.data.table(mcols(scam)), by='primerid')
    hurdle_media_live_sig$comparison <- "media-live"

    # media irrad
    subset <- subset(immune.combined, features = features, idents = c(media,irrad))
    hurdle_media_irrad <- get_mast_model(subset, relevel = paste0(cell_type, "_media"), contarst = paste0("condition",cell_type,"_irrad"))
    scam <- as(as.SingleCellExperiment(subset), 'SingleCellAssay')
    rowData(scam)$primerid <- rownames(rowData(scam))
    hurdle_media_irrad_sig <- merge(hurdle_media_irrad, as.data.table(mcols(scam)), by='primerid')
    hurdle_media_irrad_sig$comparison <- "media-irrad"

    # irrad - live
    subset <- subset(immune.combined, features = features, idents = c(irrad,live))
    hurdle_irrad_live <- get_mast_model(subset, relevel = paste0(cell_type, "_irrad"), contarst = paste0("condition",cell_type,"_live"))
    scam <- as(as.SingleCellExperiment(subset), 'SingleCellAssay')
    rowData(scam)$primerid <- rownames(rowData(scam))
    hurdle_irrad_live_sig <- merge(hurdle_irrad_live, as.data.table(mcols(scam)), by='primerid')
    hurdle_irrad_live_sig$comparison <- "irrad-live"


    df <- rbind(hurdle_irrad_live_sig,hurdle_media_irrad_sig, hurdle_media_live_sig )
    df$celltype <- cell_type
    #df <- data.frame( comparison = c("MEDIA-LIVE", "MEDIA-IRRAD", "IRRAD-LIVE"), de_genes = c(nrow(hurdle_media_live_sig),nrow(hurdle_media_irrad_sig),nrow(hurdle_irrad_live_sig)))
    return(df)
}

immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$cond, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
table(Idents(immune.combined))


all_genes <-setdiff(rownames(immune.combined), c(ebola_genes, "ribodepl-MSTRG.72510-unkown"))
immune.combined <- NormalizeData(immune.combined)
hurdle_summary_celltypes_all <-rbindlist(lapply(unique(immune.combined$celltype), function(celltype) calc_destats_per_celltype(immune.combined, all_genes, celltype)))
saveRDS(hurdle_summary_celltypes_all, paste0(robjectsdir,"hurdle_summary_celltypes_all_noebola.rds"))
