#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(SingleCellExperiment))
shhh(library(Seurat))
shhh(library(stringr))
shhh(library(dplyr))
shhh(library(MAST))
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

#file = file.path("/home/luisas/Desktop/cluster/proj/code/ebola/src/scripts/Analysis/results/seurat_pbmc_rhemac10_merged_aftercellandgeneqc_afterScrublet.rds")
file = file.path(args[1])
robjectsdir  = file.path(args[2])
dir.create(robjectsdir, recursive = T,  showWarnings = F)


immune.combined <- readRDS(file)
# -----------------------
#   General MAST model 
# -----------------------

get_mast_model  <- function(subset, relevel = "Not Infected", contarst = "conditionInfected"){
  # Convert to single cell experiment
  subset <- as.SingleCellExperiment(subset)
  # Extract ident column and use it as level
  cond<-factor(colData(subset)$ident)
  cond<-relevel(cond,relevel)
  colData(subset)$cond<-cond
  colData(subset)$wellKey<-colnames(subset)
  
  # Prepare SingleCellAssay Object
  scam <- as(subset, 'SingleCellAssay')
  rowData(scam)$primerid <- rownames(rowData(scam))
  
  # Calculate number of genes per cell
  cdr2 <-colSums(assay(scam)>0)
  colData(scam)$cngeneson <- scale(cdr2)
  
  # We’ll fit a hurdle model
  # modeling the condition and (centered) ngeneson factor + fresh frozen, thus adjusting for the cellular detection rate.
  # In order to have more interpretable coefficients, we’ll set the reference level of the factor to be the “baseline” cells.
  cond<-factor(colData(scam)$cond)
  cond<-relevel(cond,relevel)
  colData(scam)$condition<-cond
  
  # Set Up the model 
  zlmCond <- zlm(~condition + cngeneson , scam)
  
  # Summarize data and get the contrast 
  summaryCond <- summary(zlmCond, doLRT=contarst) 
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast==contarst & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast==contarst & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  return(fcHurdle)
}


# --------------------------------------
#   MAST each celltype stage vs baseline 
# --------------------------------------


myeloids <- subset(immune.combined, idents= "Myeloid")
Idents(myeloids) <-  myeloids$infection
dim(myeloids)
mast_model <- get_mast_model(myeloids, relevel = "Not Infected", contarst = paste0("condition", "Infected"))
saveRDS(mast_model,file.path(robjectsdir, paste0("MAST_exvivo_model_infected_vs_notinfected_myeloids.rds")))


# Only on live Myeloids
myeloid_live <- myeloids[,myeloids$cond == "live"  ]
mast_model <- get_mast_model(myeloid_live, relevel = "Not Infected", contarst = paste0("condition", "Infected"))
saveRDS(mast_model,file.path(robjectsdir, paste0("MAST_exvivo_model_infected_vs_notinfected_myeloids_LIVE.rds")))

# Only on live Myeloids 24h 
myeloid_live_24 <- myeloid_live[,myeloid_live$hpi == "H024"  ]
dim(myeloid_live_24)
mast_model <- get_mast_model(myeloid_live_24, relevel = "Not Infected", contarst = paste0("condition", "Infected"))
saveRDS(mast_model,file.path(robjectsdir, paste0("MAST_exvivo_model_infected_vs_notinfected_myeloids_LIVE_24only.rds")))




