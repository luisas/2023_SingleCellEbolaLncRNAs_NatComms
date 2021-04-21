#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(SingleCellExperiment))
shhh(library(Seurat))
shhh(library(stringr))
shhh(library(dplyr))
shhh(library(MAST))
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

file = file.path(args[1])
robjectsdir  = file.path(args[2])

dir.create(robjectsdir, showWarnings = FALSE)
immune.combined <- readRDS(file)

# -----------------------
#   General MAST model 
# -----------------------
get_mast_model  <- function(subset, relevel = "bystander", contarst = "conditioninfected"){
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
  zlmCond <- zlm(~condition + cngeneson+ freshfrozen , scam)
  
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
mast_stateVSbaseline_perCelltype <- function(immune.combined , celltype, stages, robjectsdir, prefix = "MAST_invivo_model_"){
  # Create the right subset (Only cells of selectedcelltype)
  subset <- subset(immune.combined, idents = celltype)
  # And cells of the comparison groups
  Idents(subset) <- subset$group
  subset <- subset(subset, idents = stages)
  
  # Get the model (Already set up with right contrast )
  mast_model <- get_mast_model(subset, relevel = stages[1], contarst = paste0("condition", stages[2]))
  saveRDS(mast_model,file.path(robjectsdir, paste0(prefix, celltype,"_", stages[2],".rds")))
  return(mast_model)
}


# Prepare which are the comparisons we are interested in 
celltypes <- unique(Idents(immune.combined))
ref_comparison <- "baseline"
comparisons <- list( "late", "middle", "early")


# Run the comparison 
for( celltype in celltypes){
  for(comparison in comparisons){
    mast_stateVSbaseline_perCelltype(immune.combined, celltype, c(ref_comparison, comparison),robjectsdir, "MAST_invivo_model_")
  }
}


