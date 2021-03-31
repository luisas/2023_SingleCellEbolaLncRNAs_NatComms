#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load libraries ####
suppressMessages(library(edgeR))
suppressMessages(library(limma))
#suppressMessages(library(variancePartition))

options(stringsAsFactors = F)

# 0. Select Tissue
tissue <- c("Lymph node")
lateonly <- TRUE

outpath <- paste0("/home/luisas/Desktop/cluster/data/99_BroadAnnotation_Feb2021/05_DEA/Tissues/",tissue,"/")
dir.create(outpath, recursive = T,  showWarnings = F)

# 1.Reading in count-data
counts <- readRDS(file.path(outpath, paste(tissue, "_genecounts.rds", sep = "")))
tpm <- readRDS(file.path(outpath, paste(tissue, "_tpms.rds", sep = "")))

# 2. Read in metadata
metadata <-  readRDS("/home/luisas/Desktop/cluster/data/99_BroadAnnotation_Feb2021/metadata_full.rds")
# Only keep samples for which i have both counts and metadata
metadata_tissue <- metadata[metadata$tissue %in% tissue, ]

# Only retain samples with more than 3M reads
metadata_tissue <- metadata_tissue[(metadata_tissue$hostReadCount > 3* 10^6),]


# Only keep samples that are n both metadata and counts
metadata_tissue <- metadata_tissue[metadata_tissue$X %in% colnames(counts),]
rownames(metadata_tissue) <- metadata_tissue$X
counts <- counts[, colnames(counts) %in% metadata_tissue$X]
metadata_tissue <- metadata_tissue[colnames(counts),]


# 2. Remove genes that are lowly expressed ###
# For good results, the counts matrix should be filtered to remove remove rows with very low counts before running voom().
# Genes expressed per tissue
# 1. TPM>=0.1 in at least 20% of the tissue samples
exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=1) ) >= 0.6*ncol(tpm)  ]  

# 2. Count >=6 in at least 20% of the tissue samples
exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=5) ) >= 0.2*ncol(counts)  ]  

# 3. Intersect gene lists
exprs_genes <- intersect(exprs_genes.tpm,
                         exprs_genes.counts) 

dim(counts[exprs_genes,])
saveRDS(exprs_genes,paste0(outpath,"/",tissue,".SelectedSamples.Expressed_genes.rds"))


# Variable to hold objects
my_data <- list()


# 3. Normalising gene expression distributions
# Create DGEList object
dge <- DGEList(counts[exprs_genes,])
# Calculate normalization factors (does not do the normalization, only computes the factors)
dge <- calcNormFactors(dge)
saveRDS(dge, file.path("/home/luisas/Desktop/cluster/data/99_BroadAnnotation_Feb2021/05_DEA/Tissues/", tissue, "dge.rds"))


# 4. Voom + Limma

# Set up the model 
# Expression ~ covariates + traits
covariates <- c("sex","weight", "batch.extraction")
covariates <- c("sex", "weight")
traits <- c("stage")

# Check dataframe 
str(metadata_tissue[,c("X",covariates, traits, "id.individual")])


fml_args_mod <- paste(c( traits, covariates), collapse = " + ")
mod <- model.matrix( as.formula(paste(" ~  0+", paste(fml_args_mod,collapse = " "))), data =  metadata_tissue)


# Voom
v_temp <- voom(dge, design = mod, normalize="quantile", save.plot=F, plot = F) 
# Duplicate correlations for controlling for individuals
dupcor <- duplicateCorrelation(v_temp,mod,block=metadata_tissue$id.individual)
v = voom( dge, design = mod, plot=FALSE, block=metadata_tissue$id.individual, correlation=dupcor$consensus)
dupcor <- duplicateCorrelation(v, mod, block=metadata_tissue$id.individual)
# But this step uses only the genome-wide average for the random effect
fit <- lmFit(v, mod, block=metadata_tissue$id.individual, correlation=dupcor$consensus)


colnames(mod) <- gsub("stage", "", colnames(mod))

# Create contrasts 
contr.matrix <- makeContrasts(
  LatevsBaseline = Late-Baseline, 
  MiddlevsBaseline = Middle-Baseline, 
  EarlyvsBaseline = Early-Baseline, 
  levels = colnames(mod))

#v = voom( dge, design = mod , plot=FALSE)
#fit <- lmFit(v, mod)
vfit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(vfit)
summary(decideTests(efit))

topTable(efit)

# Add objects to data
my_data[["dge"]] <- dge
my_data[["v"]] <- v
my_data[["fit"]] <- fit


# Save objects
saveRDS(my_data,
        paste0(outpath,tissue,".voom_limma.data_objects.rds"))


