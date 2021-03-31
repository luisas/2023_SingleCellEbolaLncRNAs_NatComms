#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load libraries ####
suppressMessages(library(edgeR))
suppressMessages(library(limma))
#suppressMessages(library(variancePartition))

options(stringsAsFactors = F)
tissue <- args[1]
#tissue <- "Vagina"
outpath <- paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/03_DEA/Tissues/",tissue,"/")
dir.create(outpath, recursive = T)
  
# Function ####
limma_lm <- function(fit, covariate, covariate_data){
  v.contrast <- rep(0,ncol(fit$design))
  if(is.factor(covariate_data[,covariate])){
    if(covariate == "Sex"){
      if(table(covariate_data$Sex)["1"] > 2 & table(covariate_data$Sex)["2"] > 2 ){      #if( table(covariate_data$Gender)[1]/sum(table(covariate_data$Gender)) > 0.1 & table(covariate_data$Gender)[1]/sum(table(covariate_data$Gender)) < 0.9){
        v.contrast[ which( colnames(fit$design) == "Sex2") ] <- 1 # contrast againts basal level which is male
        contrast.matrix <- cbind( "C1" = v.contrast)
      }else{
        return(NA)
      }
    }else if(covariate=="Ancestry"){
      #EUR_AFR <-  sum(table(covariate_data$Ancestry)[c("EUR","AFR")])
      if( table(covariate_data$Ancestry)["AFR"] > 2 & table(covariate_data$Ancestry)["EUR"] > 2 ){
        #if( table(covariate_data$Ancestry)["EUR"]/EUR_AFR > 0.1 & table(covariate_data$Ancestr)["EUR"]/EUR_AFR < 0.9){
        v.contrast[ which( colnames(fit$design) == "AncestryEUR") ] <- 1 # contrast againts basal level which is male
        contrast.matrix <- cbind( "C1" = v.contrast)
      }else{
        return(NA)
      }
    # }else if( covariate=="Age_Category"){
    #   #young_old <- sum(table(covariate_data$Age_Category)[c("1","2")])
    #   if( table(covariate_data$Age_Category)["1"] > 2 & table(covariate_data$Age_Category)["2"] > 2){
    #     #if( table(covariate_data$Age_Category)["1"]/young_old > 0.1 & table(covariate_data$Age_Category)["1"]/young_old < 0.9){
    #     v.contrast[ which( colnames(fit$design) == "Age_Category2") ] <- 1 # contrast againts basal level which is male
    #     contrast.matrix <- cbind( "C1" = v.contrast)
    #   }else{
    #     return(NA)
    #   }
    # }else if( covariate=="BMI_Category"){
    #   #young_old <- sum(table(covariate_data$Age_Category)[c("1","2")])
    #   if( table(covariate_data$BMI_Category)["1"] > 9 & table(covariate_data$BMI_Category)["2"] > 9){
    #     #if( table(covariate_data$Age_Category)["1"]/young_old > 0.1 & table(covariate_data$Age_Category)["1"]/young_old < 0.9){
    #     v.contrast[ which( colnames(fit$design) == "BMI_Category2") ] <- 1 # contrast againts basal level which is male
    #     contrast.matrix <- cbind( "C1" = v.contrast)
    #   }else{
    #     return(NA)
    #   }
    }else{
      return(NULL)
    }
  }else{ # Continuous variable
    v.contrast[ which(colnames(fit$design) == covariate)] <- 1 #
    contrast.matrix <- cbind( "C1" = v.contrast)
  }
  fitConstrasts <- contrasts.fit(fit,contrast.matrix)
  eb = eBayes(fitConstrasts)
  tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
  return(tt.smart.sv)
}  


# 1.Reading in count-data ####
counts <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/00_Data/Tissues/",tissue,"/",tissue,".SelectedSamples.counts.PC_lincRNA.rds"))
#counts <- readRDS(paste0("~/GTEx_v8/Raquel/00_Data/Tissues/",tissue,"/",tissue,".SelectedSamples.counts.PC_lincRNA.rds"))
tpm <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/00_Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.PC_lincRNA.rds"))
#tpm <- readRDS(paste0("~/GTEx_v8/Raquel/00_Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.PC_lincRNA.rds"))
metadata <- readRDS( paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/00_Data/Tissues/",
                             tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
#metadata <- readRDS( paste0("~/GTEx_v8/Raquel/00_Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
tissue_info <- readRDS("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/00_Data/Tissue_Info.rds")
#tissue_info <- readRDS("~/GTEx_v8/Raquel/00_Data/Tissue_Info.rds")
tissue_id <- tissue_info[tissue_info$tissue_ID==tissue,]$tissue_id

# 2. Remove genes that are lowly expressed ###
# For good results, the counts matrix should be filtered to remove remove rows with very low counts before running voom().
# Genes expressed per tissue
# 1. TPM>=0.1 in at least 20% of the tissue samples
exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=0.1) ) >= 0.2*ncol(tpm)  ]  

# 2. Count >=6 in at least 20% of the tissue samples
exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=6) ) >= 0.2*ncol(counts)  ]  

# 3. Intersect gene lists
exprs_genes <- intersect(exprs_genes.tpm,
                         exprs_genes.counts) 

# Exclude chrY genes in female-only tissues
# Gene Annotation
gene_annotation  <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.bed", header = F)
#gene_annotation  <- read.delim("~/bsc83_Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.bed", header = F)
colnames(gene_annotation) <- c("chr","start","end","strand","feature","Ensembl_ID","gene_name",
                               "Biotype","source")
# Select protein-coding and lincRNA
gene_annotation <- gene_annotation[gene_annotation$Biotype %in% c("protein_coding","lincRNA"),]
# Remove pseudo-autosomal genes
PAR_genes <- sapply(grep(".Y", gene_annotation$Ensembl_ID, value = T), function(gene)
  unlist(strsplit(gene, split = "_"))[[1]]
)
gene_annotation <- gene_annotation[-unlist(lapply(PAR_genes, function(gene)
  grep(gene, gene_annotation$Ensembl_ID))),]
# 26,676 genes
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$Ensembl_ID
if(tissue %in% c("Uterus","Ovary","Vagina","BreastMammaryTissue_Female")){
  exprs_genes <- exprs_genes[!exprs_genes %in% Y_genes]
}

# Save gene list
saveRDS(exprs_genes,
        paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/00_Data/Tissues/",tissue,"/",tissue,".SelectedSamples.Expressed_genes.rds"))

# Variable to hold objects
my_data <- list()
# 3. Normalising gene expression distributions
library(edgeR)

# Create DGEList object
dge <- DGEList(counts[exprs_genes,])

# Calculate normalization factors (does not do the normalization, only computes the factors)
dge <- calcNormFactors(dge)

# 4. voom + limma. Model expression ~ PEERs
library(limma)
# Voom
v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) # samples are treated as replicates

# Expression ~ covariates + traits
covariates <- c("HardyScale", "IschemicTime", "RIN", "ExonicRate","PEER1","PEER2")
if(! tissue %in% c("Vagina","Uterus","Ovary","Prostate","Testis","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
  individual_traits <- c("Age","Ancestry", "Sex","BMI")
}else{
  individual_traits <- c("Age","Ancestry", "BMI")
}
fml_args_mod <- paste(c(covariates, individual_traits), collapse = " + ")
mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)

# Limma fit
fit <- lmFit(v, mod)

# Add objects to data
my_data[["dge"]] <- dge
my_data[["v"]] <- v
my_data[["fit"]] <- fit

# Limma test
if(tissue %in% c("Vagina","Uterus","Ovary","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
  dea_res <- list()
  dea_res[[1]] <- NA
  dea_res <- c(dea_res, 
               lapply(individual_traits, function(phenotype) limma_lm(fit,phenotype,metadata ) ))
  names(dea_res) <- c("Sex",individual_traits)
  dea_res <- dea_res[c("Age","Ancestry", "Sex","BMI")]
}else{
  dea_res <- lapply(individual_traits, function(phenotype) limma_lm(fit,phenotype,metadata ) )
  names(dea_res) <- individual_traits
}


# Save table with results
saveRDS(dea_res,
        paste0(outpath,tissue,".voom_limma.covariates_and_traits.results.rds"))

# Expression ~ covariates + traits + PEER1 ####
#fml_args_mod <- paste0("+ ",c(covariates, "PEER1", individual_traits ))
#fml_args_mod <- paste(c(covariates, "PEER1", individual_traits), collapse = " + ")
#mod <- model.matrix( as.formula(paste(" ~ 1 ", paste(fml_args_mod,collapse = ""))), data =  metadata)
#mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)

# Limma fit
#fit <- lmFit(v, mod)

# Limma test
#if(tissue %in% c("Vagina","Uterus","Ovary","Testis","Prostate")){
#  dea_res <- list()
#  dea_res[[1]] <- NA
#  dea_res <- c(dea_res, 
#               lapply(individual_traits, function(phenotype) limma_lm(fit,phenotype,metadata ) ))
#  names(dea_res) <- c("Sex",individual_traits)
#  dea_res <- dea_res[c("Age","Ancestry", "Sex","BMI")]
#}else{
#  dea_res <- lapply(individual_traits, function(phenotype) limma_lm(fit,phenotype,metadata ) )
#  names(dea_res) <- individual_traits
#}

#my_data[["fit_PEER1"]] <- fit

# Save tables with results
#saveRDS(dea_res,
#        paste0(outpath,tissue,".voom_limma.covariates_and_traits_and_PEER1.results.rds"))

# Sade objects
saveRDS(my_data,
        paste0(outpath,tissue,".voom_limma.data_objects.rds"))


