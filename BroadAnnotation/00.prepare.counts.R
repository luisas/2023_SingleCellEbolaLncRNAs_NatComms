#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load libraries ####
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(dplyr))
suppressMessages(library(reshape))
suppressMessages(library(stringr))
#suppressMessages(library(variancePartition))

options(stringsAsFactors = F)
#tissue <- args[1]
#tissue <- "Lymph node"
tissue <- args[1]
quantification_dir <- "/home/luisas/Desktop/cluster/data/99_BroadAnnotation_Feb2021/06_quantification/"
# ---------------------------
#       Pre-selection 
# ---------------------------

iterate_files <- function(inpath, pattern_string){
  files <- list.files(path=inpath, pattern= pattern_string, full.names=TRUE, recursive=TRUE)
  return(files)
}

# 0. Select Tissue
outpath <- paste0("/home/luisas/Desktop/cluster/data/99_BroadAnnotation_Feb2021/05_DEA/Tissues/",tissue,"/")
dir.create(outpath, recursive = T,  showWarnings = F)



# 1. Read in METADATA summary and only select quantification files of one tissue
metadata_complete <- readRDS("/home/luisas/Desktop/cluster/data/99_BroadAnnotation_Feb2021/metadata_full.rds")


# Update metadata : separate days post infection into stages
metadata_complete$stage <- "Baseline"
metadata_complete[metadata_complete$dpi_time_numeric  %in% c(3,4), ]$stage <- "Early"
metadata_complete[metadata_complete$dpi_time_numeric  %in% c(5,6), ]$stage <- "Middle"
metadata_complete[metadata_complete$dpi_time_numeric  %in% c(7,8), ]$stage <- "Late"
metadata_complete$stage  <- factor(metadata_complete$stage, levels = c("Baseline", "Early", "Middle", "Late"))
metadata_complete$subtissue <- unlist(lapply(metadata_complete$biosample,  function(x) str_split(x, "-")[[1]][2]))

# Correct factors 
metadata_complete$infection_status <- factor(metadata_complete$infection_status, levels = c("Not Infected", "Infected"))
metadata_complete$sex <- factor(metadata_complete$sex)
metadata_complete$batch.extraction <- factor(metadata_complete$batch.extraction)
metadata_complete$sample_name <- metadata_complete$X


saveRDS(metadata_complete, "/home/luisas/Desktop/cluster/data/99_BroadAnnotation_Feb2021/metadata_full.rds")


# 2. Parse counts
quantification <- iterate_files(quantification_dir, "*gene_counts.csv$")
patterns <- paste(metadata_complete[metadata_complete$tissue == tissue,]$sample_name, collapse = "|")
quantification_tissue<- (quantification[grepl(patterns, quantification)])


get_df_count <- function(filepath){
  file <- read.csv(filepath, row.names = 1)
  rownames(file) <- unlist(lapply(rownames(file), function(x) strsplit(x,"|", fixed = T )[[1]][1]))
  file$rownames <- rownames(file)
  return(file)
}

count_list <- lapply(quantification_tissue, get_df_count)
counts_df <- merge_recurse(count_list, by = "rownames" )
counts <- as.matrix(counts_df[,-1])
rownames(counts) <- counts_df$rownames 
colnames(counts) <- gsub("_dedup", "", colnames(counts))
saveRDS(counts, file.path(outpath, paste(tissue, "_genecounts.rds", sep = "")))


# 2. Parse TPMs
gtfs <- iterate_files(quantification_dir, "*tsv$")
patterns <- paste(metadata_complete[metadata_complete$tissue == tissue,]$sample_name, collapse = "|")
gtfs_tissue<- (gtfs[grepl(patterns, gtfs)])


# Create Map for Transcript names to Gene names 
abundances <- list()
for(i in 1:length(gtfs_tissue)) {
  file <- readr::read_tsv(gtfs_tissue[i])
  # Sum up values for double entries (https://github.com/gpertea/stringtie/issues/192)
  file <- file %>% dplyr::group_by(`Gene ID`)%>% dplyr::summarise(TPM = sum(TPM))
  colname <- strsplit(basename(gtfs_tissue[i]), "_dedup")[[1]][1]
  df <-  data.frame(file$TPM, rownames =file$`Gene ID` )
  names(df) <- c(colname, "rownames")
  abundances[[i]] <- df
}


abundances_df <- merge_recurse(abundances, by = "rownames" )
tpm <- as.matrix(abundances_df[,-1])
rownames(tpm) <- abundances_df$rownames 
saveRDS(tpm, file.path(outpath, paste(tissue, "_tpms.rds", sep = "")))



