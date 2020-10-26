
library(rtracklayer)
library(dplyr)
library(tximport)
library(readr)
library(matrixStats)


# Read command line
args = commandArgs(trailingOnly=TRUE)
dir_counts_ref <- args[1]
novel_concordant <- import(args[2])
outfile <- args[3]

# ----------------------
#   Get Expression 
# ----------------------
iterate_files <- function(inpath, pattern_string){
  files <- list.files(path=inpath, pattern= pattern_string, full.names=TRUE, recursive=TRUE)
  return(files)
} 

#dir_counts_ref <- "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/02_quantification_for_filtering/"
#novel_concordant <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/00_all_novels/novel_rhemac10_concordant_ribodepleted.gtf")
abundance_files <- iterate_files(dir_counts_ref, "t_data.ctab")
# Create Map for Transcript names to Gene names 
tmp <- read_tsv(abundance_files[1])
tx2gene <- tmp[,c("t_name", "gene_id")]
txi <- tximport(abundance_files, type = "stringtie", tx2gene = tx2gene)

# ---------------------------
#   get the expressed ones from the concordant set
# ---------------------------
# Remove non expressed genes 
expression <- txi$abundance
mask<- rowSums(as.data.frame(log(expression)) > 1) > 2
expression <- expression[mask,]
novel_expressed_ribodepleted <- novel_concordant[novel_concordant$gene_id %in% rownames(expression)]

# Save 
export(novel_expressed_ribodepleted,outfile)


