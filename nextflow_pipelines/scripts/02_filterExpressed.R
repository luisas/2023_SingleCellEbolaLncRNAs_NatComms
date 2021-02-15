
library(rtracklayer)
library(dplyr)
library(tximport)
library(readr)
library(matrixStats)

# Read command line
args = commandArgs(trailingOnly=TRUE)
dir_counts_ref <- args[1]
ref <- import(args[2])
outfile <- args[3]

# ----------------------
#   Get Expression 
# ----------------------
iterate_files <- function(inpath, pattern_string){
  files <- list.files(path=inpath, pattern= pattern_string, full.names=TRUE, recursive=TRUE)
  return(files)
} 

#dir_counts_ref <- "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/01_quantification_for_filtering/"
#ref <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/00_prefilter_candidates/prefilter_candidates.gtf")

abundance_files <- iterate_files(dir_counts_ref, "*.tsv")

# Read all files
abundances <- list()
for(i in 1:length(abundance_files)) {
  file <- readr::read_tsv(abundance_files[i])
  # Sum up values for double entries (https://github.com/gpertea/stringtie/issues/192)
  file <- file %>% dplyr::group_by(`Gene ID`)%>% dplyr::summarise(TPM = sum(TPM))
  abundances[[i]] <- data.frame(file$TPM, row.names =file$`Gene ID` )
}

# Summarize all TPMs from all quantification files 
rn <- rownames(abundances[[1]])
dat <- abundances[[1]]
for(i in 2:length(abundances)) {
  dat <- merge(dat, abundances[[i]],  by= "row.names", all.x= F, all.y= F) [,-1]
  rownames(dat) <- rn
}



# ---------------------------
#  Only retain genes expressed at log(TPM) > 0.5 in at least 3 samples. 
# ---------------------------
expression <- dat

mask<- rowSums(as.data.frame(log(expression) > 0.5)) > 2
table(mask)
expression <- expression[mask,]
novel_expressed <- ref[ref$gene_id %in% rownames(expression)]


# Save 
export(novel_expressed,outfile)


