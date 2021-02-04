
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

#dir_counts_ref <- "/home/luisas/Desktop/cluster/data/99_BroadAnnotation/03_novel_lncRNAs_list/01_quantification_for_filtering/"
#novel_concordant <- import("/home/luisas/Desktop/cluster/data/99_BroadAnnotation/03_novel_lncRNAs_list/00_prefilter_candidates/prefilter_candidates.gtf")

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

# Create Map for Transcript names to Gene names 
#tmp <- read_tsv(abundance_files[1])
#tx2gene <- tmp[,c("t_name", "gene_id")]
#txi <- tximport(abundance_files, type = "stringtie", tx2gene = tx2gene)

a <- novel_concordant[!(novel_concordant$gene_id %in% rownames(expression))]
length(unique(a$gene_id))
# ---------------------------
#   Get the expressed ones from the concordant set
# ---------------------------
# Remove non expressed genes 
#expression <- txi$abundance
expression <- dat
expression_new <- expression[rownames(expression) %in% novel_concordant$gene_id,]
nrow(expression_new)

empty <- expression[rowSums(expression) == 0,]
length(unique(rownames(empty)))
mask<- rowSums(as.data.frame(expression > 0.0)) > 0
table(mask)
expression <- expression[mask,]
novel_expressed_ribodepleted <- novel_concordant[novel_concordant$gene_id %in% rownames(expression)]
length(unique(novel_expressed_ribodepleted$transcript_id))


# Save 
export(novel_expressed_ribodepleted,outfile)


