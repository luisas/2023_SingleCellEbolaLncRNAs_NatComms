#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(GenomicRanges))
shhh(library(rtracklayer))
shhh(library(reshape2))
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

# Prep Ref
ref <- import("/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames_addedgenesNovel.gtf")
ref[substr(ref$gene_id,1,3)=="MST", ]$gene_biotype <- "lncRNA"
ref_genes <- ref[ref$type == "gene", ]

# Check which are lnc and which pc 
all_lnc <- (unique(ref[!is.na(ref$gene_biotype) & ref$gene_biotype == "lncRNA",]$gene_name))
all_pc <- (unique(ref[!is.na(ref$gene_biotype) & ref$gene_biotype == "protein_coding",]$gene_name))

# Check for DE status
de_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_stages.rds")
de_lnc <- all_lnc[all_lnc %in% de_genes$primerid]
de_pc <- all_pc[all_pc %in% de_genes$primerid]

expand.grid.unique <- function(x, y, include.equals=FALSE){
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}


f <- function(x,y){
  return(GenomicRanges::distance(ref_genes[ref_genes$gene_name == x, ], ref_genes[ref_genes$gene_name == y, ], ignore.strand = TRUE))
}
f_v <- Vectorize(f)

length(unique(all_lnc))
print("------------------")
print(de_pc)
# Create all the pairs
pairs <- (expand.grid.unique(all_lnc, de_pc))
pairs_df <- as.list(data.frame(t(pairs), stringsAsFactors = F))

# Compute the distamces
distances <- mclapply(pairs_df, function(pair) f_v(pair[1],pair[2]), mc.cores = 16)

res <- cbind(pairs, unname(unlist(distances)))
matrix_result <- acast(as.data.frame(res), V1~V2, value.var="V3")

outfile <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/07_colocation/colocation_all.rds"
outfile2 <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/07_colocation/colocation_all_df.rds"


dir.create(dirname(file.path(outfile)), showWarnings = F)
saveRDS(matrix_result, outfile)
saveRDS(res, outfile2)
