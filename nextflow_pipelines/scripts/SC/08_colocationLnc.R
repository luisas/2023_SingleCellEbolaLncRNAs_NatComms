#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(GenomicRanges))
shhh(library(rtracklayer))
shhh(library(reshape2))
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

ref <- import("/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf")
de_all_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_stages.rds")

all_pc <- (unique(ref[!is.na(ref$gene_biotype) & ref$gene_biotype == "protein_coding",]$gene_name))
de_pc <- all_pc[all_pc %in% de_all_genes$primerid]


all_lncrnas <- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/all_lncrnas.rds")
#genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/genes_isg.rds")
#print(de_all_genes)

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

get_gene_range <- function(missing_gene_line, lnc_ranges){
  ranges_exons <- lnc_ranges[lnc_ranges$gene_name == missing_gene_line]
  new_start <- min(start(ranges(ranges_exons))); new_stop <- max(end(ranges(ranges_exons)))
  # Create New Gene Entry
  new_entry<- ranges_exons[1]
  new_entry$type <- "gene"
  ranges(new_entry) <- IRanges(new_start, width =new_stop-new_start+1)
  new_entry$transcript_id <- NA; new_entry$exon_id <- NA; new_entry$transcript_name <-NA
  return(new_entry)
}

# Only keep gene entries
ref_genes <- ref[ref$type == "gene", ]
lnc <- ref_genes[ref_genes$gene_name %in% all_lncrnas, ]$gene_name
# Only check colocation of DE genes
gene_ranges <- ref_genes[ref_genes$gene_name %in% de_all_genes$primerid, ]
genes <- gene_ranges$gene_name

# ADD Missing entries
missing_gene_lines <- unique(de_all_genes$primerid)[!(unique(de_all_genes$primerid) %in% genes)]
gr_list <- lapply(missing_gene_lines, function(x) get_gene_range(x, ref))
missing_gene_ranges <- do.call(base::c,gr_list)
# Add missing gene entries back to reference
ref_genes <-c(gene_ranges, missing_gene_ranges)

#genes <- genes[1:10]
#genes <- unique(ref_genes$gene_id)

f <- function(x,y){
  print(x)
  print("--")
  print(y)
  print("****")
  return(GenomicRanges::distance(ref_genes[ref_genes$gene_name == x, ], ref_genes[ref_genes$gene_name == y, ], ignore.strand = TRUE))
}
f_v <- Vectorize(f)


# Create all the pairs
pairs <- (expand.grid.unique(lnc, de_pc))
pairs_df <- as.list(data.frame(t(pairs), stringsAsFactors = F))
print(length(unique(pairs_df)))

# Compute the distamces
distances <- mclapply(pairs_df, function(pair) f_v(pair[1],pair[2]), mc.cores = 4)

res <- cbind(pairs, unname(unlist(distances)))
print(head(as.data.frame(res)))
matrix_result <- acast(as.data.frame(res), V1~V2, value.var="V3")
print(head(matrix_result))

outfile <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/07_colocation/colocationLnc.rds"
outfile2 <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/07_colocation/colocationLnc_df.rds"

saveRDS(matrix_result, outfile)
saveRDS(res, outfile2)
