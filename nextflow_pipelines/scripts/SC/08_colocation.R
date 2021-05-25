#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(GenomicRanges))
shhh(library(rtracklayer))
shhh(library(reshape2))
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

ref <- import("/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf")
de_lnc<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_lnc.rds")
de_all_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_stages.rds")
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
# Only check colocation of DE genes
gene_ranges <- ref_genes[ref_genes$gene_name %in% de_all_genes$primerid, ]
genes <- gene_ranges$gene_name

# ADD Missing entries
missing_gene_lines <- unique(de_all_genes$primerid)[!(unique(de_all_genes$primerid) %in% genes)]
gr_list <- lapply(missing_gene_lines, function(x) get_gene_range(x, ref))
missing_gene_ranges <- do.call(base::c,gr_list)
# Add missing gene entries back to reference
ref_genes <-c(gene_ranges, missing_gene_ranges)
genes <- ref_genes$gene_id
#genes <- genes[1:10]
#genes <- unique(ref_genes$gene_id)
lnc <- ref_genes[ref_genes$gene_name %in% de_lnc, ]$gene_id



f <- function(x,y){
  return(GenomicRanges::distance(ref_genes[ref_genes$gene_id == x, ], ref_genes[ref_genes$gene_id == y, ], ignore.strand = TRUE))
}
f_v <- Vectorize(f)


# Create all the pairs
pairs <- (expand.grid.unique(lnc, genes))
pairs_df <- as.list(data.frame(t(pairs), stringsAsFactors = F))

# Compute the distamces
distances <- mclapply(pairs_df, function(pair) f_v(pair[1],pair[2]), mc.cores = 48)


res <- cbind(pairs, unname(unlist(distances)))
matrix_result <- acast(as.data.frame(res), V1~V2, value.var="V3")


outfile <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/07_colocation/colocation.rds"
outfile2 <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/07_colocation/colocation_df.rds"

#outfile <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/07_colocation/colocation_lncvsDE.rds"
#outfile2 <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/07_colocation/colocation_lncvsDE_df.rds"
dir.create(dirname(file.path(outfile)), showWarnings = F)
saveRDS(matrix_result, outfile)
saveRDS(res, outfile2)
