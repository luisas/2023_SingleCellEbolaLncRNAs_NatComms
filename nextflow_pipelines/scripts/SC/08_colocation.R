#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(GenomicRanges))
shhh(library(rtracklayer))
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

ref <- import("/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf")
de_all_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjectsOLD/04_DE/de_all_genes.rds")

# Only keep gene entries
ref_genes <- ref[ref$type == "gene", ]

# Only check colocation of DE genes 
genes <- unique(ref_genes$gene_id)
genes <- ref_genes[ref_genes$gene_name %in% de_all_genes, ]$gene_id

f <- function(x,y){
  return(GenomicRanges::distance(ref_genes[ref_genes$gene_id == x, ], ref_genes[ref_genes$gene_id == y, ], ignore.strand = TRUE))
}
f_v <- Vectorize(f)

res <- (outer(genes[1:20], genes[1:20], function(x,y) f_v(x,y) ))
rownames(res) <- colnames(res) <- genes[1:20]
saveRDS(res, "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjectsOLD/07_colocation/colocation_small.rds")


res <- (outer(genes, genes, function(x,y) f_v(x,y) ))

rownames(res) <- colnames(res) <- genes
saveRDS(res, "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjectsOLD/07_colocation/colocation.rds")



