
library(rtracklayer)
library(dplyr)

# Read command line
args = commandArgs(trailingOnly=TRUE)
ribodepl_novel <- import(args[1])
ref_gtf <- import(args[2])
outfile <- args[3]
outfile_name <- args[4]

#ribodepl_novel_compared_to_polya <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/03_polyA_vs_ribodepl/ribodeplVSpolyA.annotated.gtf")
#polyA_novel <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/02_novel_expressed/novel_expressed_polya.gtf")
#ribodepl_novel <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/02_novel_expressed/novel_expressed_ribodepleted.gtf")
#ref_gtf <-  import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/gene_annotations/rheMac10_EBOV-Kikwit.gtf")


# ---------------------------
#   Add Names 
# ---------------------------
add_gene_name <- function(gtf){
  missing_gene_name <- gtf[gtf$type == "exon"][is.na(gtf[gtf$type == "exon"]$gene_name)]                                
  transcripts <- gtf[gtf$type == "transcript"]
  # To the ones that have no gene name, find the gene name associated with the transcript ID and use it
  test <-  lapply(seq_along(missing_gene_name), function(x) transcripts[transcripts$transcript_id == missing_gene_name[x]$transcript_id]$gene_name)
  gtf[gtf$type == "exon"][is.na(gtf[gtf$type == "exon"]$gene_name)]$gene_name <- unlist(test)
  gtf[is.na(gtf$gene_name)]$gene_name <- paste(gtf[is.na(gtf$gene_name)]$gene_id, "-unknown", sep = "")
  return(gtf)
}
add_prefix <- function(gtf, prefix){
  gtf$gene_id <- paste(prefix, gtf$gene_id, sep = "") 
  gtf$transcript_id <- paste(prefix, gtf$transcript_id, sep = "") 
  return(gtf)
}


#gtf_ribodepl_newids <- add_prefix(ribodepl_novel, "ribodepl-")

# ---------------------------------------
# Merge Novel and Ref  together ( just append ) 
# ---------------------------------------
glst_ref_and_novel <-list(ref_gtf,ribodepl_novel)
gtf_ref_and_novel <- do.call(c, as(glst_ref_and_novel, "GRangesList"))


gtf_ref_and_novel_sort <- sort(gtf_ref_and_novel)
export(gtf_ref_and_novel_sort, outfile)  


mask <- is.na(gtf_ref_and_novel_sort$gene_name)
gtf_ref_and_novel_sort[mask]$gene_name <- paste(gtf_ref_and_novel_sort[mask]$gene_id, "-unknown", sep = "") 
mask_t <- is.na(gtf_ref_and_novel_sort$transcript_name)
gtf_ref_and_novel_sort[mask_t]$transcript_name <-paste(gtf_ref_and_novel_sort[mask_t]$transcript_id, "-unknown", sep = "")  
export(gtf_ref_and_novel_sort, outfile_name)  


