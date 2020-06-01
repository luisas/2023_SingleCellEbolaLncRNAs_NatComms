
library(rtracklayer)
library(dplyr)

# Read command line
args = commandArgs(trailingOnly=TRUE)
ribodepl_novel_compared_to_polya <- import(args[1])
ribodepl_novel <- import(args[2])
polyA_novel <- import(args[3])
ref_gtf <- import(args[4])
outfile <- args[5]

#ribodepl_novel_compared_to_polya <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list_2/03_polyA_vs_ribodepl/ribodeplVSpolyA.annotated.gtf")
#polyA_novel <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list_2/02_novel_expressed/novel_expressed_polya.gtf")
#ribodepl_novel <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list_2/02_novel_expressed/novel_expressed_ribodepleted.gtf")
#ref_gtf <-  import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/gene_annotations/rheMac10_EBOV-Kikwit.gtf")

# ----------------------
#   Extract Novel 
# ----------------------
# Keep only novel (u) and antisense (x)
mask_ux <- ribodepl_novel_compared_to_polya$class_code == "u" | ribodepl_novel_compared_to_polya$class_code == "x" ;
mask_ux[is.na(mask_ux)] <- FALSE
# comparing the predicted one with the reference.
novel_transcript_ids <- ribodepl_novel_compared_to_polya[mask_ux]$transcript_id
novel_ribodepl <- ribodepl_novel[ribodepl_novel$transcript_id %in% novel_transcript_ids]
length(unique(novel_ribodepl$transcript_id))
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



gtf_polya_newids <- add_prefix(polyA_novel, "polya-")
gtf_ribodepl_newids <- add_prefix(novel_ribodepl, "ribodepl-")


# Identify ribodepleted transcripts, whose gene name was already present in the poly(A) ones.

#gtf_polya_newids_names <- add_gene_name(gtf_polya_newids)
#gtf_ribodepl_newids_names <- add_gene_name(gtf_ribodepl_newids)

gtf_polya_newids_names <- gtf_polya_newids
gtf_ribodepl_newids_names <- gtf_ribodepl_newids
unique(gtf_ribodepl_newids_names$gene_name[gtf_ribodepl_newids_names$gene_name %in% gtf_polya_newids_names$gene_name])

# ---------------------------------------
# Merge them together ( just append ) 
# ---------------------------------------
glst <-list(gtf_polya_newids_names, gtf_ribodepl_newids_names)
gtf <- do.call(base::c, as(glst, "GRangesList"))

# ---------------------------------------
# Merge Novel and Ref  together ( just append ) 
# ---------------------------------------
glst_ref_and_novel <-list(ref_gtf,gtf)
gtf_ref_and_novel <- do.call(c, as(glst_ref_and_novel, "GRangesList"))

gtf_ref_and_novel_sort <- sort(gtf_ref_and_novel)
export(gtf_ref_and_novel_sort, outfile)  

