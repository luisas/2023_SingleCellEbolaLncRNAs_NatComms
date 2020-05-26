
library(rtracklayer)
library(dplyr)

# Read command line
args = commandArgs(trailingOnly=TRUE)
ribodepl_novel_compared_to_polya <- import(args[1])
polyA_novel <- import(args[2])
ref_gtf <- import(args[3])
outfile <- args[4]


# ----------------------
#   Extract Novel 
# ----------------------
# Keep only novel (u) and antisense (x)
mask_ux <- ribodepl_novel_compared_to_polya$class_code == "u" | ribodepl_novel_compared_to_polya$class_code == "x" ;
mask_ux[is.na(mask_ux)] <- FALSE
# comparing the predicted one with the reference.
novel_transcript_ids <- ribodepl_novel_compared_to_polya[mask_ux]$transcript_id
novel_ribodepl <- ribodepl_novel_compared_to_polya[ribodepl_novel_compared_to_polya$transcript_id %in% novel_transcript_ids]

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
gtf_polya_newids_names <- add_gene_name(gtf_polya_newids)
gtf_ribodepl_newids_names <- add_gene_name(gtf_ribodepl_newids)

unique(gtf_ribodepl_newids_names$gene_name[gtf_ribodepl_newids_names$gene_name %in% gtf_polya_newids_names$gene_name])

# ---------------------------------------
# Merge them together ( just append ) 
# ---------------------------------------
glst <-list(gtf_polya_newids_names, gtf_ribodepl_newids_names)
gtf <- do.call(c, as(glst, "GRangesList"))

# ---------------------------------------
# Merge Novel and Ref  together ( just append ) 
# ---------------------------------------
glst_ref_and_novel <-list(ref_gtf,gtf)
gtf_ref_and_novel <- do.call(c, as(glst_ref_and_novel, "GRangesList"))

gtf_ref_and_novel[is.na(gtf_ref_and_novel$gene_name)]$gene_name <- paste(gtf_ref_and_novel[is.na(gtf_ref_and_novel$gene_name)]$gene_id, "-unkown", sep = "") 
gtf_ref_and_novel[is.na(gtf_ref_and_novel$transcript_name)]$transcript_name <-paste(gtf_ref_and_novel[is.na(gtf_ref_and_novel$transcript_name)]$transcript_id, "-unkown", sep = "")  

export(gtf_ref_and_novel, outfile)  

