
library(rtracklayer)
library(dplyr)

# Read command line
args = commandArgs(trailingOnly=TRUE)

lnc_novel_compared_to_ref_file <- args[1]
mrnas_predicted_file <- args[2]
outfile <- args[3]

# ----------------------
#   Extract Novel 
# ----------------------
lnc_novel_compared_to_ref <- import(lnc_novel_compared_to_ref_file)
# Keep only novel (u) and antisense (x)
mask_ux <- lnc_novel_compared_to_ref$class_code == "u" | lnc_novel_compared_to_ref$class_code == "x" ;
mask_ux[is.na(mask_ux)] <- FALSE
# comparing the predicted one with the reference.
novel_transcript_ids <- lnc_novel_compared_to_ref[mask_ux]$transcript_id
novel <- lnc_novel_compared_to_ref[lnc_novel_compared_to_ref$transcript_id %in% novel_transcript_ids]

# ---------------------------
#   Extract non-concordant
# ---------------------------
mRNAs_out <- import(mrnas_predicted_file)
df <- data.frame()
df <- rbind(df, data.frame(gene_id = novel$gene_id, transcript_id = novel$transcript_id, pred = "lncrna" ))
df <- rbind(df, data.frame(gene_id = mRNAs_out$gene_id, transcript_id = mRNAs_out$transcript_id, pred = "mrna" ))
unconcordant_prediction <- df %>%  dplyr::group_by(gene_id) %>% dplyr::summarise(Unique_Elements =  dplyr::n_distinct(pred)) %>%  dplyr::filter( Unique_Elements > 1)
novel_concordant <- novel[!(novel$gene_id %in% unconcordant_prediction$gene_id)]

export(novel_concordant,outfile)


