
library(rtracklayer)
library(dplyr)

# Read command line
args = commandArgs(trailingOnly=TRUE)

lnc_novel_compared_to_ref_file <- args[1]
lnc_pred_file <- import(args[2])
mrnas_predicted_file <- args[3]
ref <- import(args[4])
outfile <- args[5]


#lnc_novel_compared_to_ref_file <- "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_external/05_feelNC_prediction/feelnc_gencode_linc/01_gffcompare/merged.annotated.gtf"
#lnc_pred_file <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_external/05_feelNC_prediction/feelnc_gencode_linc/feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf")
#mrnas_predicted_file <- "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_external/05_feelNC_prediction/feelnc_gencode_linc/feelnc_codpot_out/candidate_lncRNA.gtf.mRNA.gtf"
#ref <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/gene_annotations/rheMac10_EBOV-Kikwit.gtf")
#outfile <- args[5]

ref_mrna <- ref[!is.na(ref$gene_biotype) & ref$gene_biotype == "protein_coding",]
ref_mrna_exons <- ref_mrna[ref_mrna$type == "exon",]

# ----------------------
#   Extract Novel 
# ----------------------
lnc_novel_compared_to_ref <- import(lnc_novel_compared_to_ref_file)
# Keep only novel (u) and antisense (x)
mask_ux <- lnc_novel_compared_to_ref$class_code == "u" | lnc_novel_compared_to_ref$class_code == "x" ;
mask_ux[is.na(mask_ux)] <- FALSE
# comparing the predicted one with the reference.
novel_transcript_ids <- lnc_novel_compared_to_ref[mask_ux]$transcript_id
novel <- lnc_pred_file[lnc_pred_file$transcript_id %in% novel_transcript_ids]

overlap <- findOverlaps(novel, ref_mrna )
remove <- novel[unique(queryHits(overlap))]$transcript_id
novel <- lnc_novel_compared_to_ref[!(lnc_novel_compared_to_ref$transcript_id %in% remove)]

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


