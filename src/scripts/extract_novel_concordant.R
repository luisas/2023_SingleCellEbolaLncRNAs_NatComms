library(rtracklayer)
library(GenomicRanges)

datadir <- "/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all"

merged_annotated_ribodepl <- import(file.path(datadir, "02_RNA-Seq_BatchZyagen/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc/01_gffcompare/merged.annotated.gtf"))
mRNAs_out_ribodepl <- import(file.path(datadir, "02_RNA-Seq_BatchZyagen/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc/feelnc_codpot_out/candidate_lncRNA.gtf.mRNA.gtf"))

mRNAs_out_polya <- import(file.path(datadir, "02_RNA-Seq_external/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc/feelnc_codpot_out/candidate_lncRNA.gtf.mRNA.gtf"))
merged_annotated_polya <- import(file.path(datadir, "02_RNA-Seq_external/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc/01_gffcompare/merged.annotated.gtf"))

rhemac_ebov_ref <- import(file.path(datadir, "01_PreliminaryFiles_rheMac10/gene_annotations/rheMac10_EBOV-Kikwit.gtf"))

extract_novel_concordant <- function(merged_annotated, mRNAs_out){
  # Selects only u and x 
  mask_uix <- merged_annotated$class_code == "u" | merged_annotated$class_code == "x" ;  mask_uix[is.na(mask_uix)] <- FALSE
  # comparing the merged one with the reference. unkown, intergenic and intronic transcipts are kept.
  novel_transcript_ids <- merged_annotated[mask_uix]$transcript_id
  novel <- merged_annotated[merged_annotated$transcript_id %in% novel_transcript_ids]
  # -----------------
  # Remove Non-Concordant
  df <- data.frame()
  df <- rbind(df, data.frame(gene_id = novel$gene_id, transcript_id = novel$transcript_id, pred = "lncrna" ))
  df <- rbind(df, data.frame(gene_id = mRNAs_out$gene_id, transcript_id = mRNAs_out$transcript_id, pred = "mrna" ))
  unconcordant_prediction <- df %>%  dplyr::group_by(gene_id) %>% dplyr::summarise(Unique_Elements =  dplyr::n_distinct(pred)) %>%  dplyr::filter( Unique_Elements > 1)
  novel_concordant <- novel[!(novel$gene_id %in% unconcordant_prediction$gene_id)]
  return(novel_concordant)
}

# Extract the novel concordant for boths and create a new dataset
novel_concordant_ribodepl <- extract_novel_concordant(merged_annotated_ribodepl,mRNAs_out_ribodepl)
novel_concordant_polya <- extract_novel_concordant(merged_annotated_polya,mRNAs_out_polya)

# Merge them with the reference
glst <-list(rhemac_ebov_ref,novel_concordant_ribodepl, novel_concordant_polya)
gtf_ref_all_novel <- do.call(c, as(glst, "GRangesList"))

export(novel_concordant_ribodepl, (file.path(datadir,"03_novel_lncrnas/00_all_novels/novel_rhemac10_concordant_ribodepleted.gtf")))
export(novel_concordant_polya, (file.path(datadir,"03_novel_lncrnas/00_all_novels/novel_rhemac10_concordant_polyA.gtf")))
export(gtf_ref_all_novel, (file.path(datadir,"03_novel_lncrnas/00_all_novels/rheMac10_EBOV-Kikwit_and_both_novel.gtf")))


