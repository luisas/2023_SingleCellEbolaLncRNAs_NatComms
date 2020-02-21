library(rtracklayer)
library(GenomicRanges)

datadir <- "/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all"
novel_concordant <- import(file.path(datadir,"03_novel_lncrnas/00_all_novels/novel_rhemac10_concordant_ribodepleted.gtf"))

dir_counts_ref = file.path(datadir, "04_Quantification_correct/")
htseq_files_ref <- iterate_files(dir_counts_ref, ".gene_counts.tab")
batch <- htseq_files_ref[str_detect(htseq_files_ref, "Batch01")]
zyagen <- htseq_files_ref[str_detect(htseq_files_ref, "Zyagen")]
htseq_files_ref <- c(batch,zyagen)

extract_expressed_only <- function(){
  dds_ref <- create_dds(htseq_files_ref,dir_counts_ref)
  dge_ref <- DGEList(counts = assays(dds_ref)$counts, group = dds_ref$dpo, genes = rownames(dds_ref), samples = colData(dds_ref))
  # Calc logcpm 
  assays(dds_ref)$logCPM <- edgeR::cpm(dge_ref, log = TRUE,prior.count = 0.5)
  # Create dataframe w/ logcpms
  expression_logcpm <-  as.matrix(assays(dds_ref)$logCPM)
  ## ADD only if we want to select a certain tissue
  #expression_logcpm <- expression_logcpm[,str_detect(colnames(expression_logcpm), "Liver")]
  # select only when > than 1 at least 3 tissues
  mask<- rowSums(as.data.frame(expression_logcpm) > 1) > 2
  expression_logcpm <- expression_logcpm[mask,]
  max_expression <- data.frame(id=rownames(expression_logcpm), expr=rowMax(expression_logcpm))
  
}
