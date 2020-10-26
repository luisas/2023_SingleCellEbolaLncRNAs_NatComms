ref <- import(file.path("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf"))
ref[substr(ref$gene_id, 1,3 ) == "MST",]$exon_id <- paste(ref[substr(ref$gene_id, 1,3 ) == "MST",]$transcript_id, ref[substr(ref$gene_id, 1,3 ) == "MST",]$exon_number, sep = ".")
export(ref,file.path("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames_exonid.gtf") )
