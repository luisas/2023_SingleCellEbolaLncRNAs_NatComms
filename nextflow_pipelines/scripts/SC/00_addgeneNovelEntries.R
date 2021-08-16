library(rtracklayer)

# Gene annotation
#ref <- import(file.path("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf"))
ref <- import("/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf")

# I should add before the 
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

# Only check colocation of DE genes
gene_entries <- ref[ref$type == "gene", ]$gene_name
# ADD Missing entries
missing_gene_lines <- unique(ref$gene_name )[!(unique(ref$gene_name ) %in% gene_entries)]
gr_list <- lapply(missing_gene_lines, function(x) get_gene_range(x, ref))
missing_gene_ranges <- do.call(base::c,gr_list)
# Add missing gene entries back to reference
ref_genes <-c(ref, missing_gene_ranges)
export(ref_genes, "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames_addedgenesNovel.gtf")