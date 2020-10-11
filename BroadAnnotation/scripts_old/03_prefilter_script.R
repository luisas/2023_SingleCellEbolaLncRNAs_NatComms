
library(rtracklayer)
library(dplyr, warn.conflicts = FALSE)

exon_nr_summary <- function(gr){
  get_nr_exons_mod<- function(gr){
    df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id,"exon_number" = as.numeric(gr$exon_number))
    number_exons <- df %>% dplyr::group_by(gene_id,transcript_id) %>%dplyr::summarize(max_exon = max(exon_number))
    return(number_exons)
  }
  gr <- gr[gr$type =="exon",]
  #gr <- get_only_max_transcript(gr)
  df_l <- data.frame(get_nr_exons_mod(gr))
  return(df_l)
}

calc_transcript_length <- function(gr){
  gr <- gr[gr$type =="exon",]
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "range_width" = width(ranges(gr)))
  collapsed <- df%>% dplyr::group_by(gene_id,transcript_id) %>% dplyr::summarize("range" = sum(range_width))
  return(collapsed)
}

bed12_entry <- function(transcript){
  
  entry <- gtf_no_gene[gtf_no_gene$transcript_id == transcript, ]
  exons <- entry[entry$type == "exon",]
  
  # Extract all informations 
  chrom <- as.character(seqnames(entry))[1]
  start <- min(start(ranges(exons)))-1
  stop <- max(end(ranges(exons)))
  name <- exons$transcript_id[1]
  score <- "-"
  strand <- as.character(strand(exons))[1]
  thickstart<- start
  thickend <- stop
  rgb <- 0 
  blockcount <- length(exons) 
  blocksizes <- paste(width(ranges(exons)), collapse = ",")
  blockstarts <- paste(start(ranges(exons))-1-start , collapse= ",")
  line <- data.frame(chrom, start, stop, name, score, strand, thickstart, thickend, rgb, blockcount, blocksizes, blockstarts)
  
  return(line)
}


# Read command line
args = commandArgs(trailingOnly=TRUE)

lnc_novel_compared_to_ref_file <- args[1]
assembly <- import(args[2])
ref <- import(args[3])
outfile <- args[4]

# ---------------------------------
#        ABOUT THE ASSEMBLY
# ---------------------------------

# How many transcripts do we assemble ? 
#assembly <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_ribodepl/04_stringtie/stringtie_merged_reference_guided.gtf")
#lnc_novel_compared_to_ref_file <- ("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_ribodepl/04_stringtie_gffcompare/01_gffCompare/merged.annotated.gtf")

# ------------------------------- INSPECT ----------------------------------
assembly_stringtie <- assembly[assembly$source == "StringTie",]
assembly_ensembl <- assembly[assembly$source == "ensembl",]

# How many transcripts do we assemble in total ?
length(unique(assembly_stringtie$transcript_id))
length(unique(assembly_stringtie$gene_id))

# Assess the number of exons - how many monoexonic 
summary_exons <- exon_nr_summary(assembly_stringtie)
monoexonic_transcript_ids <- summary_exons[summary_exons$max_exon == 1,]$transcript_id
#plot_assembly_stats(summary_exons)
# ----------------------------------------------------------------------------

print("After inspection ")

# ------------------------------- FILTER 1: retain only novel ones  --------------

# Get only novel with GFF compare 
lnc_novel_compared_to_ref <- import(lnc_novel_compared_to_ref_file)

# Keep only novel (u) and antisense (x)
mask_ux <- lnc_novel_compared_to_ref$class_code == "u" | lnc_novel_compared_to_ref$class_code == "x" ;
mask_ux[is.na(mask_ux)] <- FALSE

# comparing the predicted one with the reference.
novel_transcript_ids <- lnc_novel_compared_to_ref[mask_ux]$transcript_id
assembly_stringtie <- assembly_stringtie[assembly_stringtie$transcript_id %in% novel_transcript_ids]
print("After filter 0:")
length(unique(assembly_stringtie$gene_id))
length(unique(assembly_stringtie$transcript_id))
# ------------------------------------------------------------------------------------------------------------------

print("After filter1 ")



# ------------------------------- FILTER 1: remove monoexonic transcripts ----------------------------------
# Remove monoexonic entries 
assembly_stringtie_multiple_exons <- assembly_stringtie[!(assembly_stringtie$transcript_id %in% monoexonic_transcript_ids),]
print("After filter 1:")
length(unique(assembly_stringtie_multiple_exons$gene_id))
length(unique(assembly_stringtie_multiple_exons$transcript_id))
# ------------------------------------------------------------------------------------------------------------------

print("After filter1b ")

#export(assembly_stringtie_multiple_exons,"/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_ribodepl/04_stringtie/stringtie_merged_reference_guided_multipleexons.gtf")


# ------------------------------- FILTER 2: Remove transcripts shorter than 200 nucleotides --------------
transcript_length <- calc_transcript_length(assembly_stringtie_multiple_exons)
short_transcripts_ids <- transcript_length[transcript_length$range < 200,]$transcript_id
assembly_filtered <- assembly_stringtie_multiple_exons[!(assembly_stringtie_multiple_exons$transcript_id %in% short_transcripts_ids),]
print("After filter 2:")
length(unique(assembly_filtered$gene_id))
length(unique(assembly_filtered$transcript_id))
# ------------------------------------------------------------------------------------------------------------------
print("After filter2 ")



# ------------------------------- FILTER 3: Remove anything overlapping protein coding genes --------------
# Remove anything overlapping a protein coding genes 
ref_mrna <- ref[!is.na(ref$gene_biotype) & ref$gene_biotype == "protein_coding",]
ref_mrna_exons <- ref_mrna[ref_mrna$type == "exon",]
overlap <- findOverlaps(assembly_filtered, ref_mrna )
remove <- assembly_filtered[unique(queryHits(overlap))]$transcript_id
length(unique(remove))
assembly_filtered <- assembly_filtered[!(assembly_filtered$transcript_id %in% remove)]
print("After filter 3:")
length(unique(assembly_filtered$gene_id))
length(unique(assembly_filtered$transcript_id))
# ------------------------------------------------------------------------------------------------------------------
print("After filter3 ")


#export(assembly_filtered, "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list_CPC2/candidates.gtf" ) 
candidates <- assembly_filtered

export(candidates,outfile)


