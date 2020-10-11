shhh <- suppressPackageStartupMessages
shhh(library(rtracklayer))

args = commandArgs(trailingOnly=TRUE)
file = args[1]
#file = "/home/luisas/Desktop/cluster/gene_annotations/gencode/release_26/gencode.v26.GRCh38.annotation.gtf"
output = args[2]

file <- "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel.gtf"
gtf<- import(file)
output <- "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.bed12"


file <- ("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_ribodepl/05_feelNC_prediction/feelnc_gencode_linc/candidate_lncRNA.gtf")
gtf<-import(file)


basename <- str_split(rev(str_split(file,"/")[[1]])[1], "\\.")[[1]][1]
output <- "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_ribodepl/05_feelNC_prediction/feelnc_gencode_linc/candidate_lncRNA.bed12"

# Remove gene entries ( not needed )
gtf_no_gene <- gtf[gtf$type != "gene", ]

# Extract all transcripts
transcript_ids <- (unique(gtf_no_gene[!is.na(gtf_no_gene$transcript_id), ]$transcript_id))

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

transcript_ids_novel <- (unique(gtf_no_gene[substr(gtf_no_gene$transcript_id,1,4) == "MSTR",]$transcript_id))
length(unique(transcript_ids_novel))

dfset <- lapply(transcript_ids_novel, function(x) (bed12_entry(x)))
df <- do.call(rbind, dfset)
# Save the file
write.table(df, output,sep="\t",row.names=FALSE, col.names = FALSE, quote = FALSE)
