shhh <- suppressPackageStartupMessages
shhh(library(rtracklayer))

args = commandArgs(trailingOnly=TRUE)
file = args[1]
#file = "/home/luisas/Desktop/cluster/gene_annotations/gencode/release_26/gencode.v26.GRCh38.annotation.gtf"
output = args[2]

file <- "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf"
gtf<- import(file)




# Remove gene entries ( not needed )
gtf_no_gene <- gtf[gtf$type != "gene", ]

# Extract all transcripts
transcript_ids <- (unique(gtf_no_gene[!is.na(gtf_no_gene$transcript_id), ]$transcript_id))

bed12_entry <- function(transcript){
  
  entry <- gtf_no_gene[gtf_no_gene$transcript_id == transcript, ]
  transcript_gr <- entry[entry$type == "transcript",]
  exons <- entry[entry$type == "exon",]
  
  # Extract all informations 
  chrom <- as.character(seqnames(entry))[1]
  start <- start(ranges(transcript_gr))-1
  stop <- end(ranges(transcript_gr))
  name <- transcript_gr$transcript_id
  score <- "-"
  strand <- as.character(strand(transcript_gr))[1]
  thickstart<- start
  thickend <- stop
  rgb <- 0 
  blockcount <- length(exons) 
  blocksizes <- paste(width(ranges(exons)), collapse = ",")
  blockstarts <- paste(start(ranges(exons))-1 , collapse= ",")
  
  line <- paste(chrom, start, stop, name, score, strand, thickstart, thickend, rgb, blockcount, blocksizes, blockstarts, collapse = "\n")
  
  return(line)
}


bed12_entry(transcript)


# Save the file
write.table(df, output,sep="\t",row.names=FALSE, col.names = FALSE, quote = FALSE)
