
library(rtracklayer)
library(dplyr)
library(tximport)
library(readr)
library(matrixStats)
library(forcats)


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
gtf <- args[1]
output <- args[2]

# Import file 
candidates<- import(gtf)
gtf_no_gene <- candidates[candidates$type != "gene", ]

# Extract all transcripts
transcript_ids <- (unique(gtf_no_gene[!is.na(gtf_no_gene$transcript_id), ]$transcript_id))
transcript_ids_novel <- (unique(gtf_no_gene[substr(gtf_no_gene$transcript_id,1,4) == "MSTR",]$transcript_id))
length(unique(transcript_ids_novel))
dfset <- lapply(transcript_ids_novel, function(x) (bed12_entry(x)))
df <- do.call(rbind, dfset)

# Save
write.table(df, output,sep="\t",row.names=FALSE, col.names = FALSE, quote = FALSE)


