shhh <- suppressPackageStartupMessages
shhh(library(rtracklayer))

args = commandArgs(trailingOnly=TRUE)
file = args[1]
#file = "/home/luisas/Desktop/cluster/gene_annotations/gencode/release_26/gencode.v26.GRCh38.annotation.gtf"
output = args[2]
filter = args[3]

gtf<- import(file)

# If specified in the input, only select non coding rnas from set
if( filter == "lncrna"){
  print("!!!!!!!!!!!!!!! Filtering only lncrnas")
  gtf[substr(gtf$gene_id,1,3) == "pol",]$gene_biotype <- "lncRNA"
  gtf[substr(gtf$gene_id,1,3) == "rib",]$gene_biotype <- "lncRNA"
  gtf <- gtf[!is.na(gtf$gene_biotype) &gtf$gene_biotype == "lncRNA"]
}

if( filter == "mrna"){
  print("!!!!!!!!!!!!!!! Filtering only proteincoding")
  gtf <- gtf[!is.na(gtf$gene_biotype) &gtf$gene_biotype == "protein_coding"]
}



# Select only exons
#gtf_exons <- gtf[gtf$type == "exon",]
gtf_exons <- gtf

# Add an id containing gene and exon id
#gtf_exons$complete_id <- paste(gtf_exons$gene_id, gtf_exons$exon_id, sep  = "_")
gtf_exons$complete_id <- gtf_exons$gene_id

# Change name to chromosomes back to ucsc annotation
seqlevelsStyle(gtf_exons) <- "UCSC"

# Create bed file format (n = 4)
df <- NULL
df$seq <-seqnames(gtf_exons)
df$start <- start(ranges(gtf_exons))
df$end <- end(ranges(gtf_exons))
df$name <- gtf_exons$complete_id
df$score  <- rep("-", length(df$name ))
df$strand <- strand(gtf_exons)



if(filter == "addbiotype"){
  print("Adding gene biotype")
  df$biotype <- gtf_exons$gene_type
}


df <- as.data.frame(df)

# Save the file
write.table(df, output,sep="\t",row.names=FALSE, col.names = FALSE, quote = FALSE)
