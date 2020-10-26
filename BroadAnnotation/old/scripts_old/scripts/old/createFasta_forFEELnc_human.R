knitr::opts_chunk$set(echo = TRUE)
library(seqinr)
library(rtracklayer)
library(Biostrings)
library(BSgenome)

extract_transcript_sequences <- function(transcript_id, fasta, gtf_exons, type, output_file){
  
  transcript_exons <- gtf_exons[gtf_exons$transcript_id == transcript_id, ]
  gene_id <- transcript_exons$gene_id[1]
  strand <- as.vector(strand(transcript_exons))[1]
  chr <- as.vector(seqnames(transcript_exons)[1])
  mygrs <- GRanges(seqnames = chr,
                   ranges=IRanges::reduce(ranges(transcript_exons)),
                   strand=strand) 
  # Create ID 
  pos <- paste(chr, paste(start(mygrs)[1],rev(end(mygrs))[1], sep ="-"), sep = ":")
  
  
  id <- paste(transcript_id, gene_id,strand,pos, type, sep = "_")
  seq <- paste(as.character(getSeq(fasta, mygrs)), collapse = "")
  write.fasta(seq, id, nbchar = 60,output_file, open = "a")
  
}
prefix <- "/gpfs/projects/bsc83/Data/"
output_file_lncrna_human <-  file.path(prefix,'Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/Homo_sapiens.GRCh38.100_known_lncrna.fa')
output_file_mrna_human <-  file.path(prefix, 'Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/Homo_sapiens.GRCh38.100_known_proteincoding.fa')

file.create(output_file_lncrna_human)
file.create(output_file_mrna_human)
gtf_human <- import(file.path(prefix, "gene_annotation/ensembl_release100/homo_sapiens/Homo_sapiens.GRCh38.100.gtf"))
fasta_human <- readDNAStringSet(file.path(prefix, "assemblies/ensembl/release-100/homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa"))
saveRDS(fasta_human, file.path(prefix,"Ebola/RObjects/fasta_human.rds"))
print("File Read")
fasta_human <- readRDS(file.path(prefix,"Ebola/RObjects/fasta_human.rds"))

gtf_human


names(fasta_human) <- unlist(lapply(names(fasta_human), function(x) strsplit(x, " ")[[1]][1]))

# --------------------------
#   Prepare files Human 
# --------------------------
known_lnc_human <- gtf_human[gtf_human$gene_biotype =="lncRNA", ]
known_lnc_human <- known_lnc_human[known_lnc_human$type == "transcript", ]

known_pc_human <- gtf_human[gtf_human$gene_biotype =="protein_coding", ]
known_pc_human <- known_pc_human[known_pc_human$type == "transcript", ]

gtf_exons_human <- gtf_human[gtf_human$type %in% c("exon", "CDS", "start_codon", "stop_codon"),]

length(unique(known_lnc_human$transcript_id))

lapply(unique(known_pc_human$transcript_id), function(x) extract_transcript_sequences(x, fasta_human, gtf_exons_human, "protein_coding", output_file_mrna_human))
lapply(unique(known_lnc_human$transcript_id), function(x) extract_transcript_sequences(x,fasta_human,gtf_exons_human, "lncRNA", output_file_lncrna_human))
