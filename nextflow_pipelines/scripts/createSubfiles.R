#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(rtracklayer)
library(readr)

# Rhemac10 from UCSC
#gtf <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames_UCSC.gtf")
gtf <- import(args[1])
sourceTable <- read.table(args[2])
outdir <- args[3]
prefix <- args[4]
ref_novel_full<- import(args[5])
chrom<- read_csv(args[6])
# Get the biotype from the official table 
#sourceTable <- read.table("/home/luisas/Desktop/slncky-master/annotations/ensemblSource.txt")
#chrom<- read_tsv("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/slncky-master/annotations/chromosomes_in_fasta.tsv")
table(sourceTable$V2)
gtf_transcripts <- gtf[gtf$type == "transcript", ]

# Select Protein Coding and pseudogenes 
pc_ids <- sourceTable[sourceTable$V2 == "protein_coding",]
pseudogene_ids <- sourceTable[sourceTable$V2 %in% c("pseudogene", "processed_pseudogene"),]
lnc_ids <- sourceTable[sourceTable$V2 == "lncRNA",]
mi_ids <- sourceTable[sourceTable$V2 == "miRNA",]
sno_ids <- sourceTable[sourceTable$V2 == "snoRNA",]

pc_ids$V1 <- unlist(lapply(as.character(pc_ids$V1), function(x) strsplit(x, ".", fixed = T)[[1]][1]))
pseudogene_ids$V1 <- unlist(lapply(as.character(pseudogene_ids$V1), function(x) strsplit(x, ".", fixed = T)[[1]][1]))
lnc_ids$V1 <- unlist(lapply(as.character(lnc_ids$V1), function(x) strsplit(x, ".", fixed = T)[[1]][1]))
mi_ids$V1 <- unlist(lapply(as.character(mi_ids$V1), function(x) strsplit(x, ".", fixed = T)[[1]][1]))
sno_ids$V1 <- unlist(lapply(as.character(sno_ids$V1), function(x) strsplit(x, ".", fixed = T)[[1]][1]))


#outdir <- "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/slncky-master/annotations/"
#prefix <- "rheMac10"

print(head(pc_ids))
print(head(gtf$transcript_id))
print(gtf[gtf$transcript_id %in% pc_ids$V1, ])

# Prepare separate gtfs for running slnky 
export(gtf[gtf$transcript_id %in% pc_ids$V1, ], file.path(outdir, paste(prefix,".ensGene.PC.gtf", sep = ""))) 
export(gtf[gtf$transcript_id %in% pseudogene_ids$V1, ], file.path(outdir, paste(prefix,".ensGene.Pseudogenes.gtf", sep = "")))
export(gtf[gtf$transcript_id %in% lnc_ids$V1, ], file.path(outdir, paste(prefix,".ensGene.lnc.gtf", sep = ""))) 
export(gtf[gtf$transcript_id %in% mi_ids$V1, ], file.path(outdir, paste(prefix,".ensGene.mirna.gtf", sep = ""))) 
export(gtf[gtf$transcript_id %in% sno_ids$V1, ], file.path(outdir, paste(prefix,".ensGene.snorna.gtf", sep = ""))) 

# Same chrnaming - remove not coincidin contigs
#ref_novel_full<- import(file.path(outdir, "rheMac10_EBOV_and_novel_genenames_UCSC.gtf"))
ref_novel_full <- ref_novel_full[seqnames(ref_novel_full) %in% unlist(chrom),]
#export(ref_novel_full[ref_novel_full$type == "exon", ], "/home/luisas/Desktop/slncky-master/rheMac10_EBOV_and_novel_genenames_UCSC.exons.gtf")
export(ref_novel_full, file.path(outdir, "rheMac10_EBOV_and_novel_genenames_UCSC.contigsfiltered.gtf"))
#export(ref_novel_full[ref_novel_full$type == "gene", ], "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/slncky-master/rheMac10_EBOV_and_novel_genenames_UCSC.gene.gtf")



