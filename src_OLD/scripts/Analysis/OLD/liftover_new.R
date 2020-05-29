library(rtracklayer)
library(IRanges)
library(GenomicRanges)
library(ggplot2);

datadir <- "/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq"
dir_assemblies<- "/home/luisas/Desktop/cluster/assemblies"


human_annotation <- import(file.path(dir_assemblies, "ensembl/release-96/Homo_sapiens_hg38/Homo_sapiens.GRCh38.96.gtf"))
seqlevelsStyle(human_annotation) <- "UCSC"

#Run Liftover - outputs a granges list 
Human_Orthologs <- import("/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq/02_RNA-Seq_rheMac10/07_liftover/new.gff")



# Extract overlaps from corresponding human annotation 
overlaps <- lapply( Human_Orthologs, function(x) findOverlaps(x,human_annotation))

# See what is in the overlapping regions
hits_in_human <-lapply(overlaps, function(x) human_annotation[subjectHits(x)])

# Obtain the classes 
orthologs_class <- (unlist(lapply(hits_in_human, function(x) (unique(x$gene_biotype)))))


p1 <- ggplot(data.frame(orthologs_class), aes(x=orthologs_class, fill = orthologs_class)) +
  geom_bar()+
  coord_flip()

plot_path <- file.path(datadir,"/02_RNA-Seq_rheMac10/07_liftover/plot_orthologs.png")
ggsave(plot_path, plot = p1,height=30,width=30)




