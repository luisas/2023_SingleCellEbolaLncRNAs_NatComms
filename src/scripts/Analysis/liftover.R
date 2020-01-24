library(rtracklayer)
library(IRanges)
library(GenomicRanges)
library(ggplot2);

datadir <- "/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq"
dir_assemblies<- "/home/luisas/Desktop/cluster/assemblies"
human_annotation <- import(file.path(dir_assemblies, "ensembl/release-96/Homo_sapiens_hg38/Homo_sapiens.GRCh38.96.gtf"))
seqlevelsStyle(human_annotation) <- "UCSC"

#Chain File
chainfile <- file.path(datadir,"/01_PreliminaryFiles_rheMac10/gene_annotations/chains/rheMac10ToHg38.over.chain")
chain <- import.chain(chainfile) # Macaque!

#Novel lncrnas in macaque
macaque <- import(file.path(datadir,"/02_RNA-Seq_rheMac10/07_liftover/novel_rhemac10.gtf"))
macaque <- import(file.path(datadir,"/02_RNA-Seq_rheMac10/07_liftover/novel_rhemac10.bed"))
head(macaque)
seqlevelsStyle(macaque) <- "UCSC"
head(macaque)
export(macaque, file.path(datadir,"/02_RNA-Seq_rheMac10/07_liftover/novel_rhemac10_ucsc.gtf"))
head(macaque)
#Run Liftover - outputs a granges list 
Human_Orthologs <- liftOver(macaque, chain)

# Quick stats: For how many we did not find the orthologs
MISSED <- length(unique(macaque$gene_id)) - length(unique(unlist(Human_Orthologs)$gene_id))
FOUND <-  length(unique(unlist(Human_Orthologs)$gene_id))
print(paste0("#Genes where no ortholog was identified: ",MISSED))
print(paste0("#Genes where ortholog was identified: ",FOUND))

export()
# Only select the whole region 
gr_toHuman_collapsed <- lapply(Human_Orthologs, function(x) {range(split(x, ~exon_number))})

# Extract overlaps from corresponding human annotation 
overlaps <- lapply( gr_toHuman_collapsed, function(x) findOverlaps(x,human_annotation))

# See what is in the overlapping regions
hits_in_human <-lapply(overlaps, function(x) human_annotation[subjectHits(x)])

# Obtain the classes 
orthologs_class <- (unlist(lapply(hits_in_human, function(x) (unique(x$gene_biotype)))))


p1 <- ggplot(data.frame(orthologs_class), aes(x=orthologs_class, fill = orthologs_class)) +
  geom_bar()+
  coord_flip()

plot_path <- file.path(datadir,"/02_RNA-Seq_rheMac10/07_liftover/plot_orthologs.png")
ggsave(plot_path, plot = p1,height=30,width=30)




