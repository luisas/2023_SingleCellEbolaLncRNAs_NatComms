# ---- Testing
library(DESeq2)
library(rtracklayer)
library(edgeR)
library(zeallot)

# Load Utils
scripts_dir = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/code/ebola/src/scripts/"
source(paste0(scripts_dir,"de_utils.R"))

# GTF read
gtf = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/01_PreliminaryFiles/gene_annotations/rheMac8.ensembl_release97.gtf"
annotation <- readGFF(gtf)
annotation$size <- annotation$end -annotation$start

#Extract only genes from the annotation - this is as well what htseq does 
annotation_genes <- annotation[annotation$type == "gene",]
gene_lengths<- data.frame(annotation_genes$gene_id, annotation_genes$size)
colnames(gene_lengths) <- c('gene_id','basepairs')

# Temp 
directory <- "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/03_hisat"
# UMI_DEDUP 
sampleFiles <- iterate_files(directory, "UMI.f3.q60.umi_dedup.HTseq.gene_counts.tab")
sampleFiles_md <- iterate_files(directory, "UMI.f3.q60.md.HTseq.gene_counts.tab")

# --------------------------------------
dds <- create_dds(sampleFiles)
#dds <- add_gene_sizes(dds,gene_lengths)
assays(dds)$fpkm <- fpkm(dds)
# --------------------------------------
dds_md <- create_dds(sampleFiles_md)
dds_md <- add_gene_sizes(dds_md,gene_lengths)
assays(dds_md)$fpkm <- fpkm(dds_md)
# --------------------------------------
dds_filter <- create_dds(sampleFiles_filter)
dds_filter <- add_gene_sizes(dds_md,gene_lengths)
assays(dds_filter)$fpkm <- fpkm(dds_filter)
# --------------------------------------

umi_total_count <- sum(assays(dds)$counts)
md_total_count <- sum(assays(dds_md)$counts)
filter_total_count <- sum(assays(dds_filter)$counts)
dev.off()
getOption("device")
umi_total_count
md_total_count
filter_total_count
barplot(as.matrix(infer_stat_df),yaxt="n", ,las=2,yaxt='n', main="Inferred strandness")


barplot(c(umi_total_count,md_total_count,filter_total_count),las=2, col=brewer.pal(3,"Paired"))
# Plotting
out_dir = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/plots/03_hisat/Zyagen/"
plot_fpkm_heatmap(dds, paste(out_dir,"heatmap_umi_dedup_pearson.png"), "pearson")
plot_fpkm_heatmap(dds, paste(out_dir,"heatmap_umi_dedup_spearman.png"), "spearman")

# ----------------# ----------------# ----------------# ----------------# ----------------
# Marked Duplicates
# Plotting
out_dir = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/plots/03_hisat/Zyagen/"
plot_fpkm_heatmap(dds_md, paste(out_dir,"heatmap_md_pearson.png"), "pearson")
plot_fpkm_heatmap(dds_md, paste(out_dir,"heatmap_md_spearman.png"), "spearman")
# ---------------# ---------------# ---------------# ---------------# ---------------
# Marked Duplicates
sampleFiles_filter <- iterate_files(directory, "UMI.f3.q60.HTseq.gene_counts.tab")
# Plotting
out_dir = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/plots/03_hisat/Zyagen/"
plot_fpkm_heatmap(dds_filter, paste(out_dir,"heatmap_filter_pearson.png"), "pearson")
plot_fpkm_heatmap(dds_filter, paste(out_dir,"heatmap_filter_spearman.png"), "spearman")


