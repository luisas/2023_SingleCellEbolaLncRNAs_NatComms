## QC Plotting 
# Libraries loading
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ComplexHeatmap)
library(dendextend)
library(RColorBrewer)
# Options
options(stringsAsFactors = F)



# Output Path 
scripts_dir = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/code/ebola/src/scripts/"
out_dir = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/plots/"


# -----------------------
#         FASTQC
# -----------------------
source(paste0(scripts_dir,"fastqc_results.heatmap.R"))

inpath_fqc="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/02_fastqc/Zyagen"
outpath_fqc= paste0(out_dir,"02_fastqc/Zyagen")

plot_heatmap_fast_zyagen(inpath_fqc,outpath_fqc)

inpath_fqc_trimmomatic="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/02_fastqc_filtered/Zyagen"
outpath_fqc_trimmomatic= paste0(out_dir,"02b_fastqc/Zyagen")

plot_heatmap_fast_zyagen(inpath_fqc_trimmomatic,outpath_fqc_trimmomatic)

inpath_fqc_batch="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/02_fastqc/Batch01"
outpath_fqc_batch= paste0(out_dir,"02_fastqc/Batch01")

plot_heatmap_fast_zyagen(inpath_fqc_batch,outpath_fqc_batch)

# -----------------------
#         MAPPING
# -----------------------
source(paste0(scripts_dir,"mapping_qa.R"))

inpath_hisat="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/03_hisat/Zyagen/"
outpath_hisat= paste0(out_dir,"03_hisat/Zyagen")
inpath_hisat_trimmed="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/03b_hisat/Zyagen/"
outpath_hisat_trimmed= paste0(out_dir,"03b_hisat/Zyagen")
inpath_hisat_batch="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/03_hisat/Batch01/"
outpath_hisat_batch= paste0(out_dir,"03_hisat/Batch01")

# Mapped reads
barplot_mapped_reads(inpath_hisat,outpath_hisat)
barplot_mapped_reads(inpath_hisat_trimmed,outpath_hisat_trimmed)

barplot_mapped_reads(inpath_hisat_batch,outpath_hisat_batch)

# Infer Strandness
plot_infer_strandness(inpath_hisat,outpath_hisat)

# Barplot counts_stats
barplot_counts_stats(inpath_hisat,outpath_hisat)

# Barplot reads distribution  
plot_read_distribution(inpath_hisat,outpath_hisat)

