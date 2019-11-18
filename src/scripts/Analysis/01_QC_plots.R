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
projdir = "/home/luisas/Desktop/cluster/proj"
datadir <- "/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq"

scripts_dir = file.path(projdir,"code/ebola/src/scripts/Analysis/utils")
source(file.path(scripts_dir,"fastqc_results.heatmap.R"))
source(file.path(scripts_dir,"mapping_qa.R"))

out_dir = file.path(projdir,"data/plots")


# -----------------------
#         FASTQC
# -----------------------

inpath_fqc=file.path(baseDir,"data/02_RNA-Seq/02_fastqc/Zyagen")
outpath_fqc= paste0(out_dir,"02_fastqc/Zyagen")

plot_heatmap_fast_zyagen(inpath_fqc,outpath_fqc)

inpath_fqc_trimmomatic=file.path(baseDir,"data/02_RNA-Seq/02_fastqc_filtered/Zyagen")
outpath_fqc_trimmomatic= paste0(out_dir,"02b_fastqc/Zyagen")

plot_heatmap_fast_zyagen(inpath_fqc_trimmomatic,outpath_fqc_trimmomatic)

inpath_fqc_batch=file.path(baseDir,"data/02_RNA-Seq/02_fastqc/Batch01")
outpath_fqc_batch= paste0(out_dir,"02_fastqc/Batch01")

plot_heatmap_fast_zyagen(inpath_fqc_batch,outpath_fqc_batch)

inpath_fqc_batch=file.path(baseDir,"data/02_RNA-Seq/02_fastqc_filtered/Batch01")
outpath_fqc_batch= paste0(out_dir,"02b_fastqc/Batch01")
plot_heatmap_fast_zyagen(inpath_fqc_batch,outpath_fqc_batch)

# -----------------------
#         MAPPING
# -----------------------


inpath_hisat=file.path(baseDir,"data/02_RNA-Seq/03_hisat/Zyagen/")
outpath_hisat= paste0(out_dir,"03_hisat/Zyagen")

inpath_hisat_trimmed=file.path(baseDir,"data/02_RNA-Seq/03b_hisat/Zyagen/")
outpath_hisat_trimmed= paste0(out_dir,"03b_hisat/Zyagen")

inpath_hisat_batch=file.path(baseDir,"data/02_RNA-Seq/03_hisat/Batch01/")
outpath_hisat_batch= paste0(out_dir,"03_hisat/Batch01")

inpath_hisat_batch_trimmed=file.path(baseDir,"data/02_RNA-Seq/03b_hisat/Batch01/")
outpath_hisat_batch_trimmed= paste0(out_dir,"03b_hisat/Batch01")

# Mapped reads
barplot_mapped_reads(inpath_hisat,outpath_hisat)
barplot_mapped_reads(inpath_hisat_trimmed,outpath_hisat_trimmed)
barplot_mapped_reads(inpath_hisat_batch,outpath_hisat_batch)
barplot_mapped_reads(inpath_hisat_batch_trimmed,outpath_hisat_batch)


barplot_mapped_reads_two(inpath_hisat,inpath_hisat_batch,outpath_hisat)
barplot_counts_stats_2(inpath_hisat,inpath_hisat_batch,outpath_hisat)



# Infer Strandness
plot_infer_strandness(inpath_hisat,outpath_hisat)

# Barplot counts_stats
barplot_counts_stats(inpath_hisat,outpath_hisat)

# Barplot reads distribution  


pattern_distr = ".read_distribution.txt"
inpath_hisat=file.path(datadir,"02_RNA-Seq_old/03_hisat/Batch01")

distr_files_un <- distr_files[!str_detect(distr_files,"chrUn")]
outpath_hisat= file.path(out_dir,"03_hisat/Batch01")
plot_read_distribution(distr_files_un,outpath_hisat,"filtered_out_readdistribution.png", yli = 20000000)




distr_files <- iterate_files(inpath_hisat, pattern_distr)
distr_files <- distr_files[!str_detect(distr_files,"nvsf")]
distr_files <- distr_files[!str_detect(distr_files,"filter")]
distr_files_2 <- iterate_files(inpath_hisat_batch, pattern_distr)
distr_files_2 <- distr_files_2[!str_detect(distr_files_2,"nvsf")]
distr_files_2 <- distr_files_2[!str_detect(distr_files_2,"filter")]
plot_read_distribution_two(distr_files,distr_files_2,outpath_hisat,"filtered_out_readdistribution_2_tog.png", yli = 80000000)



distr_files_nvsf <-distr_files[str_detect(distr_files,"read_distribution_nvsf")]
plot_read_distribution(distr_files_nvsf,outpath_hisat,"filtered_out_readdistribution.png", yli = 20000000)






