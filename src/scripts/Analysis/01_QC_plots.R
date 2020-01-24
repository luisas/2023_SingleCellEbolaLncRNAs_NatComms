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

dataset <- "02_RNA-Seq_external"
subset <- "SRP016501_rhesus"

#dataset <- "02_RNA-Seq_rheMac10"
#subset <- "Zyagen"

baseDir = "/home/luisas/Desktop/cluster/"
projdir = "/home/luisas/Desktop/cluster/proj"
datadir <- file.path("/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq/", dataset)

scripts_dir = file.path(projdir,"code/ebola/src/scripts/Analysis/utils")
source(file.path(scripts_dir,"fastqc_results.heatmap.R"))
source(file.path(scripts_dir,"mapping_qa.R"))

out_dir = file.path(projdir,"data/plots", dataset)
dir.create(out_dir)

# -----------------------
#         FASTQC
# -----------------------

inpath_fqc=file.path(datadir,"02_fastqc", subset)
outpath_fqc= file.path(out_dir,"/02_fastqc/", subset)
# Create dirs
dir.create(file.path(out_dir,"/02_fastqc/"))
dir.create(outpath_fqc)


#inpath_fqc_batch=file.path(datadir,"02_fastqc/Batch01")
#outpath_fqc_batch= paste0(out_dir,"02_fastqc/Batch01")

plot_heatmap_fast_zyagen(inpath_fqc,outpath_fqc, ".")
#plot_heatmap_fast_zyagen(inpath_fqc_batch,outpath_fqc_batch)


# -----------------------
#         MAPPING
# -----------------------


inpath_hisat=file.path(datadir,"03_hisat", subset)
outpath_hisat= file.path(out_dir,"/03_hisat/", subset)
# Create dirs
dir.create(file.path(out_dir,"/03_hisat/"))
dir.create(outpath_hisat)


#inpath_hisat_batch=file.path(datadir,"data/02_RNA-Seq/03_hisat/Batch01/")
#outpath_hisat_batch= paste0(out_dir,"03_hisat/Batch01")

# Mapped reads
source(file.path(scripts_dir,"mapping_qa.R"))
barplot_mapped_reads(inpath_hisat,outpath_hisat, lanes = FALSE)
barplot_mapped_reads(inpath_hisat_batch,outpath_hisat_batch)


#barplot_mapped_reads_two(inpath_hisat,inpath_hisat_batch,outpath_hisat)
#barplot_counts_stats_2(inpath_hisat,inpath_hisat_batch,outpath_hisat)

# Infer Strandness
#plot_infer_strandness(inpath_hisat,outpath_hisat)

# Barplot counts_stats
#barplot_counts_stats(inpath_hisat,outpath_hisat)
# Barplot reads distribution  


# -----------------------------
#         Read distribution 
# -----------------------------
pattern_distr = ".read_distribution.txt"
inpath_hisat=file.path(datadir,"/03_hisat/", subset)
distr_files <- iterate_files(inpath_hisat, pattern_distr)
distr_files_un <- distr_files[!str_detect(distr_files,"chrUn")]
outpath_hisat= file.path(out_dir,"03_hisat/", subset)

distr_files <- iterate_files(inpath_hisat, pattern_distr)
distr_files <- distr_files[!str_detect(distr_files,"nvsf")]
distr_files <-distr_files[str_detect(distr_files,"new")]
distr_files <- distr_files[!str_detect(distr_files,"filter")]
plot_read_distribution(distr_files,outpath_hisat,"filtered_out_readdistribution.png", yli = 20000000, lanes = TRUE)


  


distr_files <- iterate_files(inpath_hisat, pattern_distr)

distr_files_2 <- iterate_files(inpath_hisat_batch, pattern_distr)
distr_files_2 <- distr_files_2[!str_detect(distr_files_2,"nvsf")]
distr_files_2 <- distr_files_2[!str_detect(distr_files_2,"filter")]

plot_read_distribution_two(distr_files,distr_files_2,outpath_hisat,"filtered_out_readdistribution_2_tog.png", yli = 80000000)



distr_files_nvsf <-distr_files[str_detect(distr_files,"read_distribution_nvsf")]
plot_read_distribution(distr_files_nvsf,outpath_hisat,"filtered_out_readdistribution.png", yli = 20000000)






