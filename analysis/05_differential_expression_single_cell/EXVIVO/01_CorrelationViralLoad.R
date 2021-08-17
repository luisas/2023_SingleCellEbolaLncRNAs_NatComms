library(dplyr)
library(Seurat)
library(Matrix)
library(SingleCellExperiment)
library(stringr)
library(rtracklayer)
library(RColorBrewer)
library(scales)
library(ggthemes)
library(ggplot2)


# Define paths for data
source("../utils/00_datapaths.R")
source("../utils/02_sc_utils.R")
source("../utils/04_utils_graph.R")

# Load needed files EX VIVO 
orthologs <- readRDS(file.path(data_path, "/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/orthologs_geneid.rds"))
mono_live_h24_inf <- readRDS( file.path(data_path,"02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/03_prep/03_monocytes_infected_late.rds"))
robjectsdir_stats <- file.path(data_path,"02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/05_stats/")
all_lncrnas <- readRDS(file.path(robjectsdir_stats, "all_lncrnas.rds"))
annotated_mrnas <- readRDS(file.path(robjectsdir_stats,"annotated_mrnas.rds"))
ebola_genome_percentage_df <- readRDS(file.path(data_path, "02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/03_prep/df_viralpercentage.rds"))


# --------------------------------
# 1. LOAD and prep CORRELATIONS
RHOTHRESHOLD = 0.0
FDRTHRESHOLD<- 0.05
correlations_infected_24 <- readRDS(file.path(data_path, "02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/03_viralload_infected_late/ViralLoad_spearman-correlations_standard.rds"))
correlations_infected_24 <- correlations_infected_24[!is.na(correlations_infected_24$pval),]
sig_correlations_pc_infected_24 <- correlations_infected_24[!is.na(correlations_infected_24$pval) & correlations_infected_24$pval < FDRTHRESHOLD & abs(correlations_infected_24$rho) > RHOTHRESHOLD, ]

# Add Gene Type
sig_correlations_pc_infected_24$type <- "-"
sig_correlations_pc_infected_24[sig_correlations_pc_infected_24$gene %in% all_lncrnas ,]$type <- "lnc"
sig_correlations_pc_infected_24[sig_correlations_pc_infected_24$gene %in% annotated_mrnas ,]$type <- "pc"
sig_cor_lnc <- sig_correlations_pc_infected_24[sig_correlations_pc_infected_24$type == "lnc",]
sig_cor_lnc$gene_name <- unlist(lapply(sig_cor_lnc$gene, get_orthologname_))
sig_cor_lnc <- sig_cor_lnc[abs(sig_cor_lnc$rho) > 0.1,]
# remove ebola gene 
sig_cor_lnc <- sig_cor_lnc[sig_cor_lnc$gene != "MSTRG.252489-unknown",]


# --------------------------------
# 2. PLOT MAIN 
cols <- c('highlight' = 'red', 'other' = 'grey')
sig_cor_lnc
genes_expression <- Reduce("rbind", lapply(sig_cor_lnc$gene, function(gene) get_expression_summary_gene(gene, mono_live_h24_inf)))
pdf(file.path(plots, "05/main_corr.pdf"), width = 7, height = 5)
# 2. Plot 
ggplot(genes_expression, aes(x = percentage_viral_reads, y = value, col = gene))+stat_smooth(method = "loess", formula = y ~ x, size = 0.4, se = F, n = 7, span = 0.5)+
  theme_classic()+
  theme( legend.position = "", text = element_text(size = 17))+ 
  xlab("viral load")+ylab("gene expression (logCP10K)")+scale_x_log10()+
  scale_color_manual(values = c(rep("#848484", length(unique(genes_expression$gene)))))
dev.off()




