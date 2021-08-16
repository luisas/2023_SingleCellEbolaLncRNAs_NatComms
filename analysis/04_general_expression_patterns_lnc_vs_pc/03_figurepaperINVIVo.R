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

source(file.path("../utils/00_datapaths.R"))
source("../utils/02_sc_utils.R")

# 0. Palette used throughout the scripts
col_lnc = "#f72631"
col_mrna = "#3153a3"
palette_plot_percentage <- c(col_lnc, col_mrna)

# 1. Theme for plots
theme_sc_paper <- theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(legend.text = element_text(size=10), legend.title = element_blank())+theme_paper

theme_matching <- theme(panel.background = element_rect(fill = "white", color = "black"))+
  theme(legend.position ="",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length=unit(.2, "cm"))+
  theme(axis.text = element_text(size = 15, color ="black"), axis.title = element_text(size = 17))

# 2. Directories
# Where stats about lnc and pc computed on cluster are stored
robjectsdir<- file.path(data_path, "/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/")

# 3. Load stats
df_complete_celltype <- readRDS(file.path(robjectsdir, "df_celltype.rds"))
df_complete <- readRDS(file.path(robjectsdir, "df_complete.rds"))
df_lnc <- readRDS(file.path(robjectsdir, "df_lnc.rds"))
df_mrna <- readRDS(file.path(robjectsdir, "df_mrna.rds"))

# Apply same filters
df_lnc <- df_lnc[df_lnc$n_cells > 30, ]
df_mrna <- df_mrna[df_mrna$n_cells > 30, ]
df_lnc$type <- "lncRNA"
df_mrna$type <- "Protein Coding"


# Import Seurat object 
file <- file.path(data_path, "02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds")
orthologs <- readRDS(file.path(data_path, "/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/orthologs_geneid.rds"))
immune.combined <- readRDS(file)

candidates  <- df_lnc[ df_lnc$medianexpr>1.5 & df_lnc$n_cells > 50& df_lnc$n_cells < 1200,]
gene <- candidates$gene_id[11]
print(gene)
# ----------------------------------------------
# 1 . Find the matching (median expr) Gene
# ----------------------------------------------
df_lnc_filtered <- df_lnc[df_lnc$gene_id == gene,]
df_complete<- rbind(df_lnc_filtered, df_mrna)
c <- df_complete[complete.cases(df_complete),]
c$type = ifelse(c$type == "lnc", 1, 0)
set.seed(125)
mi <- matchit(type ~ medianexpr,c)
matches <- get_matches(mi, data = c)
rownames(matches) <- gsub("-unknown", "", matches$gene_id)
matches$type = ifelse(matches$type == 0, "mrna", "lnc")
genes <- as.character(matches$gene_id)
matches[,c("medianexpr", "n_cells", "type")]


# ----------------------------------------------
#  2. Visualize that the expression of the 2 genes BOXPLOT
# ----------------------------------------------
expr_matrix <- immune.combined@assays$RNA@data
df_genes <- data.frame(t(expr_matrix[genes, ]))
colnames(df_genes) <-  gsub(".unknown", "", colnames(df_genes))
df_genes_melted <- melt(df_genes)
df_genes_melted <- df_genes_melted[df_genes_melted$value > 0, ]
bp <- ggboxplot(df_genes_melted, x = "variable", y = "value", fill = "variable", col = "black",alpha=1)+theme()+scale_fill_manual(values = rev(palette_plot_percentage))+scale_color_manual(values = rev(palette_plot_percentage))+stat_compare_means(comparisons = list(c(colnames(df_genes))))+xlab("")+ylab("logCP10K")+theme(text = element_text(size=15))+scale_y_continuous(limits = c(0,5), expand = c(0,0))+theme(axis.line.x = element_blank(), legend.position = "", axis.ticks.x = element_blank(), axis.text.x = element_blank())

pdf(file.path(plots, "02/genes_boxplot.pdf"), width = 6, height = 6)
bp
dev.off()
# ----------------------------------------------
#  2. Visualize that the expression of the 2 genes FEATURE PLOT
# ----------------------------------------------
f1 <- FeaturePlot(immune.combined,label = FALSE, features = c(genes[1]), pt.size = 0.6, order = TRUE,max.cutoff = 3.5,  cols = c("lightgrey", "navy"))+theme(axis.title.x.top = element_text())+theme_minimal()+ theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title = element_text(size=15), plot.title = element_text(size =22, hjust = 0.5))+labs(title = genes[1])+theme_void()+theme(title = element_blank(), text = element_text(size = 20))
f1
f2 <- FeaturePlot(immune.combined,label = FALSE, features = c(genes[2]), pt.size = 0.6, order = TRUE, max.cutoff = 3.5,cols = c("lightgrey", "#B0052D"))+theme(axis.title.x.top = element_text())+theme_minimal()+ theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title = element_text(size=15), plot.title = element_text(size =22, hjust = 0.5))+labs(title = gsub("-unknown", "", genes[2]))+theme_void()+theme(title = element_blank(), text = element_text(size = 20))

CombinePlots(list(f1,f2))

pdf(file.path(plots, "02/gene_1.pdf"), width = 6, height = 6)
f1
dev.off()


pdf(file.path(plots, "02/gene_2.pdf"), width = 6, height = 6)
f2
dev.off()