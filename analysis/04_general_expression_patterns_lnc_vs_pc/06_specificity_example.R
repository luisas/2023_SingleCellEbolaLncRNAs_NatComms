library(Seurat)
library(ggplot2)
library(ggthemes)


# 0. Palette used throughout the scripts
col_lnc <- "#B0052D"
col_mrna = "navy"
palette_plot_percentage <- c(col_lnc, col_mrna)


# 1. READ IN FILES 

robjectsdir <- file.path(data_path,"02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/")
all_lncrnas <- readRDS(file.path(robjectsdir, "all_lncrnas.rds"))
annotated_mrnas <- readRDS(file.path(robjectsdir,"annotated_mrnas.rds"))
marker.genes <- readRDS(file.path(data_path, "02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/marker.genes.rds"))
immune.combined <- readRDS(file.path(data_path,"02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds"))
specificity_scores <- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/00_specificity/04_specificity_alternativescore.rds"))
table_summary_scores <- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/00_specificity/table_summary_alternative.rds"))
table_summary_scores$subtype <- "annotated"
table_summary_scores[substr(table_summary_scores$gene,1,4) == "MSTR", ]$subtype <- "novel"
table_summary_scores$type_subtype <- paste(table_summary_scores$type, table_summary_scores$subtype, sep = "_")
annotated_lnc <- table_summary_scores[table_summary_scores$type == "lnc" & table_summary_scores$subtype == "annotated",]
novel_lnc <- table_summary_scores[table_summary_scores$type == "lnc" & table_summary_scores$subtype == "novel",]



# SELECT GENES TO PLOT
# 

annotated_lnc_celltypespec <- rownames(annotated_lnc[order(annotated_lnc$score, decreasing = T),][1,])
novel_lnc_celltypespec <- rownames(novel_lnc[order(novel_lnc$score, decreasing = T),][7,])
# Matching pc genes in cell-type specificity score 
table_summary_scores[novel_lnc_celltypespec,]
pc_celltype_spec <- "RORA"
table_summary_scores[table_summary_scores$score > 0.89 & table_summary_scores$score < 0.91 & table_summary_scores$n_cells > 2500,]
novel_lnc_ubiq <- rownames(novel_lnc[order(novel_lnc$score, decreasing = F),][9,])
annot_lnc_ubiq <- rownames(annotated_lnc[order(annotated_lnc$score, decreasing = F),][2,])


t_spec <- marker.genes[3]

# Function to plot genes 
fplot <- function(gene, col){
  f3 <- FeaturePlot(immune.combined,label = FALSE, features = c(gene), pt.size = 0.6, order = TRUE, cols = c("lightgrey", col))+theme(axis.title.x.top = element_text())+theme_minimal()+ theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title = element_text(size=15), plot.title = element_text(size =22, hjust = 0.5))+
    labs(title = get_orthologname_(gene))+theme_void()+theme(text = element_text(size = 20), plot.title = element_text(size =22, hjust = 0.5))
  return(f3)
}
# ----------------------------------------------
#  2. Visualize that the expression of the 2 genes FEATURE PLOT
# ----------------------------------------------



pdf(file.path(plots,"02/GENE_annot_spec.pdf"), width = 6, height = 6)
fplot(annotated_lnc_celltypespec, col_lnc)
dev.off()


pdf(file.path(plots,"02/GENE_novel_spec.pdf"), width = 6, height = 6)
fplot(novel_lnc_celltypespec, col_lnc)
dev.off()
table_summary_scores[novel_lnc_celltypespec,]


pdf(file.path(plots,"02/GENE_PC_spec.pdf"), width = 6, height = 6)
fplot(pc_celltype_spec, col_mrna)
dev.off()
table_summary_scores[pc_celltype_spec,]

pdf(file.path(plots,"02/GENE_novel_ubi.pdf"), width = 6, height = 6)
fplot(novel_lnc_ubiq, col_lnc)
dev.off()


pdf(file.path(plots,"02/GENE_annot_ubi.pdf"), width = 6, height = 6)
fplot(annot_lnc_ubiq, col_lnc)
dev.off()

pdf(file.path(plots,"02/GENE_pc_spec.pdf"), width = 6, height = 6)
fplot(t_spec, col_mrna)
dev.off()


