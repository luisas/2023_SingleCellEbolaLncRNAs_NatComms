library(Seurat)
library(ggplot2)
library(ggthemes)


# 0. Palette used throughout the scripts
col_lnc = "#B0052D"
col_mrna = "navy"
palette_plot_percentage <- c(col_lnc, col_mrna)


# 1. READ IN FILES 

robjectsdir <- file.path(data_path,"02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/")
all_lncrnas <- readRDS(file.path(robjectsdir, "all_lncrnas.rds"))
annotated_mrnas <- readRDS(file.path(robjectsdir,"annotated_mrnas.rds"))

immune.combined <- readRDS(file.path(data_path,"02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds"))
specificity_scores <- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/00_specificity/04_specificity_alternativescore.rds"))
table_summary_scores <- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/00_specificity/table_summary_alternative.rds"))
table_summary_scores$subtype <- "annotated"
table_summary_scores[substr(table_summary_scores$gene,1,4) == "MSTR", ]$subtype <- "novel"
table_summary_scores$type_subtype <- paste(table_summary_scores$type, table_summary_scores$subtype, sep = "_")

# Select genes : TO BE CHANGED
gene1 <- as.character(table_summary_scores[2,1])
gene2 <- as.character(table_summary_scores[3,1])
genes <- c(gene1, gene2)


# Lnc annotetad CELL TYPE SPEC
annotated_lnc <- table_summary_scores[table_summary_scores$type == "lnc" & table_summary_scores$subtype == "annotated",]
novel_lnc <- table_summary_scores[table_summary_scores$type == "lnc" & table_summary_scores$subtype == "novel",]
annotated_lnc_celltypespec <- rownames(annotated_lnc[order(annotated_lnc$score, decreasing = T),][1,])
novel_lnc_celltypespec <- rownames(novel_lnc[order(novel_lnc$score, decreasing = T),][7,])
novel_lnc_ubiq <- rownames(novel_lnc[order(novel_lnc$score, decreasing = F),][9,])
annot_lnc_ubiq <- rownames(annotated_lnc[order(annotated_lnc$score, decreasing = F),][2,])



fplot <- function(gene, col){
  f3 <- FeaturePlot(immune.combined,label = FALSE, features = c(gene), pt.size = 0.6, order = TRUE, cols = c("lightgrey", "#B0052D"))+theme(axis.title.x.top = element_text())+theme_minimal()+ theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(axis.title = element_text(size=15), plot.title = element_text(size =22, hjust = 0.5))+
    labs(title = get_orthologname(gene))+theme_void()+theme(text = element_text(size = 20), plot.title = element_text(size =22, hjust = 0.5))
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

pdf(file.path(plots,"02/GENE_novel_ubi.pdf"), width = 6, height = 6)
fplot(novel_lnc_ubiq, col_lnc)
dev.off()


pdf(file.path(plots,"02/GENE_annot_ubi.pdf"), width = 6, height = 6)
fplot(annot_lnc_ubiq, col_lnc)
dev.off()



# Check the distribution of cell-type specificity scores


ggplot(table_summary_scores, aes(y = score, fill = type_subtype, color  = type_subtype))+geom_boxplot(alpha = 0.5)+theme_classic2()

wilcox.test(table_summary_scores[table_summary_scores$type_subtype == "lnc_novel",]$score , table_summary_scores[table_summary_scores$type_subtype == "lnc_annotated",]$score)
wilcox.test(table_summary_scores[table_summary_scores$type_subtype == "lnc_novel",]$score , table_summary_scores[table_summary_scores$type_subtype == "pc_annotated",]$score)
wilcox.test(table_summary_scores[table_summary_scores$type_subtype == "lnc_annotated",]$score , table_summary_scores[table_summary_scores$type_subtype == "pc_annotated",]$score)


length(unique(table_summary_scores[table_summary_scores$type_subtype == "lnc_novel",]$gene))
length(unique(table_summary_scores[table_summary_scores$type_subtype == "lnc_annotated",]$gene))
length(unique(table_summary_scores[table_summary_scores$type_subtype == "pc_annotated",]$gene))



table_summary_scores_celltype <- merge(table_summary_scores, specificity_scores[,c("gene", "identity")], by = "gene")

# Total number of cells per celltype
celltypes_ncells <- as.data.frame(table((Idents(immune.combined))))

table_summary_scores_celltype$identity <- factor(table_summary_scores_celltype$identity, levels  = c("T", "B", "Monocyte", "Neutrophil"))
ggplot(table_summary_scores_celltype, aes(y = score, col = identity, fill = identity))+geom_boxplot(alpha = 0.5)+theme_void()

p <- ggplot(table_summary_scores_celltype, aes(y = score, col = type, fill = type))+geom_boxplot(alpha = 0.5)+theme_void()

facet(p, facet.by = "identity")

