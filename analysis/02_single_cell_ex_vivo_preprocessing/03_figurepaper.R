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
source(file.path("../utils/00_datapaths.R"))
source("../utils/02_sc_utils.R")

theme_umap <- theme(panel.background = element_rect(fill = "white"),
                    panel.grid.major = element_blank(),
                    legend.position = "right", 
                    panel.grid.minor = element_blank(),
                    text = element_text(size=18))
# Gene annotation 
ref <- import(file.path(data_path,"01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf"))
ebola_ref <- import(file.path(data_path,"00_RawData/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.gtf"))

immune.combined <- readRDS(file.path(data_path, "02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10//05_RObjects/03_prep/immune.combined.infectionstatus.rds"))

# Import marker genes
marker.genes_red <- readRDS(file.path(data_path,  "02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/03_prep/marker.genes.rds"))
# Define palettes
pal_celltypes <- c("#FD6467","#F1BB7B","#1F78B4","#AE4E4E")
pal4 <- c("#EBCC2A","#E55039","#3B9AB2")
pal3 <- c( "#3B9AB2","#78B7C5","#EBCC2A")

# 1. Celltypes 
pdf(file.path(plots, "05/A_celltypes.pdf"), width = 7, height = 5)
DimPlot(immune.combined, reduction = "umap", label = TRUE, cols =pal_celltypes, label.size = 9)+theme_void()+theme_umap
dev.off()

# 2. Dotplot
pdf(file.path(plots, "05/B_dotplot.pdf"), width = 8, height = 5)
DotPlot(immune.combined, features = unique(marker.genes_red), cols = c("grey", "dark red", "white"), dot.scale = 11) + RotatedAxis()+theme(axis.text = element_text(size = 20), axis.title = element_blank())+ theme(panel.background = element_rect(fill = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(legend.text = element_text(size=21))
dev.off()

# 3. SUPPL : DPI 
pdf(file.path(plots, "05/SUPPL_DPI.pdf"), width = 7, height = 5)
# Day post infection 
colfunc <- colorRampPalette(c("black", "white"))
DimPlot(immune.combined, red = "umap", group.by = "dpi", pt.size = 0.1, cols = rev(brewer.pal(8, "Paired")))+theme_umap
dev.off()

# 3. SUPPL : Individual
pdf(file.path(plots, "05/SUPPL_individual.pdf"), width = 7, height = 5)
DimPlot(immune.combined, reduction = "umap", group.by = "individual_nhp",cols = c(brewer.pal(10, "Paired")[c(1,2)]), pt.size = 0.1)+theme_umap
dev.off()

# 3: SUPPL: Condition 
pdf(file.path(plots, "05/SUPPL_hours.pdf"), width = 7, height = 5)
p4 <- DimPlot(immune.combined, reduction = "umap", group.by = "cond",cols =pal4, pt.size = 0.01)+theme_umap
p4
dev.off()

# ------------- infected cells ---------

expression_matrix <- immune.combined@assays$RNA@counts
colnames(expression_matrix) <- Idents(immune.combined)

ebola_genes <- readRDS(file.path(data_path, "02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/05_stats/ebola_genes.rds"))
# Per cell, check the % of ebola reads ( not normalized - just calculate the % of raw counts )
ebola_reads <- expression_matrix[rownames(expression_matrix) %in% ebola_genes, ]
ebola_genome_reads <-colSums(ebola_reads)

total_reads <- colSums(expression_matrix)
ebola_gene_percentage <- (ebola_reads/total_reads)*100
ebola_genome_percentage<- (ebola_genome_reads/total_reads)*100

ebola_genome_percentage_df <- data.frame(percentage_viral_reads = ebola_genome_percentage,
                                         celltype = names(ebola_genome_percentage),
                                         cond =immune.combined$cond, 
                                         hour = immune.combined$dpi)


# Define the threshold for the percentage of viral reads identified per cell for it to be called Infected 
#threshold <- quantile(ebola_genome_percentage_df[ebola_genome_percentage_df$celltype != "Myeloid", ]$percentage_viral_reads, probs = c(0.99))
threshold <- max(unlist(lapply(setdiff(unique(ebola_genome_percentage_df$celltype), "Monocyte"), function(x) quantile(ebola_genome_percentage_df[ebola_genome_percentage_df$celltype == x, ]$percentage_viral_reads, probs = c(0.99)))))



# Load file with information about infected cells
ebola_genome_percentage_df <- readRDS(file.path(data_path, "02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/03_prep/df_viralpercentage.rds"))



pdf(file.path(plots, "05/SUPPL_density.pdf"), width = 7, height = 5)
ggplot(ebola_genome_percentage_df, aes( x = percentage_viral_reads , col = classification, fill = classification))+geom_density( alpha = 0.6)+theme_paper+xlab("")+theme(axis.text.y = element_text(size = 8))+theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.y = element_text(size = 15), axis.text.y = element_text(size =15))+scale_x_log10()+scale_fill_manual(values = c( "red", "grey"))+scale_color_manual(values = c( "red", "grey"))+xlab("Viral Load")+geom_vline(xintercept = threshold, lty =  2)
dev.off()


pdf(file.path(plots, "05/Infectedcells.pdf"), width = 5, height = 5)
# Plot infected cells
infected_cells <- rownames(ebola_genome_percentage_df[ebola_genome_percentage_df$classification == "Infected",])
# Or with labels (Same plot just with cell-type labels)
DimPlot(immune.combined, label=F, cells.highlight= list(infected_cells), sizes.highlight = 0.4, cols.highlight = c("#DA0202"), cols= "grey", pt.size = 0.1)+ 
  labs(title = "Infected cells" )+ NoLegend()+theme_minimal()+ theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(legend.text = element_text(size=18))+
  theme(legend.position = "none", plot.title = element_text(size = 22), axis.title = element_text(size = 18))
dev.off()



summary_celltype <- ebola_genome_percentage_df%>% dplyr::group_by(celltype) %>% tally()
ebola_genome_percentage_df$celltype_count <- unlist(apply(ebola_genome_percentage_df,1,function(x) summary_celltype[summary_celltype$celltype ==x[2],]$n))
# Visualize the percentage of infected cells per celltype
percentage_cells_infected_per_celltype <- ebola_genome_percentage_df %>% dplyr::group_by(celltype, classification, celltype_count, cond, hour) %>% tally()

percentage_cells_infected_per_celltype <- percentage_cells_infected_per_celltype[percentage_cells_infected_per_celltype$classification == "Infected", ]

percentage_cells_infected_per_celltype$perc_infected <- (percentage_cells_infected_per_celltype$n / percentage_cells_infected_per_celltype$celltype_count)*100
# reset levels in right order
percentage_cells_infected_per_celltype$cond <- factor(percentage_cells_infected_per_celltype$cond , levels = c("media", "irrad", "live"))


pdf(file.path(plots, "05/barplot.pdf"), width = 5, height = 5)
percentage_cells_infected_per_celltype$condhr <- paste(percentage_cells_infected_per_celltype$cond, percentage_cells_infected_per_celltype$hour)
ggplot(percentage_cells_infected_per_celltype, aes( x= celltype, y = perc_infected , col = condhr, fill = condhr))+geom_bar(stat="identity",position = position_dodge(), alpha = 0.8)+theme_paper+xlab("")+ylab("% of cells infected")+theme(axis.text.y = element_text(size = 8))+scale_fill_manual(values = pal3)+scale_color_manual(values =pal3)+theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.y = element_text(size = 15), axis.text.y = element_text(size =15))+geom_hline(yintercept =1, lty =2)
dev.off()

