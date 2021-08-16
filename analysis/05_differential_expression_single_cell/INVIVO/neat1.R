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



theme_paper <- theme(legend.title = element_blank())+theme(panel.background = element_rect(fill = "white", colour = "white"))+theme(panel.background = element_rect(fill = "white", colour = "grey50"))+theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 18))

# Imports
source("../../utils/00_datapaths.R")

clusterpath <- "/gpfs/projects/bsc83/Data/"
gene_annotation_path <- file.path(clusterpath, "gene_annotation")
data_path <- file.path(clusterpath, "Ebola")
plots <- file.path(data_path,"plots")
source("../../utils/02_sc_utils.R")
source("../../utils/03_de_utils.R")
source("../../utils/04_utils_graph.R")
ref <- import(file.path(data_path,"01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf"))

# Seurat object
immune.combined <- readRDS(file.path(data_path, "02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds"))

# Corrected q-value .05
FDRTHRESHOLD<- 0.05
# Fold-Change of 30%
FCTHRESHOLD <- log(1.23)

orthologs <- readRDS(file.path(data_path, "01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/orthologs_geneid_ready.rds"))
mono <- subset(immune.combined, ident = "Monocyte")

# ----------------------------------------------------
# 1. Check ortholog 
# ----------------------------------------------------

neat1 <- get_gene_id("NEAT1")

infection_info <- readRDS(file.path(data_path,"02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/df_viralpercentage.rds"))

table(colnames(immune.combined) == rownames(infection_info))
immune.combined$viralperc <- infection_info$percentage_viral_reads
immune.combined$infectionstatus <- infection_info$classification



# ----------------------------------------------------
# 2. Visualize downregulation general 
# ----------------------------------------------------
print("Here")
plot.data_up <- Dotplot_data(mono, features = c(neat1 ), group.by = "group",scale.by = "radius", dot.scale = 20,cols = c( "lightblue", "darkred"))
plot.data_up$name <- unlist(lapply(plot.data_up$features.plot, get_orthologname_))
#ggplot(plot.data_up, aes(x = id,y = name, size= pct.exp, col = avg.exp.scaled))+geom_point()+theme_paper+xlab("")+xlab("")+scale_colour_gradient(low = "lightblue", high = "dark red", na.value = NA)+theme(axis.text.x = element_text(angle = 20, vjust = 0.9, hjust=1))
#Heatmap
print("HH")
plot.data_up$id <- gsub("late", "ulate", as.character(plot.data_up$id))

pdf(file.path(plots, "05/neat1.pdf"), width = 3, height = 2)
heatmap_celltypes(plot.data_up)
dev.off()
