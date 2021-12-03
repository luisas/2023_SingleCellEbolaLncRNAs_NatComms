library(igraph)
library(stringr)
library(Seurat)

source("../../utils/00_datapaths.R")
source("../../utils/02_sc_utils.R")
source("../../utils/04_utils_graph.R")

# 0. Load needed informations 
immune.combined <- readRDS(file.path(data_path, "/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds"))
dim(immune.combined)
monocyte <- immune.combined[,Idents(immune.combined)=="Monocyte"]


monocyte$group_red <- factor(toupper(substring(monocyte$group, 1,1)), levels = c("B", "E", "M", "L"))

# select genes 
lnc_response_to_other_organism <- c("ENSMMUG00000064224-unknown", "ENSMMUG00000062255-unknown")
lnc_leukocyte_proliferation <- c("MSTRG.181325-unknown", "MIR22HG", "ENSMMUG00000058644-unknown")
go_leukocyte_proliferation <- c("FN1","CD74","CDKN1A","PSAP","ZFP36","GADD45B","DUSP6")
go_response_to_other_organism <- c("S100A9","PLAC8", "S100A8","ISG15" ,"IFIT2", "MX1")

response_to_other_organism <-c(lnc_response_to_other_organism, go_response_to_other_organism)

# 1.  Dotplot response to other organism 
#scale.func <- switch(EXPR = "radius", size = scale_size, 
#                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
#plot.data <- Dotplot_data(monocyte, features = genes, dot.scale =25,group.by = "group",scale.by = "radius", cols = c( "blue", "red"))
#plot.data$name <- unlist(lapply(plot.data$features.plot, get_orthologname_))
#dotplot_isgs <- ggplot(plot.data, aes(x = id,y = name, size= pct.exp, col = avg.exp.scaled))+geom_point()+theme_classic()+xlab("")+ylab("")+
#        scale_colour_gradient(low = "blue", high = "red", na.value = NA)+
#        scale.func(range = c(0, 13))+ggtitle("")+
#        theme(axis.text.x= element_text(size = 20, color = "black"), axis.text.y= element_text(size = 15, color = "black"))


plot.data_up <- Reduce(rbind, lapply(go_response_to_other_organism, function(gene) plot_expression_gene(gene, monocyte)[3][[1]]))
plot.data_up$name <- unlist(lapply(plot.data_up$gene, get_orthologname_))
dotplot_isgs <- ggplot(plot.data_up, aes(x = type, y = name, size= pct.exp, col = mean))+geom_point()+theme_paper+xlab("")+scale_colour_gradient(low = "blue", high = "red", na.value = NA)+scale_size(range = c(5,10))+
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank() )
pdf(file.path(plots, "04/dotplot_isgs.pdf"), width = 8, height = 6)
dotplot_isgs
dev.off()

