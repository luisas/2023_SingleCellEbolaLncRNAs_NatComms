library(igraph)
library(stringr)

source("../../utils/00_datapaths.R")
source("../../utils/02_sc_utils.R")
source("../../utils/04_utils_graph.R")

# 0. Load needed informations 
gD2 <- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/gD2.rds"))
lou2 <- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/lou.rds"))
linkList <- readRDS(file.path(data_path, "/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/linkList.rds"))
orthologs <- readRDS(file.path(data_path, "01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/orthologs_geneid_ready.rds"))
de_all<- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_stages.rds"))
de_lnc <- de_all[de_all$type == "lnc", ]
de_lnc_mono <- gsub("-unknown", "",de_all[de_all$celltype == "Monocyte" & de_all$type == "lnc",]$primerid)
graph_info <- readRDS(file.path(data_path, "/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/Graph_communitybelonging.rds"))
summary <- readRDS(file.path(data_path, "00_Metadata/final_correlation_summary.rds"))
immresp_genes <- readRDS(file.path(data_path, "/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/enrichment_genes_ImmRespRegu.rds"))
innate_genes <- readRDS(file.path(data_path, "/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/enrichment_genes_Innate.rds"))



# 1. Basic layout
V(gD2)$frame.color <- "black"
V(gD2)$community <- membership(lou2)
#l <- layout.fruchterman.reingold(gD2, niter=5000)
#saveRDS(l, "/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/layout.rds")
l <- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/layout.rds"))


# ------ General proprieties -------------------------
# 1. Dot's size 



# -------------------------------
#      GENES TO HIGHLIGHT
# -------------------------------
go_purine_metabolism <- c("GUCY1B1", "GUCY1A1", "ADA")
#go_leukocyte_proliferation <- c("FN1","CD74","CDKN1A","PSAP","ZFP36","GADD45B","DUSP6")
#go_response_to_other_organism <- c("S100A9","PLAC8", "S100A8","ISG15" ,"IFIT2", "MX1")
go_genes <- c(go_purine_metabolism,innate_genes,immresp_genes)
go_genes <- c(go_genes, "ETS1")

V(gD2)$size <- ifelse(names(as.list(V(gD2))) %in% c(de_lnc_mono,go_genes), 5, 1.5)
V(gD2)$shape <- ifelse(names(as.list(V(gD2))) %in% c(de_lnc_mono), "fcircle",shapes()[3] )
V(gD2)$bw <-ifelse(names(as.list(V(gD2))) %in% c(de_lnc_mono), 5.0,0)

# 2. Label size
V(gD2)$label <- ifelse(names(as.list(V(gD2))) %in% c(de_lnc_mono, go_genes), unlist(lapply(names(as.list(V(gD2))), get_orthologname_)), NA)
V(gD2)$font <- ifelse(names(as.list(V(gD2))) %in% c(de_lnc_mono), 1, 1 )
V(gD2)$degree <- -pi/2 
# Change position of specific genes 
# Below
V(gD2)$degree <- ifelse(names(as.list(V(gD2))) %in% c("IL7R", "S100A9", "S100A8","GADD45B", "ISG15", "MSTRG.195139", "MS4A7", "LGALS3", "ENSMMUG00000058532"), pi/2, V(gD2)$degree )
V(gD2)$degree <- ifelse(names(as.list(V(gD2))) %in% c("FCGR3" ), pi/1.4, V(gD2)$degree )

# Right
V(gD2)$degree <- ifelse(names(as.list(V(gD2))) %in% c("CDKN1A", "GADD45B", "HBEGF", "MSTRG.168139", "LINC00861", "ENSMMUG00000052424", "SRRM2", "ENSMMUG00000064224"), 0, V(gD2)$degree )
V(gD2)$dist <- ifelse(names(as.list(V(gD2))) %in% c("CDKN1A", "GADD45B","HBEGF", "LINC00861", "ENSMMUG00000052424", "SRRM2"), 1.5,0.55 )
V(gD2)$dist <- ifelse(names(as.list(V(gD2))) %in% c("ENSMMUG00000064224"), 3.1,V(gD2)$dist )

# Left
V(gD2)$degree <- ifelse(names(as.list(V(gD2))) %in% c("ENSMMUG00000058644", "MSTRG.228510", "MSTRG.17705", "GUCY1A1", "PSAP","LINC00861"), pi, V(gD2)$degree )
V(gD2)$dist <- ifelse(names(as.list(V(gD2))) %in% c("ENSMMUG00000058644","GUCY1A1", "PSAP","LINC00861" ), 3.5,V(gD2)$dist)
V(gD2)$dist <- ifelse(names(as.list(V(gD2))) %in% c( "MSTRG.228510", "MSTRG.17705", "MSTRG.168139"), 2.5,V(gD2)$dist)
V(gD2)$dist <- ifelse(names(as.list(V(gD2))) %in% c( "LINC00861","GUCY1A1", "PSAP"), 1.5,V(gD2)$dist)

V(gD2)$degree <- ifelse(names(as.list(V(gD2))) %in% c("KLF4"), pi, V(gD2)$degree )
V(gD2)$dist <- ifelse(names(as.list(V(gD2))) %in% c("KLF4" ), 1,V(gD2)$dist)



# ---- 1. Get the direction ( UP, DOWN or MIXED)
genes <- names(as.list(V(gD2)))

de_all$id <- gsub("-unknown","",de_all$primerid) 
de_all$id <- gsub("MAMU-", "MAMU.", de_all$id, fixed = T)
de_mono <- de_all[de_all$celltype == "Monocyte",]

de_graph <- de_mono %>% dplyr::group_by(id) %>%  dplyr::filter(abs(coef) == max(abs(coef)))
de_graph <- de_graph[,c("primerid", "direction", "type", "stage", "id")]
mono_de <-de_graph

#--- Compute color Up/DOWN
graph_info$id <- graph_info$gene_name
graph_info_de <- merge(graph_info, mono_de, by="id")
# Add community size 
comm_size <- graph_info_de %>% dplyr::group_by(Community) %>% dplyr::count()
names(comm_size) <- c("Community", "community_size")
graph_info_de <- merge(graph_info_de, comm_size, by="Community")
# % of UP 
graph_info_de_full <- graph_info_de %>% dplyr::group_by(Community,community_size, direction) %>% dplyr::count()
graph_info_de_up_perc <- graph_info_de_full[graph_info_de_full$direction == "up",]
graph_info_de_up_perc$perc_up <- graph_info_de_up_perc$n/ graph_info_de_up_perc$community_size
#Get color palette
f <- colorRamp(c("blue","grey","red"))
colors <- rgb(f(graph_info_de_up_perc$perc_up)/255)


# -------------------
# Change Color 

# 0. Default Color 
# 1. Up,down
V(gD2)$color <- ifelse(names(as.list(V(gD2))) %in% c(mono_de[mono_de$direction == "up", ]$id), "red", "#6989F5" )

# 2. Down, shades
mono_de$color <- "grey"
#mono_de$stage <- as.character(mono_de$stage)
#mono_de[mono_de$stage == "early", ]$color <- "white"
#mono_de[mono_de$stage == "middle", ]$color <- "#6C6B6B"
#mono_de[mono_de$stage == "late", ]$color <- "black"


#mono_de$size <- 3
#mono_de[mono_de$stage == "early", ]$size <- 1
#mono_de[mono_de$stage == "middle", ]$size <- 2
#mono_de[mono_de$stage == "late", ]$size <- 3
rownames(mono_de) <- mono_de$id
#V(gD2)$color  <- mono_de[names(as.list(V(gD2))),]$color


mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                 plot=mycircle, parameters=list(vertex.frame.color=1,
                                                vertex.frame.width=1))


# -------------------Plot the Graph 
png(filename = file.path(data_path,"plots/04/Graph_main_UpDown.png"), width = 3300, height = 2100)
plot(lou2, gD2, frame.width = 20, 
     edge.arrow.size = .1,
     col=V(gD2)$color,
     layout=l,
     edge.color= "black",
     mark.col=adjustcolor("grey", alpha = 0.3),
     mark.border=adjustcolor("grey", alpha = 0.3),
     vertex.label.color="black", 
     vertex.label.font=V(gD2)$font,
     vertex.label.cex=2.5,
     vertex.label.dist=V(gD2)$dist , 
     vertex.label.degree=V(gD2)$degree,
     vertex.label.family="Helvetica",
     vertex.label = V(gD2)$label,
     vertex.shape=V(gD2)$shape,
     vertex.frame.width=V(gD2)$bw,
     vertex.frame.color = "black")

dev.off()






