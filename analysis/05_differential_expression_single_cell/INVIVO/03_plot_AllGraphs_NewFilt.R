library(igraph)
library(stringr)
source("../../utils/04_utils_graph.R")

# 0. Load needed informations 
gD2 <- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/gD2_newFilt_0.1.rds"))
lou2 <- readRDS(file.path(data_path,"/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/lou_newFilt_0.1.rds"))
linkList <- readRDS(file.path(data_path, "/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/linkList_NewFilt_0.1.rds"))
orthologs <- readRDS(file.path(data_path, "01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/orthologs_geneid_ready.rds"))
de_all<- readRDS(paste0(robjectsdir, "04_DE/allCells_DE_table.rds"))
de_all <- de_all[abs(de_all$logFC > 0.1) | de_all$fdr < 0.05, ]
de_lnc <- de_all[de_all$type == "lnc", ]
de_lnc_mono <- gsub("-unknown", "",de_all[de_all$celltype == "Monocyte" & de_all$gene_biotype %in% c("lnc", "novel_lnc"),]$gene_id)
graph_info <- readRDS(file.path(data_path, "/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/Graph_communitybelonging_NewFilt_0.1.rds"))
#summary <- readRDS(file.path(data_path, "00_Metadata/final_correlation_summary.rds"))


#Layout 
# 0. Basic layout
V(gD2)$frame.color <- "black"
V(gD2)$community <- membership(lou2)
l <- layout.fruchterman.reingold(gD2, niter=5000)

# Genes driving the enrichemnts (to highlight )
go_immune_resp <- c("S100A8", "S100A9", "IFIT2", "IDO1", "MX1", "ISG15")
go_cell_prol <- c("FOS", "LYZ","HSPA5", "CD74", "CDKN1A", "ANXA1", "MPEG1")
go_lps <- c("PPBP", "ZFP36")


lnc_communities <- lou2$names[lou2$membership %in% c(3, 5, 7)][lou2$names %in% de_lnc_mono]

lnc_all <- de_lnc_mono

go_genes <- c(go_immune_resp, go_cell_prol, go_lps)
V(gD2)$label <- ifelse(names(as.list(V(gD2))) %in% c(go_genes, de_lnc_mono), unlist(lapply(names(as.list(V(gD2))), get_orthologname_)), NA)
#V(gD2)$label <- NA
V(gD2)$size <- ifelse(names(as.list(V(gD2))) %in% c(de_lnc_mono,go_genes), 5, 1.5)
V(gD2)$shape <- ifelse(names(as.list(V(gD2))) %in% c(de_lnc_mono), "fcircle",shapes()[3] )
V(gD2)$bw <-ifelse(names(as.list(V(gD2))) %in% c(de_lnc_mono), 5.0,0)


####### NETWORK COLORED BY UP/DOWN GENES 

# ---- 1. Get the direction ( UP, DOWN or MIXED)
genes <- names(as.list(V(gD2)))
de_all$id <- gsub("-unknown","",de_all$gene_name) 
de_all$id <- gsub("MAMU-", "MAMU.", de_all$gene_name, fixed = T)
de_mono <- de_all[de_all$celltype == "Monocyte",]

de_directionality <- de_mono %>% dplyr::group_by(id) %>%  dplyr::filter(abs(logFC) == max(abs(logFC)))
de_directionality <- de_directionality[,c("gene_name", "direction", "gene_biotype", "id")]
mono_de <-de_directionality

#--- 2. Compute color Up/DOWN
# 1. Up,down
V(gD2)$color <- NA
de_directionality$color <- NA
de_directionality[de_directionality$direction == "up", ]$color <- "red"
de_directionality[de_directionality$direction == "down", ]$color <- "#6989F5"
table(is.na(de_directionality$color))
rownames(de_directionality) <- de_directionality$id
rownames(de_directionality)  <- gsub("-unknown", "", rownames(de_directionality))
V(gD2)$color  <- de_directionality[names(as.list(V(gD2))),]$color
table(is.na(V(gD2)$color))


mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/300 * params("vertex", "size")
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
png(filename =file.path(data_path,"Graph_main.png"), width = 3300, height = 2100)
plot(lou2, gD2, frame.width = 500000, 
     edge.arrow.size = 0.1,
     col=V(gD2)$color,
     layout=l,
     edge.color= "black",
     mark.col=adjustcolor("grey", alpha = 0.3), 
     mark.border=adjustcolor("grey", alpha = 0.3),
     vertex.label.color="black", 
     vertex.label.font=V(gD2)$font,
     vertex.label.cex=2.5,
     vertex.label.dist=V(gD2)$dist*0.9, 
     vertex.label.degree=V(gD2)$degree,
     vertex.label.family="Helvetica",
     vertex.label = V(gD2)$label,
     vertex.shape=V(gD2)$shape,
     vertex.frame.width=V(gD2)$bw*0.5,
     vertex.frame.color = "black")

dev.off()


pdf(file = file.path(data_path,"Graph_main_UpDown_withNames.pdf"))
plot(lou2, gD2,
     edge.arrow.size = .0001,
     col=V(gD2)$color,
     layout=l,
     edge.color= "darkgrey",
     mark.col=adjustcolor("grey", alpha = 0.3),
     mark.border=adjustcolor("grey", alpha = 0.3),
     vertex.label.color="#333333", 
     vertex.label.font=V(gD2)$font,
     vertex.label.cex=0.5,
     vertex.label.dist=V(gD2)$dist , 
     vertex.label.degree=V(gD2)$degree,
     vertex.label.family="Helvetica",
     vertex.label = V(gD2)$label,
     vertex.shape=V(gD2)$shape,
     vertex.frame.width=0.1,
     vertex.frame.color = "black")

dev.off()


####### NETWORK COLORED BY STAGES 

# ---- 1. Get the stage with the highest logFC 

genes <- names(as.list(V(gD2)))
de_all$ids <- gsub("-unknown","",de_all$gene_name) 
de_all$ids <- gsub("MAMU-", "MAMU.", de_all$gene_name, fixed = T)
de_mono <- de_all[de_all$celltype == "Monocyte",]


de_stages <- de_mono %>% dplyr::group_by(id) %>%  dplyr::filter(abs(logFC) == max(abs(logFC)))
de_stages <- de_stages[,c("gene_name", "stage","id")]
mono_de <-de_stages
mono_de$id <- gsub("-unknown","", mono_de$id)


#--- 2. Compute color early/mid/late 
# 1. Up,down
V(gD2)$color <- "grey"
mono_de$color <- "grey"
mono_de$stage <- as.character(mono_de$stage)
mono_de[mono_de$stage == "early", ]$color <- "yellow"
mono_de[mono_de$stage == "middle", ]$color <- "orange"
mono_de[mono_de$stage == "late", ]$color <- "red"


rownames(mono_de) <- mono_de$id
V(gD2)$color  <- mono_de[names(as.list(V(gD2))),]$color


mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/300 * params("vertex", "size")
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
png(filename = file.path(data_path,"/graph_cell.png"), width = 3300, height = 2100)
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


pdf(file = file.path(data_path,"graph_SUPPL_stages_withNames.pdf"))
plot(lou2, gD2,  
     edge.arrow.size = .0001,
     col=V(gD2)$color,
     layout=l,
     edge.color= "darkgrey",
     mark.col=adjustcolor("grey", alpha = 0.3),
     mark.border=adjustcolor("grey", alpha = 0.3),
     vertex.label.color="#333333", 
     vertex.label.font=V(gD2)$font,
     vertex.label.cex=0.5,
     vertex.label.dist=V(gD2)$dist , 
     vertex.label.degree=V(gD2)$degree,
     vertex.label.family="Helvetica",
     vertex.label = V(gD2)$label,
     vertex.shape=V(gD2)$shape,
     vertex.frame.width=0.1,
     vertex.frame.color = "black")

dev.off()



####### NETWORK COLORED BY NUMBER OF CELLS

# ---- 1. Get the cell with the highest logFC 

genes <- names(as.list(V(gD2)))

de_all$ids <- gsub("-unknown","",de_all$gene_name) 
de_all$ids <- gsub("MAMU-", "MAMU.", de_all$gene_name, fixed = T)
de_mono <- de_all[de_all$celltype == "Monocyte",]
# Calc number of celltypes
de_reduced <- unique(de_all[,c("ids", "celltype")])
de_n_celltypes <- de_reduced %>% dplyr::group_by(ids) %>% dplyr::count()
rownames(de_n_celltypes) <- de_n_celltypes$ids

# -------------------
# Change Color 
V(gD2)$color <- "grey"
de_n_celltypes$color <- "grey"
de_n_celltypes$n <- as.character(de_n_celltypes$n)
de_n_celltypes[de_n_celltypes$n == "1", ]$color <- "#ffff66"
de_n_celltypes[de_n_celltypes$n == "2", ]$color <- "#00cc99"
de_n_celltypes[de_n_celltypes$n == "3", ]$color <- "#3366ff"

de_n_celltypes <- as.data.frame(de_n_celltypes)
rownames(de_n_celltypes) <- de_n_celltypes$id
V(gD2)$color  <- de_n_celltypes[names(as.list(V(gD2))),]$color


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
png(filename = file.path(data_path,"graph_SUPPL_spec_withNames.png"), width = 3300, height = 2100)
pdf(file = file.path(data_path,"graph_SUPPL_spec_withNames.pdf"))
plot(lou2, gD2,
     edge.arrow.size = .0001,
     col=V(gD2)$color,
     layout=l,
     edge.color= "darkgrey",
     mark.col=adjustcolor("grey", alpha = 0.3),
     mark.border=adjustcolor("grey", alpha = 0.3),
     vertex.label.color="#333333", 
     vertex.label.font=V(gD2)$font,
     vertex.label.cex=0.5,
     vertex.label.dist=V(gD2)$dist , 
     vertex.label.degree=V(gD2)$degree,
     vertex.label.family="Helvetica",
     vertex.label = V(gD2)$label,
     vertex.shape=V(gD2)$shape,
     vertex.frame.width=0.1,
     vertex.frame.color = "black")

dev.off()




