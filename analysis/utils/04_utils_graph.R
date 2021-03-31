set.seed(123)

# Get the graph and extract the clusters interesting for lncRNAS
prepare_hubs <- function(weightMat_myeloids, de_lnc_names = de_lnc_names){
  # Extract lncRNAS
  linkList <- getLinkList(weightMat_myeloids)
  max_n <- nrow(linkList)*0.5/100
  linkList <- getLinkList(weightMat_myeloids, reportMax = max_n)
  
  
  lnc_graphs <- linkList[as.character(linkList$regulatoryGene) %in% de_lnc_names |as.character(linkList$targetGene) %in% de_lnc_names, ]
  Gsi <- graph.data.frame(lnc_graphs,directed = F)
  
  hubs <- decompose.graph(Gsi, min.vertices = 5) 
  print(paste("#Hubs: ",length(hubs)))
  
  # Check how big they are
  print("Size of hubs: ")
  print(unlist(lapply(hubs, function(subgraph) length(V(subgraph)))))
  return(hubs)
}

# Visualize the hub and its connections

check_hub <- function(subgraph, universe_local, de_lnc_names = de_myeloids){
  # Visualize subgraph 
  V(subgraph)$color <- ifelse(names(as.list(V(subgraph))) %in% de_lnc_names, "red", "grey")
  V(subgraph)[startsWith(names(as.list(V(subgraph))), "MSTRG")]$color <- "purple"
  V(subgraph)$frame.color <- ifelse( names(as.list(V(subgraph))) %in% de_lnc_names & names(as.list(V(subgraph))) %in% orthologs$gene_id , "blue", "dark grey")
  
  plot(subgraph,  vertex.size = 8, vertex.frame.width =19 , vertex.label=NA, vertex.ltype=3)
  # Arrow size based on weigth
  #E(subgraph)$width <- E(subgraph)$weight/6
  #edge.start <- ends(subgraph, es=E(subgraph))[,1]
  #edge.col <- V(subgraph)$color[edge.start]
  #plot(subgraph, edge.color=edge.col, edge.curved=.1, vertex.label=NA, vertex.size=8)  
  
  #plot(subgraph, vertex.size=5, vertex.label.font=1, vertex.label.cex=.6)
  # Calculate enrichment
  ego <- clusterProfiler::enrichGO(gene = names(as.list(V(subgraph))),OrgDb =org.Mmu.eg.db, keyType = "SYMBOL",ont = "BP",universe =universe_local)
  p0 <- (clusterProfiler::dotplot(ego))
  p1 <- emapplot(ego)
  p3 <- cnetplot(ego)
  netm <- as_adjacency_matrix(subgraph, attr="weight", sparse=F)
  palf <- colorRampPalette(c("light grey", "red")) 
  netm_filt <- netm[rownames(netm)[rownames(netm) %in% gsub("-unknown", "",de_lnc_names)], setdiff(rownames(netm), rownames(netm)[rownames(netm) %in% gsub("-unknown", "",de_lnc_names)]), drop = FALSE]
  if(nrow(netm_filt)>1){
    p2 <- heatmap.2(netm_filt, Rowv = TRUE, Colv = TRUE, col = palf(100), 
                    scale="none", margins=c(10,10),trace='none',dendrogram='none',cexRow = 0.9)
  }else{
    p2 <- heatmap.2(rbind(netm_filt, netm_filt), Rowv = TRUE, Colv = TRUE, col = palf(100), 
                    scale="none", margins=c(10,10),trace='none',dendrogram='none',cexRow = 0.9)
  }
  
  return(list(p0,p1,p2,p3))
}


plot_graph_heatmap <- function(subgraph){
  netm <- as_adjacency_matrix(subgraph, attr="weight", sparse=F)
  palf <- colorRampPalette(c("light grey", "red")) 
  netm_filt <- netm[rownames(netm)[rownames(netm) %in% gsub("-unknown", "",de_lnc_names)], setdiff(rownames(netm), rownames(netm)[rownames(netm) %in% gsub("-unknown", "",de_lnc_names)]), drop = FALSE]
  if(nrow(netm_filt)>1){
    p2 <- heatmap.2(netm_filt, Rowv = TRUE, Colv = TRUE, col = palf(100), 
                    scale="none", margins=c(10,10),trace='none',dendrogram='none',cexRow = 0.9)
  }else{
    p2 <- heatmap.2(rbind(netm_filt, netm_filt), Rowv = TRUE, Colv = TRUE, col = palf(100), 
                    scale="none", margins=c(10,10),trace='none',dendrogram='none',cexRow = 0.9)
  }
  return(p2)
}

calc_z <- function(expression, na.remove = F){
  mean <- rowMeans(expression, na.rm = na.remove); sd <- apply(expression, 1 , function(x) sd(x, na.rm = na.remove)); 
  z <- (expression-mean)/sd
  return(z)
}
visualize_community_expression <- function(seuobj, gD2, n){
  sub <- V(gD2)[V(gD2)$community == n]
  genes <- names(V(gD2)[V(gD2)$community == n])
  print(n)

  
  get_z_scores <- function(seuratobject){
    # First computes z score per gene per cell and then calculate the mean of z score per celltype
    averaged_z_matrix <- sapply(unique(colnames(zmatrix)), function(x) rowMeans(zmatrix[,colnames(zmatrix) ==x]))
    return(averaged_z_matrix)
  }
  expression_genes <- seuobj[gsub("-unknown", "", rownames(seuobj)) %in% genes, ]@assays$RNA@data
  colnames(expression_genes) <- seuobj$group
  
  zmatrix <- as.matrix(calc_z(expression_genes, T))
  m<- sapply(unique(colnames(zmatrix)), function(x) rowMeans(zmatrix[,colnames(zmatrix) ==x]))
  rownames(m) <- gsub("-unknown", "", rownames(m))
  
  
  # Prepare dataframe for plotting
  df <- setNames(reshape2::melt(m), c('rows', 'vars', 'values'))
  names(df) <- c("gene", "group", "value")
  # Visualize boxplot with z scores 
  p1 <- ggplot(df, aes( x = group, y = value, col = group, fill = group) )+geom_boxplot(alpha = 0.7, color = "dark grey ")+theme_paper+ylim(-1,1)+scale_fill_manual(values= c(brewer.pal(8, "Pastel2")[5:8]))+ggtitle(paste("Community: ",n)) 
  # Visualize heatmap
  column_ha <- HeatmapAnnotation(comparison = anno_text( unlist(lapply(colnames(m), function(x) toupper(substr(x,1,1)))), rot = 0, gp = gpar(fontsize = 25, fill =  brewer.pal(8, "Pastel2")[5:8]), location = 0.5, just = "center"))
  
  
  df_col<- data.frame(id = names(sub), col = sub$color, stringsAsFactors = F) 
  rownames(df_col) <- df_col$id
  row_ha_type <- rowAnnotation(genetype =df_col[rownames(m),]$col, show_legend = FALSE,
                               show_annotation_name = FALSE,
                               col = list(genetype = c("grey" = "grey", "red" = "red",
                                                       "purple" = "purple"
                               )))
  
  h1 <- Heatmap(m,  top_annotation = column_ha,right_annotation =row_ha_type,   column_order = c("baseline","early", "middle", "late"))+ggtitle(paste("Community: ",n)) 
  
  p1
  h1
  return(list(p1,h1))
}


inspect_gene_expression <- function(immune.combined, monocyte, gene){
  expression_gene <- monocyte[gsub("-unknown", "", rownames(monocyte))== gene, ]@assays$RNA@data
  colnames(expression_gene) <- monocyte$group
  df <- reshape2::melt(as.matrix(expression_gene))
  names(df) <- c("gene", "group", "value")
  df <- df[df$value !=0,]
  # Visualize boxplot with z scores 
  table(df$group)
  p1 <- ggplot(df, aes( x = group, y = value, col = group, fill = group) )+geom_violin(alpha = 0.7, color = "dark grey ")+theme_paper+scale_fill_manual(values= c(brewer.pal(8, "Pastel2")[5:9]))+ geom_boxplot(width=0.1, color = "dark grey ")
  
  p2 <- plot_seperate_features(immune.combined, gene  = c(paste(gene, "-unknown", sep = "")), ident = "group")
  
  zmatrix <- as.matrix(calc_z(expression_gene, T))
  df <- reshape2::melt(as.matrix(zmatrix))
  names(df) <- c("gene", "group", "value")
  df <- df[df$value !=0,]
  p3 <- ggplot(df, aes( x = group, y = value, col = group, fill = group) )+geom_boxplot(alpha = 0.7, color = "dark grey ")+theme_paper+ geom_boxplot(width=0.1, color = "dark grey")
  
  return(list(p1, p2,p3))
}



get_df_z <- function(gene, mono, z = F ){
  expression_gene <- mono[rownames(mono) == gene, ]@assays$RNA@data
  colnames(expression_gene) <- mono$dpi
  dpis <- mono$dpi
  zmatrix <- as.matrix(calc_z(expression_gene, T))
  df <- reshape2::melt(as.matrix(zmatrix))
  names(df) <- c("gene", "group", "value")
  df_real <- reshape2::melt(as.matrix(expression_gene))
  names(df_real) <- c("gene", "group", "value")
  df_real$gr <- mono$group
  df_real <- df_real[df_real$value > 0,]
  if(z == T ){
    return(df)
  }else{
    return(df_real)
  }
  
}
calc_correlation_genes <- function(gene1_vector,gene2_vector,gene1, gene2, type = "spearman"){
  #mask <- gene1_vector !=0 & gene2_vector !=0
  # calculate the correlation coefficient
  if(T){ 
    #pc <- bayes.cor.test(c1, c2) 
    #pc <- cor.test(gene1_vector[mask], gene2_vector[mask], method = c(type))
    pc <- cor.test(gene1_vector, gene2_vector, method = c(type))
    df <- data.frame( pval=pc$p.value, rho=pc$estimate)
    df$g1 <- gene1
    df$g2 <- gene2
    return(df)
  }else{ 
    return(NA)
  }
}

#calc_correlation_genes(expressionmatrix[g1, ], expressionmatrix[g2,], g1, g2)

########################################
# Modify Seurat dotplot funciton 
########################################

PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

Dotplot_mod <- function (object, assay = NULL, features, cols = c("lightgrey", 
                                                                  "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                         group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA, 
                         scale.max = NA) 
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  data.features <- FetchData(object = object, vars = features)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  }
  else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min, 
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- gsub("-unknown", "", data.plot$features.plot)
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = rev(x = gsub("-unknown", "",features)))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
                                      split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
                         2)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
                     no = "colors")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  
  
  
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + theme_cowplot()
  if (!is.null(x = split.by)) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (is.null(x = split.by)) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}
























