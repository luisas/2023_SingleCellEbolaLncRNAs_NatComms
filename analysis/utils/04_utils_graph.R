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
  rownames(netm_filt) <- unlist(lapply(rownames(netm_filt), get_orthologname_))
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
  rownames(m) <- unlist(lapply(rownames(m), get_orthologname_))
  h1 <- Heatmap(m,  top_annotation = column_ha,right_annotation =row_ha_type,   column_order = c("baseline","early", "middle", "late"))+ggtitle(paste("Community: ",n)) 
  
  p1
  h1
  return(list(p1,h1))
}


inspect_gene_expression <- function(immune.combined, monocyte, gene){
  expression_gene <- monocyte[rownames(monocyte)== gene, ]@assays$RNA@data
  colnames(expression_gene) <- monocyte$group
  df <- reshape2::melt(as.matrix(expression_gene))
  names(df) <- c("gene", "group", "value")
  df <- df[df$value !=0,]
  # Visualize boxplot with z scores 
  table(df$group)
  p1 <- ggplot(df, aes( x = group, y = value, col = group, fill = group) )+geom_violin(alpha = 0.7, color = "dark grey ")+theme_paper+scale_fill_manual(values= c(brewer.pal(8, "Pastel2")[5:9]))+ geom_boxplot(width=0.1, color = "dark grey ")
  
  p2 <- plot_seperate_features(immune.combined, gene  = c(gene), ident = "group")
  
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
  
  
  data.plot$orth <- unlist(lapply(data.plot$features.plot, get_orthologname_))
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "orth", 
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



Dotplot_data <- function (object, assay = NULL, features, cols = c("lightgrey", 
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
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = rev(x = features))
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
  return(data.plot)
}



get_expression_summary_gene <- function(gene, mono = monocyte, ebola_genome_percentage_df_ = ebola_genome_percentage_df, zeros = F ){
  if(!(gene %in% rownames(mono) )){
    return("Gene id not found!")
  }
  expression_gene <- mono[rownames(mono)== gene, ]@assays$RNA@data
  df <- reshape2::melt(as.matrix(expression_gene))
  names(df) <- c("gene", "cell", "value")
  ebola_genome_percentage_df_$cell <- rownames(ebola_genome_percentage_df_)
  df <- merge(ebola_genome_percentage_df_, df, by="cell")
  if(zeros == F ){
    df <- df[df$value !=0,]
  }
  return(df)
}

boxplot_bystanders <- function(df, baseline = "media", live = "live"){
  # Boxplot bystander
  df$class <- NA
  df[df$cond %in% baseline, ]$class <- "Baseline" 
  df[df$cond %in% c(live) & df$classification == "Not Infected",]$class <- "Bystanders"
  df[df$cond %in% c(live) & df$classification == "Infected",]$class <- "Infected"
  df$class <- factor(df$class, levels = c("Baseline", "Bystanders", "Infected"))
  p <- ggplot(df[!is.na(df$class),], aes(x = class, y = value, fill = class))+geom_boxplot(notch = F, alpha=0.8)+theme_paper+scale_fill_manual(values = wes_palette("Zissou1", 5, type = "discrete")[c(1,4,5)])+ggtitle(get_orthologname(unique(df$gene)))
  my_comparisons <- list( c("Baseline", "Bystanders"), c("Bystanders", "Infected"), c("Baseline", "Infected"))
  p2 <- ggboxplot(df[!is.na(df$class),],x = "class", y = "value", fill = "class" )+ stat_compare_means(comparisons = my_comparisons, label = "p.signif")+theme_paper+scale_fill_manual(values = wes_palette("Zissou1", 5, type = "discrete")[c(1,4,5)])+ggtitle(get_orthologname(unique(df$gene)))
  return(list(p,p2))
}


get_o <- function(string){
  d <- orthologs[orthologs$lnc %in% unique(ref$transcript_id[(string == ref$gene_name)]), ]
  if(nrow(d)> 0){ d$gene <- string }
  return(d)
}

# ---------------------
#   Proportions test 
# ---------------------
get_prop <- function(gene_id){
  
  # How many infected and not infected cells in the dataset for real 
  n_all <- sum(table(mono_live$infection))
  p_infected <- table(mono_live$infection)["Infected"]/n_all
  p_bystander <- table(mono_live$infection)["Not Infected"]/n_all
  
  # Chi Squared test 
  df <- get_expression_summary_gene(gene_id, mono =mono_live)
  obs <- c(table(df$classification)["Infected"], table(df$classification)["Not Infected"])
  exp <- c(p_infected, p_bystander)
  test <- chisq.test(obs, p=exp)
  
  df <- Dotplot_data(mono, features = c(gene_id), group.by = "cond", split.by = "infection", dot.scale =10,scale.by = "radius")
  df$cond <- unlist(lapply(df$id, function(x) str_split(x, "_")[[1]][1]))
  df$inf <- unlist(lapply(df$id, function(x) str_split(x, "_")[[1]][2]))
  df$cond <- factor(df$cond, levels = c("media", "irrad", "live"))
  df$inf <- factor(df$inf, levels = c("Not Infected", "Infected"))
  
  p <- ggplot(df[df$cond%in% c("live", "media"),], aes(x = cond, y =pct.exp, fill = inf, col = inf))+geom_bar(stat="Identity", position="dodge", alpha = 0.7)+theme_paper+xlab("")+ylab("%Expressed")+scale_color_manual(values = c("#E5D47F", "#C10000"))+ggtitle(paste(get_orthologname(gene_id),",  Pval:",round(test$p.value,10), sep = " "))+scale_fill_manual(values = c("#E5D47F", "#C10000"))
  return(p)
}
get_orthologname <- function(string){
  d <- orthologs[orthologs$lnc %in% unique(ref$transcript_id[(string == ref$gene_name)]), ]
  if(nrow(d)==0){
    return(gsub("-unknown", "", string) )
  }
  return(unique(as.character(d$orthologGeneSymbol)))
}

heatmap_celltypes <- function(plot.data, perc_heatmap = F){
  m <- acast(plot.data, features.plot~id, value.var="avg.exp.scaled")
  m_pc<- acast(plot.data, features.plot~id, value.var="pct.exp")
  col_fun = colorRamp2(c(min(m), max(m)), c("lightblue", "dark red"))
  rownames(m) <- unlist(lapply(rownames(m), get_orthologname_))
  rownames(m_pc) <- unlist(lapply(rownames(m_pc), get_orthologname_))
  
  p1 <- Heatmap(m, name = "Average Expression", col = col_fun, cluster_rows = T, cluster_columns = F,rect_gp = gpar(col = "white", lwd = 1),column_split = unlist(lapply(colnames(m), function(x) str_split(x, "_")[[1]][1])),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  #grid.text(sprintf("%.1f", m_pc[i, j]), x, y, gp = gpar(fontsize = 10))
                  grid.circle(x = x, y = y, r = abs( m_pc[i, j])/200 * min(unit.c(width, height)), gp = gpar(fill = NA))
                })
  
  col_fun = colorRamp2(c(min(m_pc), max(m_pc)), c("lightblue", "dark red"))
  p2 <- Heatmap(m_pc, name = "Average Expression", col = col_fun, cluster_rows = F, cluster_columns = F,rect_gp = gpar(col = "white", lwd = 1),column_split = unlist(lapply(colnames(m_pc), function(x) str_split(x, "_")[[1]][1])),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  #grid.text(sprintf("%.1f", m_pc[i, j]), x, y, gp = gpar(fontsize = 10))
                  grid.circle(x = x, y = y, r = abs( m_pc[i, j])/200 * min(unit.c(width, height)), gp = gpar(fill = col_fun(m_pc[i, j])))
                })
  if(perc_heatmap == T ){
    return(list(p1, p2))
  }else{
    return(p1)
  }
  
}
get_gene_id <- function(gene){
  gene_id <- unique(ref[!is.na(ref$transcript_id) & ref$transcript_id %in% unique(orthologs[orthologs$orthologGeneSymbol == gene,]$lnc),]$gene_name)
  gene_id  <- rownames(subset(immune.combined, features = c(gene_id)))
  return(gene_id)
}
plot_region <- function(gene, range = 100000, palette_plot_percentage = palette_plot_percentage, genes = "NULL",getmax = T ){
  gene_id <- unique(ref[ref$gene_name == gene,]$gene_id)
  chr <- paste("chr",unique(as.character(seqnames(ref[ref$gene_id == gene_id,]))), sep = "")
  
  # Define region nearby 
  from  <- max(0,(min(start(ref[ref$gene_id == gene_id,])) -range))
  print(gene)
  to <- (max(end(ref[ref$gene_id == gene_id,]))+range)
  ref_small<- ref[start(ref) > from  & end(ref) < to   ,]
  if(getmax == T){
    ref_small <- get_only_max_transcript(ref_small)
  }
  if(genes  != "NULL"){
    ref_small <- ref_small[ref_small$gene_name %in% genes,]
    print("filter")
  }
  print(genes)
  # Change reference levels to UCSC
  seqlevelsStyle(ref_small) <- "UCSC"
  
  geneModels <- as.data.frame(ref_small)
  if(!(chr %in% as.character(geneModels$seqnames))){
    return(NA)
  }
  geneModels <- geneModels[as.character(geneModels$seqnames) == chr, ]
  geneModels$gene <- geneModels$gene
  geneModels <- geneModels[geneModels$type =="exon",]
  geneModels$symbol <- gsub("-unknown", "", geneModels$gene_name)
  geneModels$exon <- geneModels$exon_id
  geneModels$transcript <- geneModels$transcript_id
  geneModels$gene <- geneModels$gene_id
  
  direction <- ifelse(geneModels$strand == "+", ">", "<")
  geneModels$symbol <- paste(geneModels$symbol, direction)
  if("StringTie" %in% geneModels$source){
    geneModels[geneModels$source == "StringTie", ]$gene_biotype <- "lncRNA"
  }
  geneModels <- geneModels[geneModels$gene_biotype %in% c("protein_coding", "lncRNA"),]
  
  #geneModels_14 <- geneModels_14[geneModels_14$symbol %in%  c("ENSMMUG00000061948_unkown", "THY1"),]
  gtrack1 <- GeneRegionTrack(geneModels, genome  = "rheMac10", chr = chr, name = "Macaque rheMac10")
  atrack <- GeneRegionTrack(geneModels, genome  = "rheMac10", chr = chr, name = "Macaque rheMac10", cex.axis = 15)
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = "rheMac10", chromosome = chr)
  
  geneModels[geneModels$gene_id == gene_id,]$gene_biotype <- "lncRNA_highlight"
  #geneModels[geneModels$gene_name == "THY1",]$gene_biotype <- "mRNA_highlight"
  
  feature(atrack) <- geneModels$gene_biotype
  p <- plotTracks(list(itrack,gtrack, atrack), showId = TRUE, from = from , to = to, col = NULL, 
                  background.title = "white", col.title = "black", cex.title = 1, lncRNA = palette_plot_percentage[1], protein_coding = palette_plot_percentage[2],
                  lncRNA_highlight = palette_plot_percentage[3],fontsize=8,cex.group=1,
                  mRNA_highlight = palette_plot_percentage[4],cex.id = 2, sizes = c(0.3,0.3,1), cex.main = 1)
  return(p)
}


Plot_dotplot2 <- function(genes, mono, group= "condition", split = NA){
  if(!is.na(split)){
    plot.data <- Dotplot_data(mono, features = c(genes), group.by = group, split.by = split, scale.by = "radius", cols = c( "lightblue", "darkred", "grey", "grey"))

  }else{
    plot.data <- Dotplot_data(mono, features = c(genes), group.by = group, scale.by = "radius", cols = c( "lightblue", "darkred"))

  }
  plot.data$name <- unlist(lapply(plot.data$features.plot, get_orthologname))
  plot.data$labels <- as.character(plot.data$id)
  plot.data$labels  <- gsub("media Not Infected", "Bystander", plot.data$labels )
  plot.data$labels  <- gsub("live Not Infected", "Live Bystander", plot.data$labels )
  scale.func <- switch(EXPR = "radius", size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  p <- ggplot(plot.data, aes(x = id,y = name, size= pct.exp, col = avg.exp.scaled))+geom_point()+theme_paper+xlab("")+xlab("")+scale_colour_gradient(low = "lightblue", high = "dark red", na.value = NA)+theme(axis.text.x = element_text(angle = 20, vjust = 0.9, hjust=1))+ scale.func(range = c(0, 15), limits = c(0, max(plot.data$pct.exp)))+ylab("")
  
  p2 <- heatmap_celltypes(plot.data )
  return(list(p,p2))
}



 plot_boxplot<- function(gene_id, mono_live_h24_inf, ebola_genome_percentage_df_ = ebola_genome_percentage_df){
  df <- get_expression_summary_gene(gene_id, mono =mono_live_h24_inf, ebola_genome_percentage_df_)
  df_2 <- get_expression_summary_gene(gene_id, mono =mono_live_h24_notinf, ebola_genome_percentage_df_)
  p1 <- ggplot(df, aes( x = classification, y = value, col = classification, fill = classification))+geom_boxplot(alpha = 0.1)+theme_paper+scale_y_log10()+ylab("logCP10K")+xlab("")+theme(axis.text.x = element_blank())+theme(legend.position = "")+theme(axis.ticks.x.bottom = element_blank(), plot.margin = unit(c(2.33,-0.9,0.95,0), "cm"))
  # Visualize viral load and gene expression 
  p2 <- ggMarginal(ggplot(df, aes(x = percentage_viral_reads, y = value, col = classification))+geom_point(alpha = 0.5)+theme_minimal()+theme_paper+theme(axis.text.y = element_blank(), legend.position = "")+xlab("Viral load")+ylab("")+ggtitle(get_orthologname_(unique(df$gene)))+geom_line(stat = "summary_bin", binwidth = 0.1)+geom_smooth()+scale_x_log10()+scale_y_log10()+theme(plot.margin = unit(c(0,0,0,0), "cm")),type = "histogram", bins = 40)
  ggarrange(p1, p2, widths = c(0.15,2))
}


check_expressed_vs_infection <- function(gene_id, mono = mono, ebola_genome_percentage_df = ebola_genome_percentage_df){
  # When neat is expressed check infected not infected, status 
  gene_expressed <- as.data.frame(t(mono[c(gene_id),]@assays$RNA@data))
  colnames(gene_expressed) <- c("gene")
  gene_expressed$exp <- "Not Expressed"
  gene_expressed[gene_expressed$gene > 1, ]$exp <- "Expressed"
  gene_expressed$infection <- ebola_genome_percentage_df[rownames(gene_expressed),]$classification
  gene_expressed$infection <- factor(gene_expressed$infection, levels =c("Not Infected", "Infected"))
  
  gene_expressed$status <- dp_dn[rownames(gene_expressed),]$status
  
  
  summary_exp <- gene_expressed %>% dplyr::group_by(exp, infection) %>%dplyr::summarise(n = n())%>% dplyr::mutate(freq = n *100 / sum(n))
  p1 <- ggplot(summary_exp, aes( x = exp, fill = infection, y = freq ))+geom_bar(stat = "identity", position = "dodge")+theme_paper+scale_fill_manual(values = c("#F2D994", "#A93E00"))+ylab("%Cells")+xlab("")+ggtitle(get_orthologname_(gene_id))+ylim(c(0,100))
  
  summary_exp <- gene_expressed %>% dplyr::group_by(exp, status) %>%dplyr::summarise(n = n())%>% dplyr::mutate(freq = n *100 / sum(n))
  p2 <- ggplot(summary_exp, aes( x = exp, fill = status, y = freq ))+geom_bar(stat = "identity", position = "dodge", alpha = 0.9)+theme_paper+scale_fill_manual(values = wes_palette("Darjeeling1",4))+ylab("%Cells")+xlab("")+ggtitle(get_orthologname_(gene_id))
  return(list(p1,p2))
}



