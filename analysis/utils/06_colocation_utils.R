
# Per celltype, obtain all DE cis pairs 
get_cis_celltype <- function(cis_distances, de_all_genes, ct){
  de_ct <-  de_all_genes[de_all_genes$celltype == ct,]
  df <- cis_distances[cis_distances$lnc_name %in% de_ct$primerid  & cis_distances$gene_name %in% de_ct$primerid, ]
  df$celltype <- ct
  return(df)
}


get_orthologname_ <- function(string, orthologs_ = orthologs){
  d <-orthologs_[gsub("-unknown", "", orthologs_$gene_id) == sub("-unknown", "",string), ]
  if(nrow(d)==0){
    return(gsub("-unknown", "", string) )
  }else{
    return(unique(as.character(d$orthologGeneSymbol))[1])
  }
}


get_correlation_and_colocation_df <- function(cis_pairs, ident){
  # Calculate the correlation for the close ones
  df_reduced <- cis_pairs[,c("lnc","gene_name")]
  df_reduced$lnc <- paste(df_reduced$lnc, "-unknown", sep = "")
  pairs <- as.list(data.frame(t(df_reduced), stringsAsFactors = F))
  length(pairs)
  
  mono <- subset(immune.combined, ident = ident)
  expressionmatrix <- as.matrix(mono@assays$RNA@counts)
  correlations <- lapply(pairs, function(genes) {calc_correlation_genes(expressionmatrix[genes[1], ], expressionmatrix[genes[2],], genes[1], genes[2])})
  
  correlations_df <- Reduce("rbind",correlations[!is.na(correlations)])
  correlations_df$lnc <- gsub("-unknown","", correlations_df$g1)
  correlations_df$pair <- paste(correlations_df$lnc, correlations_df$g2, sep = "-")
  rownames(correlations_df) <- correlations_df$pair
  
  cis_pairs$pair <- paste(cis_pairs$lnc, cis_pairs$gene_name, sep ="-")
  cis_pairs_corr <- merge(cis_pairs, correlations_df, by = "pair", all = T)
  cis_pairs_corr <- cis_pairs_corr[,1:10]
  colnames(cis_pairs_corr) <- c("pair", "lnc_id", "pc_id", "distance", "pc_name", "lnc_name", "lnc_orth","celltype", "pval","rho")
  
  return(cis_pairs_corr)
}
# ---
calc_correlation_genes <- function(gene1_vector,gene2_vector,gene1, gene2, type = "spearman"){
  mask <- gene1_vector !=0 & gene2_vector !=0
  # calculate the correlation coefficient
  if(sum(mask)>10){
    #pc <- bayes.cor.test(c1, c2)
    pc <- cor.test(gene1_vector[mask], gene2_vector[mask], method = c(type))
    df <- data.frame( pval=pc$p.value, rho=pc$estimate)
    df$g1 <- gene1
    df$g2 <- gene2
    return(df)
  }else{
    return(NA)
  }
}

expand.grid.unique <- function(x, y, include.equals=FALSE){
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}



get_cis <- function(gr_gene, gr_hits =mrnas_ranges, maxgap =1000000L ){
  p <- findOverlapPairs(gr_gene, get_gene_only(gr_hits), maxgap=maxgap, ignore.strand = TRUE)
  print(gr_gene$gene_name)
  # Compute thedistances and select the first 5
  if(length(p@second) == 0 ){
    df <- data.frame(gene=gr_gene$gene_id, distances=NA, second = NA, chr = seqnames(gr_gene), stringsAsFactors = F )
  }else{
    distances <- unlist(lapply(1:length(p@second), function(index) GenomicRanges::distance(gr_gene, p@second[index], ignore.strand = TRUE) ))
    df <- data.frame(gene=gr_gene$gene_id, distances=distances, second = p@second$gene_id, stringsAsFactors = F, chr = seqnames(gr_gene) )
  }
  return(df)
}


get_orthologname_ <- function(string, orthologs_ = orthologs){
  d <-orthologs_[gsub("-unknown", "", orthologs_$gene_id) == sub("-unknown", "",string), ]
  if(nrow(d)==0){
    return(gsub("-unknown", "", string) )
  }else{
    return(unique(as.character(d$orthologGeneSymbol))[1])
  }
}