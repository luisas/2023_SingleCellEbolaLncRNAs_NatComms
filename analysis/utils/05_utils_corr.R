
get_orthologname_ <- function(string, orthologs_ = orthologs){
  d <-orthologs_[gsub("-unknown", "", orthologs_$gene_id) == sub("-unknown", "",string), ]
  if(nrow(d)==0){
    return(gsub("-unknown", "", string) )
  }else  if(nrow(d)>1){
    print("----------")
    print(d)
  }else{
    return(unique(as.character(d$orthologGeneSymbol)))
  }
  
}
# Methods
filter_significant_correlations <- function(correlations_exvivo, ref, all_lncrnas, annotated_mrnas, FDRTHRESHOLD, RHOTHRESHOLD, usep = F){
  if(usep == T){
    correlations_exvivo <- correlations_exvivo[!is.na(correlations_exvivo$pval) & correlations_exvivo$pval < FDRTHRESHOLD & abs(correlations_exvivo$rho) > RHOTHRESHOLD, ]
  }else{
    correlations_exvivo <- correlations_exvivo[!is.na(correlations_exvivo$qval) & correlations_exvivo$qval < FDRTHRESHOLD & abs(correlations_exvivo$rho) > RHOTHRESHOLD, ]
  }
  # Add Gene Type
  correlations_exvivo$type <- "-"
  correlations_exvivo[correlations_exvivo$gene %in% all_lncrnas ,]$type <- "lnc"
  correlations_exvivo[correlations_exvivo$gene %in% annotated_mrnas ,]$type <- "pc"
  # Extract lnc only
  sig_cor_lnc_exvivo <- correlations_exvivo[correlations_exvivo$type == "lnc",]
  sig_cor_lnc_exvivo$gene_name <- unlist(lapply(sig_cor_lnc_exvivo$gene, get_orthologname_))
  # Remove ebola gene 
  sig_cor_lnc_exvivo <- sig_cor_lnc_exvivo[!(sig_cor_lnc_exvivo$gene %in% unique(ref[seqnames(ref) == "EBOV_Kikwit",]$gene_name)),]
  return(list(correlations_exvivo, sig_cor_lnc_exvivo))
}
PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}
# Compute mean expression and % cells 
get_average <- function(gene, mono_live_h24_inf){
  exp <- mono_live_h24_inf[gene,]@assays$RNA@data
  v <- as.vector(exp)
  b <- data.frame(t(Rmisc::CI(v[v>1])), gene = gene)
  b$pct.exp <- PercentAbove(v,1)*100
  b$n_cells <- count(v > 1)
  return(b)
}


get_expression_average <- function(expression_matrix, viral_load,gene,  window = 10, dropzeros = F){
  # 1. Create a DF with gene expression and viral load 
  df <- data.frame(exp = expression_matrix[gene, ], viral_load)
  # 2. Remove zeros
  if(dropzeros == T){
    df <- df[df$exp != 0,]
  }
  df_ordered <- df[order(df$viral_load),]
  if(nrow(df_ordered) > 30 ){
    n <- window
    df_ordered$bin <- cut(1:length(df_ordered$viral_load), breaks =c( seq(0, length(df_ordered$viral_load), n), ceiling(length(df_ordered$viral_load)/n)*n) )
    avg_bins <- df_ordered %>% dplyr::group_by(bin) %>% dplyr::summarise(avg = mean(exp))
    avg_bins_vl <- df_ordered %>% dplyr::group_by(bin) %>% dplyr::summarise(avg = mean(viral_load))
    df_ordered$avg <- NULL
    df_ordered$avg <- unlist(lapply(df_ordered$bin, function(x) avg_bins[x== avg_bins$bin,]$avg))
    df_ordered$avg_vl <- unlist(lapply(df_ordered$bin, function(x) avg_bins_vl[x== avg_bins_vl$bin,]$avg))
    df_ordered$gene <- gene
    df_ordered$ortholog <- get_orthologname_(gene)
    
    # Calculate correlation value
    return(df_ordered)
  }
}



# A: --------- Plot on all 
visualize_correlation <- function(genes,mono_live_h24_inf_exvivo,ebola_genome_percentage_df_exvivo, keepzeros = F){
  genes_expression <-  Reduce("rbind", lapply(genes, function(gene) get_expression_summary_gene(gene, mono_live_h24_inf_exvivo, ebola_genome_percentage_df_exvivo, keepzeros)))
  genes_expression$gene <- as.character(genes_expression$gene)
  genes_expression$orth <- unlist(lapply(as.character(genes_expression$gene), get_orthologname_))
  ggplot(genes_expression, aes(x = percentage_viral_reads, y = value, col = orth))+stat_smooth(method = "loess", formula = y ~ x,  se = F, n = 15, span =1)+
    theme_bw()+theme( legend.title = element_blank(),plot.title = element_text(size = 12), text = element_text(size = 17))+ 
    xlab("viral load (log10)")+ylab("log(CP10K+1)")+scale_x_log10()
}
# C: --------- Plot on all window 
visualize_correlation_window <- function(genes, mono_live_h24_inf_exvivo, keepzeros = F){
  expression_matrix_exvivo <- mono_live_h24_inf_exvivo@assays$RNA@data
  expressions <- Reduce(rbind, lapply(genes, function(gene) get_expression_average(expression_matrix_exvivo, log10(mono_live_h24_inf_exvivo$viral_load+1), gene,dropzeros = !keepzeros )))
  ggplot(expressions, aes(x = viral_load, y = avg, col = ortholog))+stat_smooth(method = "loess", formula = y ~ x,  se = F, n = 10, span =0.6)+
    theme_bw()+theme( legend.title = element_blank(),plot.title = element_text(size = 12), text = element_text(size = 17))+ 
    xlab("viral load (log10)")+ylab("log(CP10K+1)")+scale_x_log10()
}


plot_corr_genes_stat <- function(Approach_C_significant_exvivo){
  Approach_C_significant_exvivo$direction <- ifelse(Approach_C_significant_exvivo$rho > 0 , "positive", "negative")
  df <- data.frame(table(Approach_C_significant_exvivo[,c("type", "direction")]))
  df <- df[df$type != "-",]
  p <- ggplot(df, aes( x = type , fill = direction, y = Freq, label = Freq))+geom_bar(stat = "identity", alpha = 0.7)+geom_text(size = 4, position = position_stack(vjust = 1.1))+scale_fill_manual(values = c("dark grey", "dark red"))+theme_paper+xlab("")+ylab("# genes")
  p
}

overlapplot <- function (a,b){
  x <- list(A = a, B = b)
  ggvenn(
    x, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 0
  )
}
test_broad_genes <- function(correlations){
  test_genes <- c("STAT1", "MX1", "MX2", "HSPA5", "DYNLL1")
  correlations[correlations$gene %in% test_genes, ]
  return(correlations[correlations$gene %in% test_genes, ])
}


get_expression_subset <- function(s, label, gene){
  s <- subset(s, features = gene)
  df <- data.frame(t(s@assays$RNA@data))
  df$stage <- label
  df$viral_load <- s$viral_load
  df$q <- cut(df$viral_load, breaks = 4)
  df$q <- as.numeric(df$q)
  df$q <- factor(as.character(df$q), levels = c("1","2", "3", "4"))
  colnames(df) <- c("exp", "stage", "viral_load","bin")
  return(df)
}


plot_boxplot_cells <- function(gene,baseline, bystanders, infected, sub = ""){
  df_plotready <- rbind(get_expression_subset(baseline,"baseline",gene), get_expression_subset(bystanders,"bystander",gene), get_expression_subset(infected,"infected",gene)) 
  df_plotready <- df_plotready[df_plotready$exp != 0,]
  df_plotready$p <- ifelse(df_plotready$stage == "infected", df_plotready$bin
                           ,df_plotready$stage)
  df_plotready$p  <- factor(df_plotready$p , levels = c("baseline", "bystander", "1", "2", "3", "4"))
  
  p1 <- ggboxplot(df_plotready, x = "stage", y = "exp", col = "black", fill = "stage",palette = c(pal_jco()(2),rep("red",4)))+theme_linedraw()+theme_paper+ggtitle(gene)+xlab("")+ylab("log(CP10K+1)")+theme(legend.position = "")+labs(subtitle = sub)
  p2 <- ggboxplot(df_plotready, x = "p", y = "exp", col = "black", fill = "p",palette = c(pal_jco()(2),rep("red",4)))+theme_linedraw()+theme_paper+ggtitle(gene)+xlab("")+ylab("log(CP10K+1)")+theme(legend.position = "")+labs(subtitle = sub)
  return(list(p1,p2))
}
