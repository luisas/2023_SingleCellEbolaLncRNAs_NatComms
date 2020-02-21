

remove_unconcordant_prediction <- function(lncRNAs_out,mRNAs_out){
  # remove genes, whose transcript are predicted as both lncrnas and mrnas
  # no return, directly modify the inputs
  lo <- deparse(substitute(lncRNAs_out))
  mo <- deparse(substitute(mRNAs_out))
  
  df <- data.frame()
  df <- rbind(df, data.frame(gene_id = lncRNAs_out$gene_id, transcript_id = lncRNAs_out$transcript_id, pred = "lncrna" ))
  df <- rbind(df, data.frame(gene_id = mRNAs_out$gene_id, transcript_id = mRNAs_out$transcript_id, pred = "mrna" ))
  
  unconcordant_prediction <- df %>%  dplyr::group_by(gene_id) %>% dplyr::summarise(Unique_Elements =  dplyr::n_distinct(pred)) %>%  dplyr::filter( Unique_Elements > 1)
  
  # Remove unconcordant predictions from sets 
  lncRNAs_out <-lncRNAs_out[!(lncRNAs_out$gene_id %in% unconcordant_prediction$gene_id), ]
  mRNAs_out <- mRNAs_out[!(mRNAs_out$gene_id %in% unconcordant_prediction$gene_id), ]
  assign(lo, lncRNAs_out, envir = .GlobalEnv)
  assign(mo, mRNAs_out, envir = .GlobalEnv)
}


get_only_min_transcript_old <- function(gr){
  # given a genomic range gets only the maximum transcript of a gene
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "range_width" = width(ranges(gr)))
  gene_with_multiple_isoforms <-df[!duplicated(df$transcript_id),] %>% dplyr::group_by(gene_id) %>% dplyr::summarize(number=dplyr::n()) %>% dplyr::filter(number > 1)
  collapsed <-df %>% dplyr::group_by(gene_id,transcript_id) %>% dplyr::summarize("range" = sum(range_width)) %>% dplyr::group_by(gene_id) %>% dplyr::slice(which.min(range))
  gene_with_one_isoform <-df[!duplicated(df$transcript_id),] %>% dplyr::group_by(gene_id) %>% dplyr::summarize(number=dplyr::n()) %>% dplyr::filter(number == 1) 
  gr <- gr[gr$transcript_id %in% collapsed$transcript_id ,]
  return(gr)
}

get_only_min_transcript <- function(gr){
  # given a genomic range gets only the maximum transcript of a gene
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "range_width" = width(ranges(gr)))
  gene_with_multiple_isoforms <-df[!duplicated(df$transcript_id),] %>% dplyr::group_by(gene_id) %>% dplyr::summarize(number=dplyr::n()) %>% dplyr::filter(number > 1)
  collapsed <-df %>% dplyr::group_by(gene_id,transcript_id) %>% dplyr::summarize("range" = sum(range_width)) %>% dplyr::group_by(gene_id) %>% dplyr::slice(which.min(range))
  gene_with_one_isoform <-df[!duplicated(df$transcript_id),] %>% dplyr::group_by(gene_id) %>% dplyr::summarize(number=dplyr::n()) %>% dplyr::filter(number == 1) 
  gr <- gr[gr$transcript_id %in% collapsed$transcript_id ,]
  return(gr)
}
get_only_max_transcript <- function(gr){
  # given a genomic range gets only the maximum transcript of a gene
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "range_width" = width(ranges(gr)))
  gene_with_multiple_isoforms <-df[!duplicated(df$transcript_id),] %>% dplyr::group_by(gene_id) %>% dplyr::summarize(number=dplyr::n()) %>% dplyr::filter(number > 1)
  collapsed <-df %>% dplyr::group_by(gene_id,transcript_id) %>% dplyr::summarize("range" = sum(range_width)) %>% dplyr::group_by(gene_id) %>% dplyr::slice(which.max(range))
  gene_with_one_isoform <-df[!duplicated(df$transcript_id),] %>% dplyr::group_by(gene_id) %>% dplyr::summarize(number=dplyr::n()) %>% dplyr::filter(number == 1) 
  gr <- gr[gr$transcript_id %in% collapsed$transcript_id ,]
  return(gr)
}

get_nr_exons <- function(gr){
  df <- data.frame("gene_id" = gr$gene_id,"exon_number" = as.numeric(gr$exon_number))
  number_exons <- df %>% dplyr::group_by(gene_id) %>%dplyr::summarize(max_exon = max(exon_number))
  return(number_exons)
}
add_type <- function(mean_expression, ids, name){
  mask <- mean_expression$id %in% ids
  mean_expression[mask,]$type <- name
  return(mean_expression)
}


barplot_exon_count <- function(gr, type, col){
  ## extract the number of exons
  gr <- gr[gr$type =="exon",]
  gr <- get_only_max_transcript(gr)
  df_l <- data.frame(get_nr_exons(gr))
  df_l$type <- type
  
  h <- ggplot(df_l, aes(x=max_exon, fill=type)) + 
    geom_histogram(position="identity", binwidth =1)+
    xlim(1,15)
  
  h_plotdata <- ggplot_build(h)$data[[1]]
  h_plotdata$group <- as.factor(h_plotdata$group)
  levels(h_plotdata$group) <- c(type)
  
  ## plot with geom_bar
  p1 <-ggplot(h_plotdata, aes(x=x, y=y, fill = group, palette = col)) +
    geom_bar(stat = "identity", width = 0.8) +
    theme(legend.title=element_blank())+
    labs(y = "", x = "")+
    theme(legend.title=element_blank())+ theme(legend.position = "none")+
    scale_x_continuous( labels = as.character(h_plotdata$x), breaks = (h_plotdata$x)) + 
    theme(plot.title = element_text(hjust = 0.5))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"))+
    scale_fill_manual(values=c(col))
  return(p1)
}


transcript_length_density <- function(gr,type,color){
  # Get only exons of the reference, retain only the max transcript in term of length
  #gr <- gr[gr$type =="exon",]
  gr <- get_only_max_transcript(gr)
  gr <- gr[gr$type =="exon",]
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "range_width" = width(ranges(gr)))
  collapsed <- df%>% dplyr::group_by(gene_id,transcript_id) %>% dplyr::summarize("range" = sum(range_width))
  p1 <- ggplot(collapsed, aes(x=range)) +
    geom_density(alpha=0.5, fill = color, color = "Darkgrey")+
    labs(y = "", x = "" )+
    theme(plot.title = element_text(hjust = 0.5))+
    xlim(0,10000)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"))
  return(p1)
}



exon_length_density <- function(gr,type,color){
  gr <- gr[gr$type =="exon",]
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "exon_number" = gr$exon_number, "range_width" = width(ranges(gr)))
  collapsed <- df%>% dplyr::group_by(gene_id,transcript_id, exon_number) %>% dplyr::summarize("range" = sum(range_width))
  p1 <- ggplot(collapsed, aes(x=range)) +
        geom_density(position="identity", alpha=0.5, fill = color,color = "Darkgrey", binwidth = 100)+
        labs(x = "", y = "") +
        theme(legend.title=element_blank())+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"))+
        theme(plot.title = element_text(hjust = 0.5))+xlim(0,2000)
  return(p1)
}

calc_exon_length <- function(gr, type){
  gr <- gr[gr$type =="exon",]
  gr <- get_only_max_transcript(gr)
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "exon_number" = gr$exon_number, "range_width" = width(ranges(gr)))
  collapsed <- df%>% dplyr::group_by(gene_id,transcript_id, exon_number) %>% dplyr::summarize("range" = sum(range_width))
  collapsed["type"] <- type
  return(collapsed)
}


calc_transcript_length <- function(gr, type){
  gr <- gr[gr$type =="exon",]
  gr <- get_only_max_transcript(gr)
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "range_width" = width(ranges(gr)))
  collapsed <- df%>% dplyr::group_by(gene_id,transcript_id) %>% dplyr::summarize("range" = sum(range_width))
  collapsed["type"] <- type
  return(collapsed)
}


calc_transcript_length_min <- function(gr, type){
  gr <- gr[gr$type =="exon",]
  gr <- get_only_min_transcript(gr)
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "range_width" = width(ranges(gr)))
  collapsed <- df%>% dplyr::group_by(gene_id,transcript_id) %>% dplyr::summarize("range" = sum(range_width))
  collapsed["type"] <- type
  return(collapsed)
}

remove_one_exon <- function(gr){
  df_l <- data.frame(get_nr_exons(gr[gr$type =="exon",]))
  depleted <- gr[gr$gene_id %in% df_l[df_l$max_exon > 1,]$gene_id,]
  return(depleted)
}

plot_expression <- function(max_expression){
  plot <- ggplot(max_expression, aes(x=factor(type,level = c("lncrna", "novel", "mrna")) , y=expr, fill = type )) + 
    geom_boxplot(color="darkgrey",  alpha=1.0, fill = palette_expression )+
    labs(y = "logCPM", x = "", title = "Expression" )+
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"))
  return(plot)
}

plot_stats_annotation <- function (novel_expressed,lncRNAs_ref,lncRNAs_ref_human, mrna_ref_human, mRNAs_ref,df){
  palette <-brewer.pal(5,"Paired")
  palette <- c("#F4E3ED", "#E9A3C9", "#C51B7D", "#E6F5D0", "#A1D76A")
  levels <- c('Novel LncRNAs', 'Known lncRNAs: Macaque ', 'Known lncRNAs: Human',"Known mRNAs: Macaque", "Known mRNAs: Human" )
  # --------------------
  ## EXON COUNT
  # --------------------
  
  ec1 <- barplot_exon_count(novel_expressed, "Novel", palette[1])
  ec2 <- barplot_exon_count(lncRNAs_ref, "lncRNAs - Reference Macaque", palette[2])
  ec3 <- barplot_exon_count(lncRNAs_ref_human, "lncRNAs - Reference Human", palette[3])+labs(y = "Frequency")
  ec4 <- barplot_exon_count(mrna_ref_human, "mRNAs - Reference Human", palette[4])
  ec5 <- barplot_exon_count(mRNAs_ref, "mRNAs - Reference Macaque", palette[5])+labs(x = "Number of Exons") 
  a <- ggarrange( ec1,ec2,ec3,ec4,ec5,  ncol=1, nrow=5, top = "Transcripts Lengths") 
  # Exon count generally lower in lncrnas than mrnas (Ok. same as sources)
  
  # --------------------
  ## Transcript lengths: still weird, very long ones ...
  # --------------------
  tl1 <-transcript_length_density(novel_expressed," NOVEL",palette[1]) 
  tl2 <-transcript_length_density(lncRNAs_ref," REF MACAQUE",palette[2]) 
  tl3 <-transcript_length_density(lncRNAs_ref_human,"REF HUMAN",palette[3]) + labs( y  = "Density")
  tl4 <-transcript_length_density(mrna_ref_human," REF MACAQUE",palette[4]) 
  tl5 <-transcript_length_density(mRNAs_ref,"REF HUMAN",palette[5]) + labs( x = " Transcript length")
  b <- ggarrange( tl1,tl2,tl3,tl4,tl5,  ncol=1, nrow=5, top = "Transcripts Lengths",  legend="right")
  
  
  # --------------------
  ## Transcript lengths - Boxplot
  # --------------------
  df <- data.frame()
  df <- rbind(df,data.frame(calc_transcript_length(novel_expressed, "Novel LncRNAs")))
  df <- rbind(df,data.frame(calc_transcript_length(lncRNAs_ref, "Known lncRNAs: Macaque ")))
  df <- rbind(df,data.frame(calc_transcript_length(lncRNAs_ref_human, "Known lncRNAs: Human")))
  df <- rbind(df,data.frame(calc_transcript_length(mrna_ref_human, "Known mRNAs: Human")))
  df <- rbind(df,data.frame(calc_transcript_length(mRNAs_ref, "Known mRNAs: Macaque")))
  
  p <- ggplot(df, aes(x = factor(type, level = levels),  y = range )) +
    labs( x = "", y = "" , title = "Transcript length")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10))+
    geom_boxplot(outlier.shape=NA, fill = palette,color = "Darkgrey") + ylim(0,10000)+
    theme(legend.title=element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"))+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  # --------------------
  ## Transcript lengths min - Boxplot
  # --------------------
  df <- data.frame()
  df <- rbind(df,data.frame(calc_transcript_length_min(novel_expressed, "Novel LncRNAs")))
  df <- rbind(df,data.frame(calc_transcript_length_min(lncRNAs_ref, "Known lncRNAs: Macaque ")))
  df <- rbind(df,data.frame(calc_transcript_length_min(lncRNAs_ref_human, "Known lncRNAs: Human")))
  df <- rbind(df,data.frame(calc_transcript_length_min(mrna_ref_human, "Known mRNAs: Human")))
  df <- rbind(df,data.frame(calc_transcript_length_min(mRNAs_ref, "Known mRNAs: Macaque")))
  
  plotmin <- ggplot(df, aes(x = factor(type, level = levels),  y = range )) +
    labs( x = "", y = "" , title = "Transcript length ( min transcript length ) ")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10))+
    geom_boxplot(outlier.shape=NA, fill = palette,color = "Darkgrey") + ylim(0,10000)+
    theme(legend.title=element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"))+
    theme(plot.title = element_text(hjust = 0.5))

  
  # --------------------
  ## Exon lengths - Densities 
  # --------------------
  p1 <-exon_length_density(novel_expressed," NOVEL",palette[1]) 
  p2 <-exon_length_density(lncRNAs_ref," REF MACAQUE",palette[2]) 
  p3 <-exon_length_density(lncRNAs_ref_human,"REF HUMAN",palette[3]) + labs( y  = "Density")
  p4 <-exon_length_density(mrna_ref_human," REF MACAQUE",palette[4]) 
  p5 <-exon_length_density(mRNAs_ref,"REF HUMAN",palette[5]) + labs( x = " Exon length")
  d <- ggarrange( p1,p2,p3,p4,p5,  ncol=1, nrow=5, top = "Exon Lengths",  legend="right")
  
  # --------------------
  ## Exon lengths - Boxplot
  # --------------------
  df <- data.frame()
  df <- rbind(df,data.frame(calc_exon_length(novel_expressed, "Novel LncRNAs")))
  df <- rbind(df,data.frame(calc_exon_length(lncRNAs_ref, "Known lncRNAs: Macaque ")))
  df <- rbind(df,data.frame(calc_exon_length(lncRNAs_ref_human, "Known lncRNAs: Human")))
  df <- rbind(df,data.frame(calc_exon_length(mrna_ref_human, "Known mRNAs: Human")))
  df <- rbind(df,data.frame(calc_exon_length(mRNAs_ref, "Known mRNAs: Macaque")))
  
  p1 <- ggplot(df, aes(x = factor(type, level = levels),  y = range )) +
    labs( x = "", y = "" , title = "Exon length")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=10))+
    geom_boxplot(outlier.shape=NA, fill = palette,color = "Darkgrey") + ylim(0,1000)+
    theme(legend.title=element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"))+
    theme(plot.title = element_text(hjust = 0.5))

  return(list(a,b,p,d,p1,plotmin))
}
