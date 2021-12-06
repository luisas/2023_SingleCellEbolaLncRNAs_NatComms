
theme_paper <- theme(legend.title = element_blank())+theme(panel.background = element_rect(fill = "white", colour = "white"))+theme(panel.background = element_rect(fill = "white", colour = "grey50"))+theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 18))

# Input: GenomicRanges object
# Output: GenomicRanges object containing only the maximum transcript per gene
# Given a genomic range gets only the maximum transcript of a gene
# Computed based on the sum of the lenghts of its exons.
get_only_max_transcript <- function(gr){
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "range_width" = width(ranges(gr)))
  gene_with_multiple_isoforms <-df[!duplicated(df$transcript_id),] %>% dplyr::group_by(gene_id) %>% dplyr::summarize(number=dplyr::n()) %>% dplyr::filter(number > 1)
  collapsed <-df %>% dplyr::group_by(gene_id,transcript_id) %>% dplyr::summarize("range" = sum(range_width)) %>% dplyr::group_by(gene_id) %>% dplyr::slice(which.max(range))
  gene_with_one_isoform <-df[!duplicated(df$transcript_id),] %>% dplyr::group_by(gene_id) %>% dplyr::summarize(number=dplyr::n()) %>% dplyr::filter(number == 1) 
  gr <- gr[gr$transcript_id %in% collapsed$transcript_id ,]
  return(gr)
}


# Input: GenomicRanges object
# Outout: Number of exons per gene (Data frame)
# Def: Obtain number of exons per gene
get_nr_exons <- function(gr){
  df <- data.frame("gene_id" = gr$gene_id,"exon_number" = as.numeric(gr$exon_number))
  number_exons <- df %>% dplyr::group_by(gene_id) %>%dplyr::summarize(max_exon = max(exon_number))
  return(number_exons)
}

# Util, add a lable column to a df. Custom for plotting really.  
add_type <- function(mean_expression, ids, name){
  mask <- mean_expression$id %in% ids
  mean_expression[mask,]$type <- name
  return(mean_expression)
}


# Custom: Barplot of the number of exons for novel lncRNAs 
barplot_exon_count <- function(gr, type, col, border=NA, size= 0){
  if(is.na(border)){
    border <- col
  }
  ## extract the number of exons
  gr <- gr[gr$type =="exon",]
  gr <- get_only_max_transcript(gr)
  df_l <- data.frame(get_nr_exons(gr))
  df_l$type <- type
  
  h <- ggplot(df_l, aes(x=max_exon, fill=type, col = type)) + 
    geom_histogram(position="identity", binwidth =1)+xlim(1,15)
  
  h_plotdata <- ggplot_build(h)$data[[1]]
  h_plotdata$group <- as.factor(h_plotdata$group)
  levels(h_plotdata$group) <- c(type)
  sizes <- c(rep(1.5,2),rep(0.5,4))
  ## plot with geom_bar
  p1 <-ggplot(h_plotdata, aes(x=x, y=y, fill = group, col = group, size = group)) +
    geom_bar(stat = "identity", width = 0.8) +
    theme(legend.title=element_blank())+
    labs(y = "", x = "")+
    theme(legend.title=element_blank())+ theme(legend.position = "none")+
    scale_x_continuous( labels = as.character(h_plotdata$x), breaks = (h_plotdata$x)) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size =12))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=c(col))+scale_y_continuous(expand = c(0,1),breaks = c(round_any(max(h_plotdata$y), 100, f = floor)))+
    scale_color_manual(values=c(border))+scale_size_manual(values = c(size))+
    theme(panel.grid.major = element_blank(), axis.text = element_text(colour = "black", size = 12),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.y = element_line())+theme(axis.ticks.length=unit(.2, "cm"))

  return(p1)
}

# Input: GenomicRanges object
# Outout: Exon length summary (Data frame)
# Obtain length of the exon 
calc_exon_length <- function(gr, type){
  gr <- gr[gr$type =="exon",]
  gr <- get_only_max_transcript(gr)
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "exon_number" = gr$exon_number, "range_width" = width(ranges(gr)))
  collapsed <- df%>% dplyr::group_by(gene_id,transcript_id, exon_number) %>% dplyr::summarize("range" = sum(range_width))
  collapsed["type"] <- type
  return(collapsed)
}

# Input: GenomicRanges object
# Outout: Transcript length summary (Data frame)
# Obtain length of the transcript 
calc_transcript_length <- function(gr, type){
  gr <- gr[gr$type =="exon",]
  gr <- get_only_max_transcript(gr)
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "range_width" = width(ranges(gr)))
  collapsed <- df%>% dplyr::group_by(gene_id,transcript_id) %>% dplyr::summarize("range" = sum(range_width))
  collapsed["type"] <- type
  collapsed$gene_id <- as.character(collapsed$gene_id)
  collapsed$transcript_id <- as.character(collapsed$transcript_id)
  return(collapsed)
}

# Input: GenomicRanges object
# Output: GenomicRanges object
# Def: deplete all monoexonic transcripts 
remove_one_exon <- function(gr){
  df_l <- data.frame(get_nr_exons(gr[gr$type =="exon",]))
  depleted <- gr[gr$gene_id %in% df_l[df_l$max_exon > 1,]$gene_id,]
  return(depleted)
}

# Custom: Boxplot of expression
plot_expression <- function(max_expression, palette_expression = palette_expression,   title = "Maximum Expression", level=  c("Novel lncRNAs","Annotated lncRNAs", "mRNAs"), palette_border = NA, sizes = NA){
  if(is.na(palette_border)){
    palette_border <- rep("black",3)
  }
  if(is.na(sizes)){
    sizes = rep(0.5, length(level))
  }
  print(palette_border)
  plot <- ggplot(max_expression, aes(x=factor(type,level = level) , y=expr, fill = factor(type,level = level),color =factor(type,level = level),  size = factor(type,level = level) )) + 
    geom_boxplot(outlier.shape = NA, alpha = 1)+
    scale_fill_manual(values =palette_expression)+
    scale_color_manual(values = palette_border)+
    labs(y = "median expression log(TPM)", x = "")+
    theme(legend.title=element_blank(), axis.text = element_text(color = "black"))+
    theme(plot.title = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "darkgrey"))+
    theme(panel.background = element_rect(fill = "white"))+theme(legend.position = "none")+scale_size_manual(values = sizes)
  return(plot)
}



# Custom, summary of multiple plots. 
plot_stats_annotation <- function (novel_expressed,lncRNAs_ref,lncRNAs_ref_human, mrna_ref_human, mRNAs_ref,df, palette = c("#E8C2D8", "#D4549C", "#900051", "#8195D7", "navy"), size = 20, maxquantile= 0.9){

  levels <- c('Novel LncRNAs', 'Annotated LncRNAs - Macaque ', 'Annotated LncRNAs - Human',"Annotated mRNAs - Macaque", "Annotated mRNAs - Human" )
  # --------------------
  ## EXON COUNT
  # --------------------
  
  ec1 <- barplot_exon_count(novel_expressed, "Novel", palette[1])+
                    theme(axis.ticks.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.line.x = element_blank(), text  =element_text(size = size))
  ec2 <- barplot_exon_count(lncRNAs_ref, "lncRNAs - Reference Macaque", palette[2])+
                      theme(axis.ticks.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.line.x = element_blank(), text  =element_text(size = size))
  ec3 <- barplot_exon_count(lncRNAs_ref_human, "lncRNAs - Reference Human", palette[3])+
                        theme(axis.ticks.x = element_blank(),
                              axis.text.x = element_blank(),
                              axis.line.x = element_blank(), text  =element_text(size = size))
  ec4 <- barplot_exon_count(mrna_ref_human, "mRNAs - Reference Human", palette[4])+
                  theme(axis.ticks.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.line.x = element_blank(), text  =element_text(size = size))
  ec5 <- barplot_exon_count(mRNAs_ref, "mRNAs - Reference Macaque", palette[5])+theme(text  =element_text(size = size))
                      
  a <- ggarrange( ec1,ec2,ec3,ec4,ec5,  ncol=1, nrow=5, heights = c(1,1,1,1,1.5)) 
  a <- annotate_figure(a, bottom = text_grob("number of exons per gene", size  = size), left = text_grob("number of genes", size = size, rot = 90))

  print("---Number of exons computed -----")
  # --------------------
  ## Transcript lengths - Boxplot
  # --------------------
  df <- data.frame()
  df <- rbind(df,data.frame(calc_transcript_length(novel_expressed, "Novel LncRNAs")))
  df <- rbind(df,data.frame(calc_transcript_length(lncRNAs_ref, "Annotated LncRNAs - Macaque ")))
  df <- rbind(df,data.frame(calc_transcript_length(lncRNAs_ref_human, "Annotated LncRNAs - Human")))
  df <- rbind(df,data.frame(calc_transcript_length(mrna_ref_human, "Annotated mRNAs - Human")))
  df <- rbind(df,data.frame(calc_transcript_length(mRNAs_ref, "Annotated mRNAs - Macaque")))
  
  ylab <- seq(0,10,1)
  p <- ggplot(df, aes(x = factor(type, level = levels),  y = range )) +
    labs( x = "", y = "transcript length (bp)" )+
    geom_boxplot(outlier.shape=NA, fill = alpha(palette, 1),color = "black", na.rm = TRUE) +
    scale_y_continuous(limits = c(0,ceiling(max(quantile(df$range, c(0,maxquantile))) / 1000))*1000, labels = paste0(ylab, "000"),breaks = 10^3 * ylab, expand = c(0,0))+
    theme(axis.text.x = element_blank(),axis.text.y = element_text(size = size, color = "black"), axis.title = element_text(size = size))+theme(panel.background = element_rect(fill = "white"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(), 
          axis.line.y = element_line())+theme(axis.ticks.length=unit(.2, "cm"))

  print("---Transcript length computed -----")
  # --------------------
  ## Exon lengths - Boxplot
  # --------------------
  df <- data.frame()
  df <- rbind(df,data.frame(calc_exon_length(novel_expressed, "Novel LncRNAs")))
  df <- rbind(df,data.frame(calc_exon_length(lncRNAs_ref, "Annotated LncRNAs - Macaque ")))
  df <- rbind(df,data.frame(calc_exon_length(lncRNAs_ref_human, "Annotated LncRNAs - Human")))
  df <- rbind(df,data.frame(calc_exon_length(mrna_ref_human, "Annotated mRNAs - Human")))
  df <- rbind(df,data.frame(calc_exon_length(mRNAs_ref, "Annotated mRNAs - Macaque")))
  
  ylab <- seq(0,1001,200)
  p1 <- ggplot(df, aes(x = factor(type, level = levels),  y = range )) +
    labs( x = "", y = "exon length (bp)" )+
    geom_boxplot(outlier.shape=NA, fill = alpha(palette, 1),color = "black", na.rm = TRUE) +
    scale_y_continuous(limits = c(0,ceiling(max(quantile(df$range, c(0,maxquantile))) / 200))*200, labels = ylab,breaks =ylab, expand = c(0,0))+
    theme(axis.text.x = element_blank(),axis.text.y = element_text(size = size, color = "black"), axis.title = element_text(size = size))+theme(panel.background = element_rect(fill = "white"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(), 
          axis.line.y = element_line())+theme(axis.ticks.length=unit(.2, "cm"))
  
  print("---Exon length computed -----")
  
  return(list(a,p,p1))
}

plot_stats_annotation_separated <- function (novel_expressed_poly, novel_expressed_ribo,lncRNAs_ref,lncRNAs_ref_human, mrna_ref_human, mRNAs_ref,palette, palette_border, maxquantile = 0.9, size = 20){
  levels <- c('Intergenic Novel LncRNAs', 'Antisense Novel LncRNAs', 'Annotated LncRNAs - Macaque ', 'Annotated LncRNAs - Human',"Annotated mRNAs - Macaque", "Annotated mRNAs - Human" )
  # --------------------
  ## EXON COUNT
  # --------------------
  
  ec1 <- barplot_exon_count(novel_expressed_poly, "Intergenic Novel LncRNAs", palette[1], palette_border[1], 1)+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank())
  ec1b <- barplot_exon_count(novel_expressed_ribo, "Antisense Novel LncRNAs", palette[2], palette_border[2],1)+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank())
  ec2 <- barplot_exon_count(lncRNAs_ref, "lncRNAs - Reference Macaque", palette[3], palette_border[3])+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank())
  ec3 <- barplot_exon_count(lncRNAs_ref_human, "lncRNAs - Reference Human", palette[4], palette_border[4])+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank())
  ec4 <- barplot_exon_count(mrna_ref_human, "mRNAs - Reference Human", palette[5],palette_border[5])+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank())
  ec5 <- barplot_exon_count(mRNAs_ref, "mRNAs - Reference Macaque", palette[6],palette_border[6])
  
  a <- ggarrange( ec1, ec1b,ec2,ec3,ec4,ec5,  ncol=1, nrow=6, heights = c(1,1,1,1,1, 1.5)) 
  a <- annotate_figure(a, bottom = text_grob("number of exons per gene", size  = size), left = text_grob("number of genes", size = size, rot = 90))
  

  # --------------------
  ## Transcript lengths - Boxplot
  # --------------------
  df <- data.frame()
  df <- rbind(df,data.frame(calc_transcript_length(novel_expressed_poly, "Intergenic Novel LncRNAs")))
  df <- rbind(df,data.frame(calc_transcript_length(novel_expressed_ribo, "Antisense Novel LncRNAs")))
  df <- rbind(df,data.frame(calc_transcript_length(lncRNAs_ref, "Annotated LncRNAs - Macaque ")))
  df <- rbind(df,data.frame(calc_transcript_length(lncRNAs_ref_human, "Annotated LncRNAs - Human")))
  df <- rbind(df,data.frame(calc_transcript_length(mrna_ref_human, "Annotated mRNAs - Human")))
  df <- rbind(df,data.frame(calc_transcript_length(mRNAs_ref, "Annotated mRNAs - Macaque")))
  
  sizes <- c(rep(1.5,2),rep(0.5,4))
  ylab <- seq(0,10,1)
  p <- ggplot(df, aes(x = factor(type, level = levels),  y = range, size = type)) +
    labs( x = "", y = "transcript length (bp)" )+
    geom_boxplot(outlier.shape=NA, fill = alpha(palette, 1),color = palette_border, na.rm = TRUE,size = sizes) +
    scale_y_continuous(limits = c(0,ceiling(max(quantile(df$range, c(0,maxquantile))) / 1000))*1000, labels = paste0(ylab, "000"),breaks = 10^3 * ylab, expand = c(0,0))+
    theme(axis.text.x = element_blank(),axis.text.y = element_text(size = size, color = "black"), axis.title = element_text(size = size))+theme(panel.background = element_rect(fill = "white"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(), 
          axis.line.y = element_line())+theme(axis.ticks.length=unit(.2, "cm"))+scale_size_manual(values = c(sizes))
  
  
  # Just for having the palette 
  pal <- ggplot(df, aes(x = factor(type, level = levels),  y = range, fill =factor(type, level = levels), color = factor(type, level = levels) )) +
  geom_density()+scale_fill_manual(values = palette)+scale_color_manual(values = palette_border)
  
  
  # --------------------
  ## Exon lengths - Boxplot
  # --------------------
  df <- data.frame()
  df <- rbind(df,data.frame(calc_exon_length(novel_expressed_poly, "Intergenic Novel LncRNAs")))
  df <- rbind(df,data.frame(calc_exon_length(novel_expressed_ribo, "Antisense Novel LncRNAs")))
  df <- rbind(df,data.frame(calc_exon_length(lncRNAs_ref, "Annotated LncRNAs - Macaque ")))
  df <- rbind(df,data.frame(calc_exon_length(lncRNAs_ref_human, "Annotated LncRNAs - Human")))
  df <- rbind(df,data.frame(calc_exon_length(mrna_ref_human, "Annotated mRNAs - Human")))
  df <- rbind(df,data.frame(calc_exon_length(mRNAs_ref, "Annotated mRNAs - Macaque")))
  ylab <- seq(0,1001,200)
  
  p1 <- ggplot(df, aes(x = factor(type, level = levels),  y = range, size = type )) +
    labs( x = "", y = "exon length (bp)" )+
    geom_boxplot(outlier.shape=NA, fill = alpha(palette, 1),color = palette_border, na.rm = TRUE,  size = sizes) +
    scale_y_continuous(limits = c(0,ceiling(max(quantile(df$range, c(0,maxquantile))) / 200))*200, labels = ylab,breaks =ylab, expand = c(0,0))+
    theme(axis.text.x = element_blank(),axis.text.y = element_text(size = size, color = "black"), axis.title = element_text(size = size))+theme(panel.background = element_rect(fill = "white"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(), 
          axis.line.y = element_line())+theme(axis.ticks.length=unit(.2, "cm"))+scale_size_manual(values = c(sizes))
  return(list(a,p,p1, pal))
}

# Plot tissue specificity 
barplot_tissues <- function(df, type, col, border=NA, size = 0.5){
  if(is.na(border)){
    border <- col
  }
  ## extract the number of exons
  df_l <- df[df$type == type,]
  h <- ggplot(df_l, aes(x=n_tissues_expresseing_gene, fill=type)) + 
    geom_histogram(position="identity", binwidth =1)+
    xlim(1,13)
  
  h_plotdata <- ggplot_build(h)$data[[1]]
  h_plotdata$group <- as.factor(h_plotdata$group)
  levels(h_plotdata$group) <- c(type)
  
  ## plot with geom_bar
  p1 <-ggplot(h_plotdata, aes(x=x, y=y, fill = group, col = group, palette = col , size = group)) +
    geom_bar(stat = "identity", width = 0.8) +
    theme(legend.title=element_blank())+
    labs(y = "", x = "")+
    theme(legend.title=element_blank())+ theme(legend.position = "none")+
    scale_x_continuous( labels = as.character(h_plotdata$x), breaks = (h_plotdata$x)) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size =15, color = "black"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=c(col))+scale_y_continuous(expand = c(0,1),breaks = c(round_any(max(h_plotdata$y), 100, f = floor)))+
    scale_size_manual(values = c(size))+
    scale_color_manual(values = c(border))
  return(p1)
}




exon_nr_summary <- function(gr){
  get_nr_exons_mod<- function(gr){
    df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id,"exon_number" = as.numeric(gr$exon_number))
    number_exons <- df %>% dplyr::group_by(gene_id,transcript_id) %>%dplyr::summarize(max_exon = max(exon_number))
    return(number_exons)
  }
  gr <- gr[gr$type =="exon",]
  #gr <- get_only_max_transcript(gr)
  df_l <- data.frame(get_nr_exons_mod(gr))
  return(df_l)
}


transcript_length_summary <- function(gr){
  gr <- gr[gr$type =="exon",]
  df <- data.frame("gene_id" = gr$gene_id,"transcript_id" = gr$transcript_id, "range_width" = width(ranges(gr)))
  collapsed <- df%>% dplyr::group_by(gene_id,transcript_id) %>% dplyr::summarize("range" = sum(range_width))
  return(collapsed)
}


# ---- OLD 

bed12_entry <- function(transcript){
  
  entry <- gtf_no_gene[gtf_no_gene$transcript_id == transcript, ]
  exons <- entry[entry$type == "exon",]
  
  # Extract all informations 
  chrom <- as.character(seqnames(entry))[1]
  start <- min(start(ranges(exons)))-1
  stop <- max(end(ranges(exons)))
  name <- exons$transcript_id[1]
  score <- "-"
  strand <- as.character(strand(exons))[1]
  thickstart<- start
  thickend <- stop
  rgb <- 0 
  blockcount <- length(exons) 
  blocksizes <- paste(width(ranges(exons)), collapse = ",")
  blockstarts <- paste(start(ranges(exons))-1-start , collapse= ",")
  line <- data.frame(chrom, start, stop, name, score, strand, thickstart, thickend, rgb, blockcount, blocksizes, blockstarts)
  
  return(line)
}


# ------------------------------------------------
#  Utils for lncRNA annotation stats
# ------------------------------------------------
plot_stats_benchmark <- function(pcvslnc, palette = c("orange","#994C00")){
  data <- data.frame(table(pcvslnc[pcvslnc$type == "transcript",]$class_code))
  colnames(data) <- c("group", "value")
  data <- data[data$group %in% c("u", "i", "x"), ]
  # Basic piechart
  # Compute the position of labels
  data <- data %>% 
    arrange(desc(group)) %>%
    mutate(prop = value / sum(data$value) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  data$label <- ifelse(data$group == "x", "antisense", "intergenic")
  data$label <- paste(data$label, round(data$prop, digits = 2 ) , sep = "\n")
  # Basic piechart
  plotassembly <- ggplot(data, aes(x="", y=prop, fill=group)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="none") +
    geom_text(aes(y = ypos, label = label), color = "white", size=6) +
    scale_fill_manual(values=palette)
  return(plotassembly)
}


expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}
