# ---- Testing
library(DESeq2)
library(rtracklayer)
library(edgeR)
library(zeallot)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(data.table)
library(gridExtra)
# -------------

# Load Utils
scripts_dir = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/code/ebola/src/scripts/"
source(paste0(scripts_dir,"de_utils.R"))
source(paste0(scripts_dir,"mapping_qa.R"))
# Temp 
directory <- "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/03_hisat"
directory_zyagen <- "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/03_hisat/Zyagen"
outpath <- "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/plots"

# ------------------------------------------------
#       BARPLOT UMIS AND FILTER READ COUNTS 
# ------------------------------------------------
# Get Reads raw Counts after different filtering procedures
# Y axis corresponds to the average number of reads of the samples in one tissue.
# X axis to the tissue
# Color to the dataset 
filter_count_files <- iterate_files(directory_zyagen, "f3.q60.n_reads.txt")
umi_count_files <- iterate_files(directory_zyagen, "f3.q60.umi_dedup.n_reads.txt")
md_count_files <- iterate_files(directory, "f3.q60.md.n_reads.txt")

#png(paste0(outpath,"/umis_barplot.png"),width = 1400,height = 500)
#options(scipen=999)
#par(mfrow=c(1,2))


# Barplot of the total count 
data <-rbind(extract_total(umi_count_files,"umi"),extract_total(filter_count_files,"filter"), extract_total(md_count_files,"md"))
p1 <- ggplot(data, aes(fill=method, y=total, x=id)) + 
  geom_bar( position="dodge", stat = "identity")+
  xlab("Sample ID")+
  ylab("Total PE reads left")+
  ggtitle("Number of Reads left after removing Duplicates")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Barplot the difference  
umis <- extract_total(umi_count_files,"umi")
filter <- extract_total(filter_count_files,"filter")
mds <- extract_total(md_count_files,"md")
a <- cbind(umis$id,as.numeric(filter$total - umis$total), "Filter-Umi")
b <- cbind(umis$id,as.numeric(filter$total - mds$total), "Filter-Md")
data_diff<-data.frame(rbind(a,b))
names(data_diff) <- c("id","difference", "comparison")
p2 <- ggplot(data_diff, aes(fill=comparison, y=as.numeric(difference), x=id)) + 
  geom_bar( position="dodge", stat = "identity")+
  xlab("Sample ID")+
  ylab("Difference of PE reads left")+
  ggtitle("Difference of Number of Reads left after removing Duplicates")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
png(paste0(outpath,"/umis_barplot.png"),width = 1400,height = 500)
options(scipen=999)
grid.arrange(p1, p2, nrow = 1)
dev.off()

# ------------------------------------------------
#                 BOX PLOTS UMIS
# ------------------------------------------------
# we need the normalized values 

library(Rsubread)

iterate_files <- function(inpath, pattern_string){
  files <- list.files(path=inpath, pattern= pattern_string, full.names=TRUE, recursive=TRUE)
  return(files)
} 

hisat_bam_files <- iterate_files(directory_zyagen, "f3.q60.umi_dedup.bam$")
annotation <- "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/01_PreliminaryFiles/gene_annotations/rheMac8_EBOV-Kikwit.gtf"
featureCounts(hisat_bam_files[2], annot.ext=annotation, isGTFAnnotationFile=TRUE)


featureCounts(hisat_bam_files, annot.ext=annotation, isGTFAnnotationFile=TRUE)
library(descr)
file.head(hisat_bam_files[2])


#---










# --------------------------------------# --------------------------------------
# --------------------------------------# --------------------------------------
# ------
# UMI_DEDUP 
sampleFiles <- iterate_files(directory, "UMI.f3.q60.umi_dedup.HTseq.gene_counts.tab")
sampleFiles_md <- iterate_files(directory, "UMI.f3.q60.md.HTseq.gene_counts.tab")
sampleFiles_filter <- iterate_files(directory, "UMI.f3.q60.HTseq.gene_counts.tab")

# --------------------------------------
dds <- create_dds(sampleFiles)
dds_md <- create_dds(sampleFiles_md)
dds_filter <- create_dds(sampleFiles_filter)
# Obtain the counts for each 
umi_total_count <- sum(assays(dds)$counts)
md_total_count <- sum(assays(dds_md)$counts)
filter_total_count <- sum(assays(dds_filter)$counts)

counts<-c(filter_total_count,md_total_count,umi_total_count)
datasets <- c("filter","umi","md")
data <- data.frame(x,dataset,counts)
ggplot(data, aes(fill=dataset, y=counts, x="a")) + 
  geom_bar( position="dodge", stat = "identity")
# create a dataset

counts <- c(100,90,91)
x <- "a"
data <- data.frame(x,dataset,counts)
data
# Stacked
p <- ggplot(data, aes(fill=dataset, y=counts, x="a")) + 
  geom_bar( position="identity", stat = "identity")
p + coord_cartesian(ylim=c(85,100))




barplot(table(as.matrix(data.frame(umi_total_count/filter_total_count, md_total_count/filter_total_count))))


barplot(c(filter_total_count, umi_total_count, md_total_count),function(x) x/sum(x),names.arg=c("Filter", "UMI", "DUP"), las = 1, xpd = FALSE,  col=brewer.pal(3,"Paired"), ylim=c(md_total_count - 1000000,filter_total_count+100000) )
barplot(as.matrix(apply(distr_df,2,function(x) as.numeric(x)/sum(as.numeric(x)))),las=2,col=brewer.pal(10,"Paired"))

calc_and_add_summarycounts <- function(summary, dds, method){
  for ( tissue in levels(dds$tissue)){
    dds_tissue <- dds[,dds$tissue == tissue]
    # Check if they have more than one sample
    if (length(dds_tissue$id) > 1){
      dds_tissue_sample <- dds_tissue[,dds_tissue$id == 190507]
      dds_tissue_sample2 <- dds_tissue[,dds_tissue$id == 160421]
      diff <- abs(assays(dds_tissue_sample)$counts - assays(dds_tissue_sample2)$counts)
      tissue_summary <- data.frame("dif" =  as.vector(as.integer(diff)), "tissue" = as.factor(tissue), "method" = as.factor(method))
      summary <- rbind(summary,tissue_summary)
    }
  }
  return(summary)
}
tissue_summary_1 <- list()
calc_and_add_summarycounts_removezero <- function(summary, dds_1, dds_2){
  for ( tissue in levels(dds$tissue)){
    dds_tissue_1 <- dds_1[,dds_1$tissue == tissue]
    dds_tissue_2 <- dds_2[,dds_2$tissue == tissue]
    # Check if they have more than one sample
    if (length(dds_tissue_1$id) > 1){
      dds_tissue_1 <- dds_tissue_1[rowSums(counts(dds_tissue_1)) > 0]
      dds_tissue_2 <- dds_tissue_2[rowSums(counts(dds_tissue_2)) > 0]
      dds_tissue_sample_1_1 <- dds_tissue_1[,dds_tissue_1$id == 190507]
      dds_tissue_sample_1_2 <- dds_tissue_1[,dds_tissue_1$id == 160421]
      dds_tissue_sample_2_1 <- dds_tissue_2[,dds_tissue_2$id == 190507]
      dds_tissue_sample_2_2 <- dds_tissue_2[,dds_tissue_2$id == 160421]
      
    
      diff_1 <- abs(assays(dds_tissue_sample_1_1)$counts - assays(dds_tissue_sample_1_2)$counts)
      diff_2 <- abs(assays(dds_tissue_sample_2_1)$counts - assays(dds_tissue_sample_2_2)$counts)
      if(as.numeric(diff_1) > 0 & as.numeric(diff_2 > 0)){
        tissue_summary_1 <- data.frame("dif" =  as.vector(as.integer(diff_1)), "tissue" = as.factor(tissue), "method" = "UMI")
        tissue_summary_2 <- data.frame("dif" =  as.vector(as.integer(diff_2)), "tissue" = as.factor(tissue), "method" = "MD")
        tissue_summary_total <- data.frame("dif_umi" =  as.vector(as.integer(diff_1)),
                                           "dif_md" =  as.vector(as.integer(diff_2)),
                                           "tissue" = as.factor(tissue)
                                           )
        summary <- rbind(summary,tissue_summary_total)
        }
      }
  }
  return(summary)
}

dds_1 <- dds
dds_2 <- dds_md
for ( tissue in levels(dds$tissue)){
  dds_tissue_1 <- dds_1[,dds_1$tissue == tissue]
  dds_tissue_2 <- dds_2[,dds_2$tissue == tissue]
  # Check if they have more than one sample
  if (length(dds_tissue_1$id) > 1){
    dds_tissue_1 <- dds_tissue_1[rowSums(counts(dds_tissue_1)) > 0]
    dds_tissue_2 <- dds_tissue_2[rowSums(counts(dds_tissue_2)) > 0]
    dds_tissue_sample_1_1 <- dds_tissue_1[,dds_tissue_1$id == 190507]
    dds_tissue_sample_1_2 <- dds_tissue_1[,dds_tissue_1$id == 160421]
    dds_tissue_sample_2_1 <- dds_tissue_2[,dds_tissue_2$id == 190507]
    dds_tissue_sample_2_2 <- dds_tissue_2[,dds_tissue_2$id == 160421]
    rownames(diff_1)
    diff_1 <- abs(assays(dds_tissue_sample_1_1)$counts - assays(dds_tissue_sample_1_2)$counts)
    diff_2 <- abs(assays(dds_tissue_sample_2_1)$counts - assays(dds_tissue_sample_2_2)$counts)
    diff_1 <- as.data.frame(diff_1)
    diff_2 <- as.data.frame(diff_2)
    a <- merge(diff_1,diff_2,by="row.names")
    colnames(a)
    if(as.numeric(diff_1) > 0 & as.numeric(diff_2 > 0)){
      tissue_summary_1 <- data.frame("dif" =  as.vector(as.integer(diff_1)), "tissue" = as.factor(tissue), "method" = "UMI")
      tissue_summary_2 <- data.frame("dif" =  as.vector(as.integer(diff_2)), "tissue" = as.factor(tissue), "method" = "MD")
      tissue_summary_total <- data.frame("dif_umi" =  as.vector(as.integer(diff_1)),
                                         "dif_md" =  as.vector(as.integer(diff_2)),
                                         "tissue" = as.factor(tissue)
      )
      summary <- rbind(summary,tissue_summary_total)
    }
  }
}
summary  <- DataFrame()
summary <- calc_and_add_summarycounts_removezero(summary,dds,dds_md)

summary <- summary[summary$tissue == "Spleen",]
summary <- summary[summary$dif_umi !=0 & summary$dif_md != 0, ]
boxplot(list(summary$dif_umi, summary$dif_md ), outline=F)

identical(summary$dif_umi, summary$dif_md)
wilcox.test(summary$dif_umi, summary$dif_md )


summary_1  <- DataFrame()
summary_2  <- DataFrame()
summary_1 <- calc_and_add_summarycounts(summary,dds,"UMI")
summary_2 <- calc_and_add_summarycounts(summary,dds_md,"MD")

cbind(summary_1,summary_2)



as.data.frame(summary)
levels(summary$method)
ggplot(as.data.frame(summary), aes(x=as.factor(tissue), y=dif, color= method)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("Tissue") + 
  xlab("Difference") + 
  ggtitle("UMI vs MD")+ ylim(0,5)
as.vector
umis= calc_and_add_summarycounts(summary,dds,"UMI")

boxplot(list(summary[summary$method=="UMI" & summary$tissue =="Adrenal", "dif"],
            summary[summary$method=="MD" & summary$tissue =="Adrenal", "dif"]),outline=F)
# -------------------------
#     Plotting - old
# -------------------------
out_dir = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/plots/03_hisat/Zyagen/"
plot_counts_heatmap(dds, paste(out_dir,"heatmap_umi_dedup_pearson.png"), "pearson")
plot_counts_heatmap(dds, paste(out_dir,"heatmap_umi_dedup_spearman.png"), "spearman")
plot_counts_heatmap(dds_md, paste(out_dir,"heatmap_md_pearson.png"), "pearson")
plot_counts_heatmap(dds_md, paste(out_dir,"heatmap_md_spearman.png"), "spearman")
plot_counts_heatmap(dds_filter, paste(out_dir,"heatmap_filter_pearson.png"), "pearson")
plot_counts_heatmap(dds_filter, paste(out_dir,"heatmap_filter_spearman.png"), "spearman")

# -- old 


# GTF read
# gtf = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/01_PreliminaryFiles/gene_annotations/rheMac8.ensembl_release97.gtf"
# annotation <- readGFF(gtf)
# annotation$size <- annotation$end -annotation$start
#Extract only genes from the annotation - this is as well what htseq does 
# annotation_genes <- annotation[annotation$type == "gene",]
# gene_lengths<- data.frame(annotation_genes$gene_id, annotation_genes$size)
# colnames(gene_lengths) <- c('gene_id','basepairs')