# ------------------------------------------------------------------------------------------------------------------------------------------------------
# Libraries loading
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ComplexHeatmap)
library(dendextend)
library(RColorBrewer)
# Options
options(stringsAsFactors = F)
# ------------------------------------------------------------------------------------------------------------------------------------------------------

# Variables
#inpath="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/02_fastqc/Zyagen/"
#outpath="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data/02_RNA-Seq/plots/02_fastqc/Zyagen"

plot_heatmap_fast_zyagen <- function(inpath,outpath){
  # FastqQC results column names
  hm_colnames<-c("Basic Statistics",
                 "Per base sequence quality",
                 "Per tile sequence quality",   
                 "Per sequence quality scores",
                 "Per base sequence content",
                 "Per sequence GC content",     
                 "Per base N content",
                 "Sequence Length Distribution",
                 "Sequence Duplication Levels", 
                 "Overrepresented sequences",
                 "Adapter Content","Kmer Content")
  
  res_fastqc<-list()
  files <- list.files(path=inpath, pattern="summary.txt", full.names=TRUE, recursive=TRUE)
  
  for(file in unique(gsub("*.\\d_fastqc/summary.txt", "", files))){
    info <- strsplit(tail(strsplit(file,"/", fixed = T)[[1]],1),"_", fixed = T)[[1]]
    fastq_summary_1<-paste0(file,".1_fastqc/summary.txt")
    qc_res_1<-read.delim(fastq_summary_1,sep="\t",header = F)
    qc_res_1$number<-as.numeric(gsub("FAIL",0,gsub("WARN",1,gsub("PASS",2,qc_res_1$V1))))
    sample_data_r1<-c(info[2],info[4],info[5],"R1",qc_res_1$number)
    
    fastq_summary_2<-paste0(file,".2_fastqc/summary.txt")
    qc_res_2<-read.delim(fastq_summary_2,sep="\t",header = F)
    qc_res_2$number<-as.numeric(gsub("FAIL",0,gsub("WARN",1,gsub("PASS",2,qc_res_2$V1))))
    sample_data_r2<-c(info[2],info[4],info[5],"R2",qc_res_2$number)
    
    res_fastqc[[length(res_fastqc)+1]]<-sample_data_r1
    res_fastqc[[length(res_fastqc)+1]]<-sample_data_r2
  }
  
  
  res_fastqc
  res_fastqc_df<-do.call(cbind.data.frame,res_fastqc)
  colnames(res_fastqc_df)<-""
  rownames(res_fastqc_df)<-c("Tissue","Sample","Lane","ReadPair",hm_colnames)
  
  
  
  # ----------------------------------------  Define palette colors for plot -------------------------------------------------
  
  # Legend colors
  qc_cols<-c("#008744","#ffa700","#d62d20")
  names(qc_cols)<-c(2,1,0)
  
  # Tissue Colors
  tissue_cols<-c("#2B061E","#2B061E","#775253", "#4392F1","#7180AC", "#2B4570","#8AEA92", "#F6E27F","#E2C391", "#EB5160", "#7E52A0" )
  names(tissue_cols)<-c("Blood","Adrenal","Brain","Kidney","LN", "Liver","Ovary","Skin","SpinalCord","Spleen", "Testis")
  
  # Sample colors
  samples<-unique(as.factor(as.character(res_fastqc_df["Sample",])))
  if (length(samples)  == 1){ 
    sample_cols <- c("#FFA552")
    names(sample_cols)<-levels(samples)
  }else if(length(samples)  == 2) { 
    sample_cols <- c("#0000A0", "#ADD8E6")
    names(sample_cols)<-levels(samples)
  }else{
    sample_cols<-brewer.pal(length(samples),"Paired")
    names(sample_cols)<-samples
  }
  
  # Read Pair Colors
  readpair_cols<-brewer.pal(3,"Set1")[c(2:3)]
  names(readpair_cols)<-c("R1","R2")
  
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------
  # HEATMAP FASTQC QA PLOT SUMMARY
  # ------------------------------------------------------------------------------------------------------------------------------------------------------
  ht <- Heatmap(as.matrix(res_fastqc_df[hm_colnames[-1],]),
                col=qc_cols,
                name="",
                heatmap_legend_param = list(
                  at = c(2,1,0),
                  labels = c("PASS", "WARN", "FAIL")),
                top_annotation = HeatmapAnnotation(Tissue = as.factor(as.character(res_fastqc_df["Tissue",])),
                                                   Sample = as.factor(as.character(res_fastqc_df["Sample",])),
                                                   Lane = as.factor(as.character(res_fastqc_df["Lane",])),
                                                   ReadPair = as.factor(as.character(res_fastqc_df["ReadPair",])),
                                                   col = list(Tissue = tissue_cols,
                                                              Sample = sample_cols,
                                                              ReadPair = readpair_cols),
                                                   show_annotation_name = T),
                cluster_rows=T,cluster_columns = T,#d696bb
                show_column_names = F)
  png(paste0(outpath,"/fastqc_res.png"),width = 900,height = 400)
  draw(ht)
  dev.off()
}



