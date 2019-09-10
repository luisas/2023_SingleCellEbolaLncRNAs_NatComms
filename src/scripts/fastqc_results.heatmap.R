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
inpath="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/02_fastqc/Zyagen_samples/"
sample_info=read.table("/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/00_InformationFiles/sample_info.Zyagen_test.txt")
colnames(sample_info)<-c("fastq_file","path",
                         "sample","tissue",
                         "lane","method",
                         "stranded")

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
for(row in 1:nrow(sample_info)){
  # r1
  # Creates path to fastqc folder
  fastq_summary<-paste0(inpath,sample_info[row,"tissue"],"/",sample_info[row,"sample"],"/",
                        sample_info[row,"fastq_file"],".1_fastqc/summary.txt")
  qc_res<-read.delim(fastq_summary,sep="\t",header = F)
  qc_res$number<-as.numeric(gsub("FAIL",0,gsub("WARN",1,gsub("PASS",2,qc_res$V1))))
  sample_data_r1<-c(sample_info[row,"tissue"],
                    sample_info[row,"sample"],
                    sample_info[row,"lane"],
                    "R1",
                    ifelse(sample_info[row,"tissue"]=="Blood",sample_info[row,"method"],"NA"),
                    qc_res$number)
  sample_data_r2 <- sample_data_r1
  sample_data_r1[4] <- "R2"
  res_fastqc[[length(res_fastqc)+1]]<-sample_data_r1
  res_fastqc[[length(res_fastqc)+1]]<-sample_data_r2
}

res_fastqc_df<-do.call(cbind.data.frame,res_fastqc)
colnames(res_fastqc_df)<-""
rownames(res_fastqc_df)<-c("Tissue","Sample","Lane","ReadPair","Method",hm_colnames)



# ----------------------------------------  Define palette colors for plot -------------------------------------------------

# Legend colors
qc_cols<-c("#008744","#ffa700","#d62d20")
names(qc_cols)<-c(2,1,0)

# Tissue Colors
tissue_cols<-c("#a30041","#468499")
names(tissue_cols)<-c("Adrenal","Brain")
# tissue_cols<-c("#2B061E","#775253", "#4392F1","#7180AC", "#2B4570","#8AEA92", "#F6E27F","#E2C391", "#EB5160", "#7E52A0" )
# names(tissue_cols)<-c("Adrenal","Brain","Kidney","LN", "Liver","Ovary","Skin","SpinalCord","Spleen", "Testis")

# Sample colors
samples<-unique(as.factor(as.character(res_fastqc_df["Sample",])))
if (length(samples)  == 1){ 
  sample_cols <- c("#FFA552")
  names(sample_cols)<-levels(samples)
}else if(length(samples)  == 2) { 
  sample_cols <- c("#FFA552", "381D2A")
  names(sample_cols)<-levels(samples)
}else{
  sample_cols<-brewer.pal(length(samples),"Paired")
  names(sample_cols)<-samples
}

# Read Pair Colors
readpair_cols<-brewer.pal(3,"Set1")[c(2:3)]
names(readpair_cols)<-c("R1","R2")

# Lane Colors
lanes<-unique(as.factor(as.character(res_fastqc_df["Lane",])))
lanes_cols<-c("#bada55","#ff6666")
names(lanes_cols)<-c("1","2")

# Metadata Colors
metadata<-unique(as.factor(as.character(res_fastqc_df["Method",])))
length(metadata)

if (length(metadata)  == 1){ 
  metadata_cols <- c("#EAC435")
  names(metadata_cols)<-levels(metadata)
}else if(length(metadata)  == 2) { 
  metadata_cols <- c("#EAC435", "#E40066")
  names(metadata_cols)<-levels(metadata)
}else{
  metadata_cols<-brewer.pal(length(metadata),"Dark2") # c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC")#
  names(metadata_cols)<-metadata
  }
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
                                                 Metadata = as.factor(as.character(res_fastqc_df["Method",])),
                                                 col = list( Tissue = tissue_cols,
                                                             Sample = sample_cols,
                                                             Lane = lanes_cols,
                                                             ReadPair = readpair_cols,
                                                             Metadata = metadata_cols) ),
              cluster_rows=T,cluster_columns = T,#d696bb
              show_column_names = F)
png(paste0(inpath,"fastqc_resLuisa.png"),width = 900,height = 400)
draw(ht)
dev.off()

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
                                                 Metadata = as.factor(as.character(res_fastqc_df["Method",])),
                                                 col = list( Tissue = tissue_cols,
                                                             Sample = sample_cols,
                                                             Lane = lanes_cols,
                                                             ReadPair = readpair_cols,
                                                             Metadata = metadata_cols),
                                                 show_annotation_name = T),
              cluster_rows=F,cluster_columns = F,#d696bb
              show_column_names = F)
#column_title = "No clustering")
png(paste0(inpath,"fastqc_res.png"),width = 900,height = 400)
draw(ht)
dev.off()


