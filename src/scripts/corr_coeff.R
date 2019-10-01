# libraries 
library(stringr)
library(zeallot)
library(pheatmap)

# Input files
input_dir = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/Ebola_Raquel/hisat2/Zyagen_samples"
htseq_file_pattern = "HTseq.gene_counts.tab"

# Utils
iterate_files <- function(inpath, pattern_string){
  files <- list.files(path=inpath, pattern= pattern_string, full.names=TRUE, recursive=TRUE)
  return(files)
} 


extract_info_sample_from_filename <- function(file){
  file_no_ext <- strsplit(fpkm_file,".", fixed = T)[[1]][1]
  info <- strsplit(tail(strsplit(file_no_ext,"/", fixed = T)[[1]],1),"_", fixed = T)[[1]]
  tissue <- info[2]
  sample <- info[4]
  sample_name=paste(tissue,sample,sep="_")
  return(c(tissue,sample, sample_name))
}

htseq_count_files <-iterate_files(input_dir,htseq_file_pattern) 

# FPKM files
fpkm_files<-list("/Users/luisasantus/Desktop/Zyagen_Adrenal_D000_160421.tab","/Users/luisasantus/Desktop/Zyagen_Adrenal_D000_1604212.tab","/Users/luisasantus/Desktop/Zyagen_Adrenal_D000_1000.tab")
fpkms<-list()
for (fpkm_file in fpkm_files){
   c(tissue, sample,complete_id) %<-% extract_info_sample_from_filename(fpkm_file)
   fpkms[[complete_id]]<- read.table(fpkm_file, sep="\t",header=TRUE,row.names = 1, as.is=TRUE)
}
fpkms_df<-do.call(cbind.data.frame,fpkms)
colnames(fpkms_df)<-names(fpkms)
fpkms_df
# Calc Correlation Coefficient
pheatmap(cor(fpkms_df,method="pearson"))



# MDS plots
# NOT DONE!!

library(edgeR)
d0 <- DGEList(fpkms_df)
d0 <- calcNormFactors(d0)
# Filter low expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
ebov_genes_index<-c(6410:6419)
ebov_genes_index %in% drop # ebov genes are expressed aboved threshold so they are the last 10 genes from d
#dim(d0)
#dim(d)
png(paste0("/Users/luisasantus/Desktop/","MDS.top_500_genes.png"),width = 600,height = 500)
mds_xy<-plotMDS(fpkms_df)
text(1.7,-1.5,"D7-NEC",cex=0.75)     
text(1.9,-1.2,"D6-NEC",cex=0.75)     
text(2.75,0.1,"BL_xGen-012",cex=0.75)
text(2.7,0.8,"BL_xGen-013",cex=0.75)
text(2.75,0.45,"BL_xGen-09",cex=0.75)
text(2.8,0.62,"BL_xGen-010",cex=0.75)
text(2.7,0.55,"D0",cex=0.75)

# text(mds_xy$cmdscale.out[,1][seq(1,28,by=2)],
#      mds_xy$cmdscale.out[,2][seq(1,28,by=2)],
#      labels=sapply( colnames(counts_df),
#                 function(sample) paste(unlist(strsplit(sample,split="\\."))[3:length( unlist(strsplit(sample,split="\\.")))], collapse="."  ))[seq(1,28,by=2)],
#      cex=0.75,
#      col="black")

legend("top",c(names(sample_cols),"Liver","Blood"),col=c(sample_cols,"grey","grey"),bty='n',pch=c(rep(16,9),15,16),ncol=4)
dev.off()

par(mfrow=c(1,2))
mds_xy<-plotMDS(d,
                col=sapply(colnames(counts_df), function(sample) sample_cols[ unlist(strsplit(sample,split = "\\."))[[1]] ] ),
                pch=sapply(colnames(counts_df), function(sample) ifelse( unlist(strsplit(sample,split = "\\."))[[2]]=="Liver",15,16 )),
                cex=1,
                top=500,main="Including EBOV genes")
mds_xy<-plotMDS(d[c(1:5491),],
                col=sapply(colnames(counts_df), function(sample) sample_cols[ unlist(strsplit(sample,split = "\\."))[[1]] ] ),
                pch=sapply(colnames(counts_df), function(sample) ifelse( unlist(strsplit(sample,split = "\\."))[[2]]=="Liver",15,16 )),
                cex=1,
                top=500,main="Excluding EBOV genes")



