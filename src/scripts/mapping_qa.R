# ------------------------------------------------------------------------------------------------------------------------------------------------------
# Libraries loading
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ComplexHeatmap)
library(dendextend)
library(RColorBrewer)
library(stringr)
# Options
options(stringsAsFactors = F)
# ------------------------------------------------------------------------------------------------------------------------------------------------------
# Accessory functions
parse_hisat2_summaryFile = function(filepath) {
  con = file(filepath, "r")
  summary <- list()
  PE_string <- "were paired"
  PE_conc_1 <- "aligned concordantly exactly 1 time"
  PE_conc_gr_1 <- "aligned concordantly >1 times"
  PE_disc_1 <- "aligned discordantly 1 time"
  mate_1 <- "aligned exactly 1 time"
  mate_gr_1 <- "aligned >1 times"
  name_correspondence <- list()
  name_correspondence[PE_string] <- "PE_total"
  name_correspondence[PE_conc_1] <- "PE aligned concordantly one time"
  name_correspondence[PE_conc_gr_1] <- "PE aligned concordantly more than 1 time"
  name_correspondence[PE_disc_1] <- "PE aligned disconcordantly 1 time"
  name_correspondence[mate_1] <- "mate aligned 1 time"
  name_correspondence[mate_gr_1] <- "mate aligned more than 1 time"
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    for( string_to_be_matched in names(name_correspondence) ){
      if (str_detect(line, string_to_be_matched )){
        key <- as.character(name_correspondence[string_to_be_matched])
        summary[key] <- as.numeric(strsplit(trimws(line, which ="left"),"\\s")[[1]][1])
      }
    }
  }
  return(summary)
  close(con)
}

# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
# HISAT 2 MAPPING RESULTS
# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------
## ---- BARPLOTS MAPPED READS
# ---------------------------------------------------------------------------------------

inpath="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/Ebola_Raquel/hisat2/Zyagen_samples"
outpath = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/code/ebola/src/scripts/"
mapping_stats<-list()
files <- list.files(path=inpath, pattern="summary.txt", full.names=TRUE, recursive=TRUE)

# Extract information from each summary file
for(file in files){
  # Extract tissue, sample and lane from filename of the hisat2 mapping reports
  info <- strsplit(strsplit(file,"/", fixed = T)[[1]][13],".", fixed = T)[[1]]
  tissue <- info[1]
  sample <- info[2]
  lane <- info[3]
  sample_name=paste(sample,tissue,lane,sep=".")

  res_path <- paste(inpath,tissue,sample,sep="/")
  res_file <- paste0(res_path,"/",paste(tissue,sample,lane,"hisat2_summary.txt",sep="."))
  summary <- parse_hisat2_summaryFile(res_file)
  
  # Error handling 
  if(! file.exists(res_file)){print(paste0("File does not exits ",file))}
  
  # Extract the informations
  mapping_stats[[sample_name]]<- c(as.numeric(summary[names(summary)[c(1:4)]])*2, 
                                   # PE reads -> first 4 number multiply by 2, mapped as PE
                                   # hisat1 reports fort 
                                   # PE reads, counting R1 and R1 as 1 read (count)
                                   as.numeric(summary[names(summary)[c(5:6)]])) 
                                   # only one mapped mate of the PE, here one mate (either R1 or R2) is one count
}

# Prepare data structure for plotting
mapping_stats_df<-do.call(cbind.data.frame,mapping_stats)
rownames(mapping_stats_df)<-  mapping_stats[[sample_name]]<-report <- names(summary)
# compute number of UNMAPPED reads
mapping_stats_df["Unmapped",]<-as.numeric(mapping_stats_df[1,])-apply(mapping_stats_df[c(2:6),],2,sum) 

# -------- plot 
create_barplots_reads_mapping <- function(mapping_stats_df) {
  png(paste0(outpath,"hisat2_mapping_report.png"),width = 1200,height = 500)
  options(scipen=999)
  par(mfrow=c(1,2),oma=c(12,4,0,0),mar=c(2,2,2,2),xpd=NA)
  # First Barplot
  barplot(as.matrix(mapping_stats_df[-1,]),yaxt="n", col=brewer.pal(6,"Paired"),las=2,ylim = c(0,32000000), main="Number PE reads")
  axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","),las=2)
  # Second barplot plotting the proportion
  barplot(as.matrix(apply(mapping_stats_df[-1,],2,function(x) x/sum(x) )),las=2,col=brewer.pal(6,"Paired"),ylim = c(0,1.2),yaxt='n', main="Proportion PE Reads")
  axis(2,at=seq(0,1,0.2))
  # Plot
  ?legend
  legend(-50,-0.6,rownames(mapping_stats_df)[-1],bty="n",pch=15,col=brewer.pal(6,"Paired"),ncol=2,cex=1.1)
  dev.off() 
}
create_barplots_reads_mapping(mapping_stats_df)
# -------- finish plot 
# ---------------------------------------------------------------------------------------















# Hisat2 mapping results. RefSeq gene annotation ####
inpath="~/Desktop/mn_cluster/ebola_link/hisat2/Batch_01.RefSeq.KU182905.1/"
mapping_stats.refSeq<-list()
for(row in 1:nrow(sample_info)){
  # hisat2 mappig report
  sample <- sample_info[row,"sample"]
  tissue <- sample_info[row,"tissue"]
  lane <- sample_info[row, "lane"]
  if(tissue == "Blood"){ method <- sample_info[row,"method"]}
  if(tissue == "Blood"){
    sample_name=paste(sample,tissue,method,lane,sep=".")
  }else{
    sample_name=paste(sample,tissue,lane,sep=".")
  }
  res_path <- paste(inpath,tissue,sample,sep="/")
  if(tissue=="Blood"){
    res_file <- paste0(res_path,"/",paste(sample,tissue,method,lane,"hisat2_summary.parsed.txt",sep="."))
  }else{
    res_file <- paste0(res_path,"/",paste(sample,tissue,lane,"hisat2_summary.parsed.txt",sep="."))
  }
  if(! file.exists(res_file)){print(paste0("File does not exits ",row))}
  mapping_stats.refSeq[[sample_name]]<- c(read.csv(res_file,sep="\t",header = F)[,1][c(1:4)]*2, # PE reads -> first 4 number multiply by 2, mapped as PE
                                          # hisat1 reports fort # PE reads, counting R1 and R1 as 1 read (count)
                                          read.csv(res_file,sep="\t",header = F)[,1][c(5:6)]) # only one mapped mate of the PE, here one mate (either R1 or R2) is one count
}
mapping_stats.refSeq_df<-do.call(cbind.data.frame,mapping_stats.refSeq)
rownames(mapping_stats.refSeq_df)<-  mapping_stats.refSeq[[sample_name]]<-report <- read.csv(res_file,sep="\t",header = F)[,2]
mapping_stats.refSeq_df["Unmapped",]<-as.numeric(mapping_stats.refSeq_df[1,])-apply(mapping_stats.refSeq_df[c(2:6),],2,sum) # computer number of unmmapped reads
# doubled check with samtools view -f4 
#identical(as.numeric(apply(mapping_stats.refSeq_df,2,function(x) sum(x[-1]))),as.numeric(mapping_stats.refSeq_df[1,]))
samples_ordered<-unique(sapply(1:ncol(res_fastqc_df), function(col) 
  ifelse(res_fastqc_df[1,col]=="Blood",paste(res_fastqc_df["Sample",col],res_fastqc_df["Tissue",col],res_fastqc_df["Method",col],
                                             res_fastqc_df["Lane",col],sep="."),
         paste(res_fastqc_df["Sample",col],res_fastqc_df["Tissue",col],
               res_fastqc_df["Lane",col],sep="."))))

mapping_stats.refSeq_df<-mapping_stats.refSeq_df[,samples_ordered]


png(paste0(outpath,"hisat2_mapping_report.png"),width = 1200,height = 500)
options(scipen=999)
#layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
par(mfrow=c(1,2),oma=c(12,4,0,0),mar=c(2,2,2,2),xpd=NA)
barplot(as.matrix(mapping_stats.refSeq_df[-1,]),yaxt="n", col=brewer.pal(6,"Paired"),las=2,ylim = c(0,32000000))
axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","),las=2)
#legend("top",rownames(mapping_stats.refSeq_df)[-1],bty="n",pch=15,col=brewer.pal(6,"Paired"),ncol=3,cex=0.7 )
#dev.off()

#png(paste0(inpath,"hisat2_mapping_report.proportion.png"),width = 900,height = 400)
#par(oma=c(12,2,0,0))
barplot(as.matrix(apply(mapping_stats.refSeq_df[-1,],2,function(x) x/sum(x) )),las=2,col=brewer.pal(6,"Paired"),ylim = c(0,1.2),yaxt='n')
axis(2,at=seq(0,1,0.2))

#par(mai=c(0,0,0,0))
#plot.new()
legend(-18,1.3,rownames(mapping_stats.refSeq_df)[-1],bty="n",pch=15,col=brewer.pal(6,"Paired"),ncol=2,cex=0.7)

dev.off()  

# Comparison between Ensembl and RefSeq annotation ####
png(paste0("~/Desktop/mn_cluster/ebola_link/hisat2/","no_pmPE.ensembl_refSeq_comparison.png"),
    width = 600,height = 500)
par(mfrow=c(1,1),oma=c(10,4,0,0))
barplot(as.numeric(mapping_stats_df[2,])-as.numeric(mapping_stats.refSeq_df[2,]),yaxt='n',
        ylab="No. of properly mapped PE reads ",
        main="More complete gene annotation lead to an increased number of mapped reads")
axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","))
axis(1,at=seq(0.75,33.25,length.out = ncol(mapping_stats_df)),labels = colnames(mapping_stats.refSeq_df),las=2)
dev.off()

# MDS plot using edgeR and cpm  ####
# RefSeq genes-> 6409 macaque genes; 10 EBOV genes
counts <- list()
for(row in 1:nrow(sample_info)){
  sample <- sample_info[row,"sample"]
  tissue <- sample_info[row,"tissue"]
  lane <- sample_info[row, "lane"]
  if(tissue == "Blood"){ method <- sample_info[row,"method"]}
  if(tissue == "Blood"){
    sample_name=paste(sample,tissue,method,lane,sep=".")
  }else{
    sample_name=paste(sample,tissue,lane,sep=".")
  }
  counts_path <- paste(inpath,tissue,sample,'counts/',sep="/")
  if(tissue=="Blood"){
    counts_file <- paste0(counts_path,"/",paste(sample,tissue,method,lane,"FPKM.xls",sep="."))
  }else{
    counts_file <- paste0(counts_path,"/",paste(sample,tissue,lane,"FPKM.xls",sep="."))
  }
  counts[[sample_name]] <- read.csv(counts_file,sep="\t",header =T) 
}
sample_names<-vector()
for(row in 1:nrow(sample_info)){
  sample <- sample_info[row,"sample"]
  tissue <- sample_info[row,"tissue"]
  lane <- sample_info[row, "lane"]
  if(tissue == "Blood"){ method <- sample_info[row,"method"]}
  if(tissue == "Blood"){
    sample_name=paste(sample,tissue,method,lane,sep=".")
  }else{
    sample_name=paste(sample,tissue,lane,sep=".")
  }
  sample_names<-c(sample_names,sample_name)
}
counts_df<-do.call(cbind.data.frame,
                   lapply(1:length(counts), function(i) counts[[i]][,"Frag_count"]))
colnames(counts_df)<-sample_names
# rownames(counts_df)<-counts[[1]][,"accession"] # 40 duplicated gene IDs
# RefSeq genes -> only 6409 annoated genes, wo info about gene_type

library(edgeR)
d0 <- DGEList(counts_df)
d0 <- calcNormFactors(d0)
# Filter low expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
ebov_genes_index<-c(6410:6419)
ebov_genes_index %in% drop # ebov genes are expressed aboved threshold so they are the last 10 genes from d
#dim(d0)
#dim(d)
png(paste0(inpath,"MDS.top_500_genes.png"),width = 600,height = 500)
mds_xy<-plotMDS(d,
                col=sapply(colnames(counts_df), function(sample) sample_cols[ unlist(strsplit(sample,split = "\\."))[[1]] ] ),
                pch=sapply(colnames(counts_df), function(sample) ifelse( unlist(strsplit(sample,split = "\\."))[[2]]=="Liver",15,16 )),
                cex=1,
                top=500)
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

# Expression level in FPKM of EBOV genes ####
FPKM <- list()
for(row in 1:nrow(sample_info)){
  sample <- sample_info[row,"sample"]
  tissue <- sample_info[row,"tissue"]
  lane <- sample_info[row, "lane"]
  if(tissue == "Blood"){ method <- sample_info[row,"method"]}
  if(tissue == "Blood"){
    sample_name=paste(sample,tissue,method,lane,sep=".")
  }else{
    sample_name=paste(sample,tissue,lane,sep=".")
  }
  FPKM_path <- paste(inpath,tissue,sample,'counts/',sep="/")
  if(tissue=="Blood"){
    FPKM_file <- paste0(FPKM_path,"/",paste(sample,tissue,method,lane,"FPKM.xls",sep="."))
  }else{
    FPKM_file <- paste0(FPKM_path,"/",paste(sample,tissue,lane,"FPKM.xls",sep="."))
  }
  FPKM[[sample_name]] <- read.csv(FPKM_file,sep="\t",header =T) 
}
fpkm_df<-do.call(cbind.data.frame,
                 lapply(1:length(FPKM), function(i) FPKM[[i]][,"FPKM"]))
colnames(fpkm_df)<-sample_names
ebov_fpkm<-fpkm_df[c(6410:6419),]
rownames(ebov_fpkm)<-FPKM[[1]][c(6410:6419),"accession"]


ebov_fpkm<-ebov_fpkm[,samples_ordered]

png(paste0(inpath,"FPKM.ebov_genes.png"),width = 800,height = 600)
par(mfrow=c(3,4),oma=c(4,0,0,0),mar=c(4,4,4,2),mgp=c(3,1,0))
for(i in 1:10){
  plot(1:28,ebov_fpkm[i,],main=rownames(ebov_fpkm)[i],
       ylab="FPKM",xlab="",yaxt='n',xaxt="n",
       col=sapply(colnames(ebov_fpkm), function(sample) sample_cols[ unlist(strsplit(sample,split = "\\."))[[1]] ] ),
       pch=sapply(colnames(ebov_fpkm), function(sample) ifelse( unlist(strsplit(sample,split = "\\."))[[2]]=="Liver",15,16 ))
  )
  axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","))
  axis(1,at=seq(1.5,14.5,length.out = 7),labels=unique(sapply( colnames(ebov_fpkm),
                                                               function(sample) paste(unlist(strsplit(sample,split="\\."))[3], collapse="."  )))[1:7],
       las=2,cex=0.5,tck=0.015)
}
plot.new()
legend("top",c(names(sample_cols),"Liver","Blood"),col=c(sample_cols,"grey","grey"),bty='n',pch=c(rep(16,9),15,16),ncol=2,cex=1.5)
dev.off()



# Filtered reads####
counts <- list()
for(row in seq(1,nrow(sample_info),by=2)){
  sample <- sample_info[row,"sample"]
  tissue <- sample_info[row,"tissue"]
  lane <- sample_info[row, "lane"]
  if(tissue == "Blood"){ method <- sample_info[row,"method"]}
  if(tissue == "Blood"){
    sample_name=paste(sample,tissue,method,lane,sep=".")
  }else{
    sample_name=paste(sample,tissue,lane,sep=".")
  }
  counts_path <- paste(inpath,tissue,sample,sep="/")
  if(tissue=="Blood"){
    counts_file <- paste0(counts_path,"/",paste(sample,tissue,method,"md.f3.q60.n_reads.txt",sep="."))
  }else{
    counts_file <- paste0(counts_path,"/",paste(sample,tissue,"md.f3.q60.n_reads.txt",sep="."))
  }
  counts[[sample_name]] <- read.csv(counts_file,sep="\t",header =F) 
}
counts_df<-do.call(cbind.data.frame,
                   lapply(1:length(counts), function(i) counts[[i]][,2]))

sample_names<-vector()
for(row in seq(1,nrow(sample_info),by=2)){
  sample <- sample_info[row,"sample"]
  tissue <- sample_info[row,"tissue"]
  #lane <- sample_info[row, "lane"]
  if(tissue == "Blood"){ method <- sample_info[row,"method"]}
  if(tissue == "Blood"){
    sample_name=paste(sample,tissue,method,sep=".")
  }else{
    sample_name=paste(sample,tissue,sep=".")
  }
  sample_names<-c(sample_names,sample_name)
}
colnames(counts_df)<-sample_names
rownames(counts_df)<-counts[[1]][,1]
counts_df<-counts_df[,unique(sapply(samples_ordered,
                                    function(i) paste(unlist(strsplit(i,split="\\."))[1:(length(unlist(strsplit(i,split="\\.")))-1)],collapse = ".")))]

par(mfrow=c(1,2))
par(oma=c(8,2,0,0),xpd=NA)
barplot(as.matrix(counts_df[-c(1,6),]),las=2,col=brewer.pal(5,"Dark2"),yaxt='n',ylim=c(0,38000000))
axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","),las=2)
#legend("top",rownames(counts_df)[-c(1,6)],col=brewer.pal(5,"Dark2"),pch=15,bty="n",ncol=5,cex=0.6)

barplot(as.matrix(apply(counts_df[-c(1,6),],2,function(x) x/sum(x))),las=2,col=brewer.pal(5,"Dark2"))
axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","),las=2)
legend(-17,1.2,rownames(counts_df)[-c(1,6)],col=brewer.pal(5,"Dark2"),pch=15,bty="n",ncol=5)

par(oma=c(8,2,0,0))
barplot(as.matrix(counts_df[6,]),las=2,yaxt='n',main="No. of filtered PE reads mapped to EBOV")
axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","),las=2)


# Read distribution ###
# Read distribution plots ####
inpath="~/Desktop/mn_cluster/ebola_link/hisat2/Batch_01/"
read_dist_stats<-list()
for(row in 1:nrow(sample_info)){
  # hisat2 mappig report
  sample <- sample_info[row,"sample"]
  tissue <- sample_info[row,"tissue"]
  #lane <- sample_info[row, "lane"]
  if(tissue == "Blood"){ method <- sample_info[row,"method"]}
  if(tissue == "Blood"){
    sample_name=paste(sample,tissue,method,sep=".")
  }else{
    sample_name=paste(sample,tissue,sep=".")
  }
  res_path <- paste(inpath,tissue,sample,sep="/")
  if(tissue=="Blood"){
    res_file <- paste0(res_path,"/",paste("/read_distribution/",sample_name,".md.f3.q60.read_distribution.parsed.txt",sep=""))
  }else{
    res_file <- paste0(res_path,"/",paste("/read_distribution/",sample_name,".md.f3.q60.read_distribution.parsed.txt",sep=""))
  }
  if(! file.exists(res_file)){print(paste0("File does not exits ",row))}
  read_dist_stats[[sample_name]]<- read.csv(res_file,sep="\t",header = F)
}
read_dist_stats_df<-cbind(read_dist_stats[[1]][,1],
                          do.call(cbind.data.frame,lapply(1:length(read_dist_stats), function(i) read_dist_stats[[i]][,2] )))
rownames(read_dist_stats_df) <- as.character(read_dist_stats_df[,1])
read_dist_stats_df <- read_dist_stats_df[,-1]
colnames(read_dist_stats_df) <- names(read_dist_stats)
ordered_samples <- unique(sapply(1:ncol(res_fastqc_df), 
                                 function(col) 
                                   ifelse(res_fastqc_df[1,col]=="Liver", paste0(as.character(res_fastqc_df[2,col]),".",as.character(res_fastqc_df[1,col])),
                                          paste0(res_fastqc_df[2,col],".",res_fastqc_df[1,col],".",res_fastqc_df[5,col]))) )
ordered_samples %in% colnames(read_dist_stats_df)
read_dist_stats_df <- read_dist_stats_df[,ordered_samples]
#barplot(as.matrix(read_dist_stats_df[c(4:nrow(read_dist_stats_df)),]))

library(RSkittleBrewer)
getPalette =colorRampPalette(RSkittleBrewer('original'))
#getPalette(10)
#tropical_colors = RSkittleBrewer('tropical')
png(paste0(inpath,"f3.q60.read_distribution.png"),width = 1200,height = 600)
options(scipen=999)
#layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
par(mfrow=c(1,2),oma=c(12,2,4,0),mar=c(2,6,2,2),mgp=c(5,1,0),xpd=NA)
barplot(as.matrix(read_dist_stats_df[c(4:ncol(read_dist_stats_df)),]),
        yaxt="n", col=getPalette(10),las=2,ylab="No. of assigned tags")
axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","),las=2)
#legend("top",rownames(mapping_stats_df)[-1],bty="n",pch=15,col=brewer.pal(6,"Paired"),ncol=3,cex=0.7 )
#dev.off()

#png(paste0(inpath,"hisat2_mapping_report.proportion.png"),width = 900,height = 400)
#par(oma=c(12,2,0,0))
barplot(as.matrix(apply(read_dist_stats_df[c(4:nrow(read_dist_stats_df)),],2,function(x) x/sum(x) )),
        las=2,col=getPalette(10),ylim = c(0,1.2),yaxt='n',ylab="Proportion of assigned tags")
axis(2,at=seq(0,1,0.2))

#par(mai=c(0,0,0,0))
#plot.new()
legend(-10,1.5,rownames(read_dist_stats_df)[c(4:14)],bty="n",pch=15,col=getPalette(10),ncol=3,cex=1)

dev.off()  

read_dist_stats_df<-read_dist_stats_df[,-1]
colnames(read_dist_stats_df)<-names(read_dist_stats)
mapping_stats[[sample_name]]<-report <- read.csv(res_file,sep="\t",header = F)[,2]
mapping_stats_df["Unmapped",]<-as.numeric(mapping_stats_df[1,])-apply(mapping_stats_df[c(2:6),],2,sum) # computer number of unmmapped reads
# doubled check with samtools view -f4 

png(paste0(inpath,"No_of_filtered_read_wo_duplicates_and _number_of_tags.png"),width = 1200,height = 600)
options(scipen=999)
#layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
par(mfrow=c(1,2),oma=c(12,2,4,0),mar=c(2,6,2,2),mgp=c(5,1,0),xpd=NA)
barplot(as.matrix(read_dist_stats_df[1,]),
        yaxt="n",las=2,ylab="No. of non-duplicated filtered PE reads",cex.lab=0.9)
axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","),las=2,cex.axis=0.75)

barplot(as.matrix(rbind(read_dist_stats_df[2,]-read_dist_stats_df[3,],read_dist_stats_df[3,])),
        las=2,ylab="No. of tags",yaxt='n',col=c("black","light grey"))
axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","),las=2,cex.axis=0.75)
legend("topleft",c('Not assigned',"Assigned"),bty='n',cex=0.75,pch=15,col=c("black","light grey"))

dev.off()  

read_dist_stats_df<-read_dist_stats_df[,-1]
colnames(read_dist_stats_df)<-names(read_dist_stats)
mapping_stats[[sample_name]]<-report <- read.csv(res_file,sep="\t",header = F)[,2]
mapping_stats_df["Unmapped",]<-as.numeric(mapping_stats_df[1,])-apply(mapping_stats_df[c(2:6),],2,sum) # computer number o
# DUplicated reads per sample ####
dup_reads <- rbind.data.frame(sapply(seq(1,28,by=2), function(i) sum(mapping_stats_df[2,i],mapping_stats_df[2,i+1])), read_dist_stats_df[1,])
dup_reads[3,] <- dup_reads[1,]-dup_reads[2,]
colnames(dup_reads) <- colnames(read_dist_stats_df)
rownames(dup_reads) <- c("Total","Non-duplicated","Duplicated")

png(paste0(inpath,"duplicated_reads.png"),width = 1200,height = 600)
options(scipen=999)
#layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
par(mfrow=c(1,2),oma=c(12,2,2,0),mar=c(2,6,2,2),mgp=c(5,1,0),xpd=NA)
barplot(as.matrix(dup_reads[c(2,3),]),
        yaxt="n", las=2,ylab="No. of PE reads",col=c("black","grey"))
axis(2,at=axTicks(2),labels = prettyNum(axTicks(2),big.mark = ","),las=2,cex.axis=0.75)
barplot(as.matrix(apply(dup_reads[c(2,3),],2,function(x) x/sum(x) )),
        las=2,ylim = c(0,1.2),yaxt='n',ylab="Proportion of PE reads",col=c("black","grey"))
axis(2,at=seq(0,1,0.2),las=1)

#par(mai=c(0,0,0,0))
#plot.new()
legend(-7,1.3,c("Non-duplicated","Duplicated"),bty="n",pch=15,ncol=3,col=c("black","grey"))

dev.off()  


