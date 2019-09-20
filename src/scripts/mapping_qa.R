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
library(zeallot)
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


iterate_files <- function(inpath, pattern_string){
  files <- list.files(path=inpath, pattern= pattern_string, full.names=TRUE, recursive=TRUE)
  return(files)
} 

get_info_sample <- function(file, pattern){
  info <- strsplit(strsplit(file,"/", fixed = T)[[1]][13],".", fixed = T)[[1]]
  tissue <- info[1]
  sample <- info[2]
  lane <- info[3]
  sample_name=paste(sample,tissue,lane,sep=".")
  
  res_path <- paste(inpath,tissue,sample,sep="/")
  res_file <- paste0(res_path,"/",paste(tissue,sample,lane,pattern,sep="."))
  return(c(tissue,sample,lane, sample_name, res_file))
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
  # Get mapped reads summary
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
summary
mapping_stats_df
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





# ----------------- BARPLOT INFERREDN STRANDNESS 
####### TO BE CHANGED
inpath_batch1="/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/Ebola_Raquel/hisat2/Batch_01"
pattern_infer = "inferred_strandness.txt"
infer_strandness_files <- iterate_files(inpath_batch1, pattern_infer)

# Extract info strandness
for(file in infer_strandness_files){
  # Extract tissue, sample and lane from filename of the hisat2 mapping reports
  c(tissue,sample,lane, sample_name, res_file)  %<-% get_info_sample(file,pattern_infer )
  # Error handling 
  if(! file.exists(res_file)){print(paste0("File does not exits ",file))}
  
  # PARSE INFERRED STRANDNESS FILE
  
  failed<- "Fraction of reads failed to determine:"
  stranded <- "Fraction of reads explained by \"1++,1--,2+-,2-+\""
  n_stranded <- "Fraction of reads explained by \"1+-,1-+,2++,2--\""
  strings_to_be_matched <- list(failed, stranded, n_stranded)
  con = file(res_file, "r")
  summary <- list()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    for( string_to_be_matched in strings_to_be_matched ){
      if (str_detect(line, string_to_be_matched )){
        summary[string_to_be_matched] <- as.numeric(strsplit(trimws(line, which ="left"),":")[[1]][1])
      }
    }
  }
  return(summary)
  close(con)
}
























