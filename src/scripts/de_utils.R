#  ------------------  Utils 
library(pheatmap)

iterate_files <- function(inpath, pattern_string){
  files <- list.files(path=inpath, pattern= pattern_string, full.names=TRUE, recursive=TRUE)
  return(files)
} 

extract_info_sample_from_filename <- function(file){
  file_no_ext <- strsplit(file,".", fixed = T)[[1]][1]
  info <- strsplit(tail(strsplit(file_no_ext,"/", fixed = T)[[1]],1),"_", fixed = T)[[1]]
  dataset <- info[1]
  tissue <- info[2]
  dpo <- info[3]
  sample <- info[4]
  sample_name=paste(dataset,tissue,dpo,sample,sep="_")
  return(c(dataset,tissue,dpo,sample, sample_name))
}
get_filename <- function(file){
  return(tail(strsplit(file,"/", fixed = T)[[1]],1))
}
# Only for htseq counts
get_filename_and_dirs <- function(file){
  infos %<-% extract_info_sample_from_filename(file)
  file_path = paste(infos[c(1:4)],collapse="/")
  filename <- get_filename(file)
  file_path = paste(file_path, "htseq_counts", sep ="/")
  file_path = paste(file_path, filename, sep ="/")
  return(file_path)
}


create_dds <- function(sampleFiles){
  samples_info <- sapply(sampleFiles,extract_info_sample_from_filename)
  # Create ddsHTSeq Object 
  sampleTable <- data.frame(sampleName = sampleFiles,
                            fileName = sapply(sampleFiles,get_filename_and_dirs),
                            condition = samples_info[1,],
                            tissue = samples_info[2,], 
                            dpo = samples_info[3,], 
                            id = samples_info[4,],
                            complete_id = samples_info[5,]
  )
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                    directory = directory,
                                    design= ~ 1)
  return(dds)
}


add_gene_sizes <- function(dds,gene_lengths){
  gene_lengths<-gene_lengths[match(rownames(mcols(dds, use.names=TRUE)), gene_lengths$gene_id),]
  #Checking
  all(rownames(dds) == gene_lengths$gene_id)
  sum(is.na(gene_lengths$size))
  mcols(dds) <- cbind(mcols(dds), gene_lengths)
  return(dds)
}

plot_counts_heatmap <- function(dds, outpath,method = "pearson"){
  df <- as.data.frame(assays(dds)$counts)
  colnames(df) <- dds$complete_id
  png(outpath,width = 1400,height = 500)
  pheatmap(cor(df,method=method))
  dev.off()
}