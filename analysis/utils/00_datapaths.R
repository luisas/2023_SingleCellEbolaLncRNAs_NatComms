
var = "localtorre"

if(var == "localtorre"){
  # Luisa: Local Torre 
  clusterpath <- "/home/luisa/Desktop/cluster/"
  gene_annotation_path <- file.path(clusterpath, "gene_annotation")
  data_path <- file.path(clusterpath, "data")
  plots <- file.path(data_path,"plots")
  
}else if(var == "locallaptop"){
  # Luisa: Local Laptop
  clusterpath <- "/Users/luisasantus/Desktop/mn_cluster/mount_dirs"
  gene_annotation_path <- file.path(clusterpath, "gene_annotation")
  data_path <- file.path(clusterpath, "Data")
  plots <- file.path(data_path,"plots")
  
}else if(var == "gpfs"){
  # Cluster: GPFS 
  clusterpath <- "/gpfs/projects/bsc83/Data/"
  gene_annotation_path <- file.path(clusterpath, "gene_annotation")
  data_path <- file.path(clusterpath, "Ebola")
  plots <- file.path(data_path,"plots")
}




