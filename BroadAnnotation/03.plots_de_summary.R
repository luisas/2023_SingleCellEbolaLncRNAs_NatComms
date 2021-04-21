#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(edgeR))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library("wesanderson"))
suppressMessages(library(org.Mmu.eg.db))
suppressMessages(library(rtracklayer))
options(stringsAsFactors = F)


# 0. Palettes & themes 
pal1 <- wes_palette("Darjeeling1", 5, type = "discrete")
pal2 <- wes_palette("Darjeeling2", 5, type = "discrete")
colfunc <- colorRampPalette(c("grey", "red"))

theme_paper <- theme(legend.title = element_blank())+theme(panel.background = element_rect(fill = "white", colour = "white"))+theme(panel.background = element_rect(fill = "white", colour = "grey50"))+theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text = element_text(size = 8))


# 1. Input Variables and Files
tissue <- "Lymph node"
outpath <- paste0("/home/luisas/Desktop/cluster/data/99_BroadAnnotation_Feb2021/05_DEA/Tissues/",tissue,"/")

# Sade objects
my_data <- readRDS(paste0(outpath,tissue,".voom_limma.data_objects.rds"))

dge <- my_data[["dge"]] 
v <- my_data[["v"]] 
fit <- my_data[["fit"]] 
efit <- my_data[["efit"]] 

FDRTHRESHOLD <- 0.05


# 2. Read in metadata
metadata <-  readRDS("/home/luisas/Desktop/cluster/data/99_BroadAnnotation_Feb2021/metadata_full.rds")
# Only keep samples for which i have both counts and metadata
metadata_tissue <- metadata[metadata$tissue %in% tissue, ]

# Only retain samples with more than 3M reads
metadata_tissue <- metadata_tissue[(metadata_tissue$hostReadCount > 3* 10^6),]


ref_ensembl <- import("/home/luisas/Desktop/cluster/gene_annotations/ensembl_release100/rheMac10/Macaca_mulatta.Mmul_10.100.gtf")

get_names <- function(ids, ref_e  = ref_ensembl){
  ids <- unlist(lapply(ids, function(x) strsplit(x,".", fixed = T)[[1]][1]))
  names <- ref_e[ref_e$gene_id %in% ids, ]$gene_name
  names <- names[!is.na(names)]
  return(unique(names)) 
}
do_go <- function(list, universe, keytype="SYMBOL", title = ""){
  ego <- clusterProfiler::enrichGO(gene = get_names(list),OrgDb =  org.Mmu.eg.db, keyType = keytype,ont = "BP",universe =get_names(universe))
  if(!is.null(ego)){
    return(clusterProfiler::dotplot(ego, showCategory = 5,  font.size = 4)+ggtitle(title)+theme_paper)
  }else{
    return(NA)
  }
}





plots_path <- file.path(outpath, "plots")
dir.create(plots_path, recursive = T,  showWarnings = F)
pcas_path <- file.path(plots_path, "de_summary")


# -----------------------------------------------
#       Plot summary of samples and viral reads
# -----------------------------------------------

var <- "viralloadPerSasmple"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
dir.create(pcas_path, recursive = T,  showWarnings = F)
df <- metadata_tissue[order(metadata_tissue$dpi_time_factor),]
df$X <- factor(df$X,levels = df$X)
variable <- "id.cohort"
ggplot(df, aes(x = X, y = pc.viral.reads, col  =  eval(parse(text = variable)), fill = eval(parse(text = variable))))+geom_bar(stat = "identity")+theme_paper+xlab("Samples sorted by DPI")+ylab("% Viral read per sample")+theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+ggtitle(tissue)+scale_color_brewer(palette = "Set3")+scale_fill_brewer(palette = "Set3")
dev.off()

# -----------------------------------------------
#       Plot summary INPUT GENES
# -----------------------------------------------


exprs_genes <- (rownames(dge))
ref <- import("/home/luisas/Desktop/cluster/data/99_BroadAnnotation_Feb2021/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf")
matching <- read.table("/home/luisas/Desktop/slncky-master/annotations/ensemblSource.txt")
lnc_ids <- matching[matching$V2 == "lncRNA",]$V1
pc_ids <- matching[matching$V2 == "protein_coding",]$V1
lnc_gene_ids <- ref[ref$type == "transcript" & ref$transcript_id %in% lnc_ids, ]$gene_id
pc_gene_ids <- ref[ref$type == "transcript" & ref$transcript_id %in% pc_ids, ]$gene_id


summarize_genes <- function(exprs_genes, lnc_ids = lnc_gene_ids, pc_ids = pc_gene_ids, nopc = F ){
  genes_type <- data.frame(gene_id = exprs_genes) 
  genes_type$type <- "Other"
  genes_type[genes_type$gene_id %in% pc_ids, ]$type <- "Protein Coding"
  genes_type[genes_type$gene_id %in% lnc_ids, ]$type <- "LncRNAs - Annotated"
  genes_type[grepl("MSTR",genes_type$gene_id), ]$type <- "LncRNAs -Novel"
  if(nopc == T){
    genes_type <- genes_type[!(genes_type$type %in% c("Protein Coding", "Other")), ]
  }
  
  p <- ggplot(genes_type, aes(x = type, fill = type))+geom_bar(stat= "count", position = "dodge")+scale_fill_manual(values = pal1)+ggtitle(paste("Total number of genes ", nrow(genes_type), sep = ": "))+xlab("")+theme(axis.text.x=element_text(angle = 28, hjust = 1))+ylab("")+theme_minimal()
  return(list(genes_type,p))
}

l <- summarize_genes(exprs_genes)
genes_type <- l[[1]]
l[[2]]



var <- "InputGenes_summary_noPC"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
dir.create(pcas_path, recursive = T,  showWarnings = F)
summarize_genes(exprs_genes,nopc = T)[[2]]+theme_paper
dev.off()

var <- "InputGenes_summary_PC"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
dir.create(pcas_path, recursive = T,  showWarnings = F)
summarize_genes(exprs_genes,nopc = F)[[2]]+theme_paper
dev.off()

# -----------------------------------------------
#       Plot summary DE general 
# -----------------------------------------------
var <- "DE_general"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
# General Summary: number of DE
n_de_summary <- data.frame(summary(decideTests(efit)))
names(n_de_summary) <- c("direction", "comparison", "n")
n_de_summary <- n_de_summary[n_de_summary$direction != "NotSig",]
n_de_summary$comparison <- factor(gsub("vsBaseline", "", n_de_summary$comparison) , levels = c("Early", "Middle", "Late"))
ggplot(n_de_summary, aes( x= comparison, fill = direction, y = n ) )+ geom_bar(stat = "identity", alpha = 0.7)+theme_paper+xlab("")+ylab("")+scale_fill_manual(values = rev(pal2[1:2]))+ggtitle("Number of DE genes")
dev.off()


# -----------------------------------------------
#       Plot summary DE lnc 
# -----------------------------------------------
res_late <- topTable(efit, coef = "LatevsBaseline", adjust.method = "fdr", confint = TRUE, n = Inf)
sig_de_late <- res_late[res_late$adj.P.Val < FDRTHRESHOLD, ]

res_middle <- topTable(efit, coef = "MiddlevsBaseline", adjust.method = "fdr", confint = TRUE, n = Inf)
sig_de_middle <- res_middle[res_middle$adj.P.Val < FDRTHRESHOLD, ]

res_early <- topTable(efit, coef = "EarlyvsBaseline", adjust.method = "fdr", confint = TRUE, n = Inf)
sig_de_early <- res_early[res_early$adj.P.Val < FDRTHRESHOLD, ]


# Numbers 
df_late <- summarize_genes(rownames(sig_de_late), nopc = F )[[1]]
df_late$stage <- "Late"
df_middle <- summarize_genes(rownames(sig_de_middle), nopc = F )[[1]]
df_middle$stage <- "Middle"
df_early <- summarize_genes(rownames(sig_de_early), nopc = F )[[1]]
df_early$stage <- "Early"
df_all <- rbind(df_late, df_middle, df_early)
df_all$stage <- factor(df_all$stage, levels = c("Early", "Middle", "Late"))

var <- "DE_lnc"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
ggplot(df_all[!(df_all$type %in% c("Protein Coding", "Other")), ], aes(x =stage, fill = type ))+geom_bar(stat = "count")+theme_minimal()+scale_fill_manual(values = pal1)+ggtitle(paste("Total ", nrow(df_all[!(df_all$type %in% c("Protein Coding", "Other")), ]), sep = ": "))
dev.off()


# -----------------------------------------------
#       Enrichments
# -----------------------------------------------

var <- "Up_Late"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
do_go(rownames(sig_de_late[sig_de_late$logFC > 0,]),rownames(res_late), title = var)
dev.off()

var <- "Down_Late"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
do_go(rownames(sig_de_late[sig_de_late$logFC < 0,]),rownames(res_late), title = var)
dev.off()


var <- "Up_Midlle"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
do_go(rownames(sig_de_middle[sig_de_middle$logFC > 0,]),rownames(res_middle), title = var)
dev.off()

var <- "Down_Middle"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
do_go(rownames(sig_de_middle[sig_de_middle$logFC < 0,]),rownames(res_middle), title= var)
dev.off()


var <- "Up_Early"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
do_go(rownames(sig_de_early[sig_de_early$logFC > 0,]),rownames(res_early), title = var)
dev.off()

var <- "Down_Early"
jpeg(file=file.path(pcas_path, paste(tissue,var, ".jpeg", sep = "_")), width=3.25,height=3.25,units="in", res = 600)
do_go(rownames(sig_de_early[sig_de_early$logFC < 0,]),rownames(res_early), title= var)
dev.off()







