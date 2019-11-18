library(ggplot2)
library(GenomicFeatures)
library(zeallot)
library(tidyverse)
library(tximport)
library(tximportData)
library(edgeR)
library(SummarizedExperiment)
library(DESeq2)
library("RColorBrewer")
library("dplyr")
library(plyr)


# Load Utils
baseDir = "/home/luisas/Desktop/cluster/proj"
scripts_dir = file.path(baseDir,"code/ebola/src/scripts/")
source(paste0(scripts_dir,"de_utils.R"))
# Dirs
directory_zyagen = file.path(baseDir,"data/02_RNA-Seq/04_stringtie_counts/Zyagen")

# UMIS
## Create Summarized Experiment
ctabs<- iterate_files(directory_zyagen, 't_data.ctab')
coldata = extract_coldata(ctabs)
# get tx2gene table 
tx2gene <- read_tsv(ctabs[1])[, c("t_name","gene_id", "gene_name")]
rowdata = tx2gene
tx2gene = rowdata[,1:2]
# import tx 
txi <- tximport(ctabs, type = "stringtie", tx2gene = tx2gene, txOut = TRUE)
rownames(coldata) = coldata[["completeid"]]
rowdata = rowdata[match(rownames(txi[[1]]), as.character(rowdata[["t_name"]])),]
rownames(rowdata) = rowdata[["t_name"]]
#Create summarized experiment
se = SummarizedExperiment(assays = list(counts = txi[["counts"]],
                                        abundance = txi[["abundance"]],
                                        length = txi[["length"]]),
                          colData = coldata, 
                          rowData = rowdata)
# Create summarized experiment at gene level 
if (!is.null(tx2gene)){
  gi = summarizeToGene(txi, tx2gene = tx2gene)
  growdata = unique(rowdata[,2:3])
  growdata = growdata[match(rownames(gi[[1]]), growdata[["gene_id"]]),]
  rownames(growdata) = growdata[["tx"]]
  gse = SummarizedExperiment(assays = list(counts = gi[["counts"]],
                                           abundance = gi[["abundance"]],
                                           length = gi[["length"]]),
                             colData = DataFrame(coldata),
                             rowData = growdata)
}


se <- se[sapply(rownames(se), function(str) !str_detect(str,"MSTRG")),]
se <- se[sapply(rowData(se)$gene_id, function(str) !str_detect(str,"MSTRG")),]
se <- se[sapply(rowData(se)$t_name, function(str) !str_detect(str,"MSTRG")),]

#dge <- DGEList(counts = assays(se)$counts, group = se$dataset, genes = as.data.frame(mcols(se)))

umis = se[,se$method == "umis"]
mds = se[,se$method == "mds"]
filter = se[,se$method == "filter"]
calc_diff <- function(se_1, method){
  summary <- data.frame()
  for(tissue in unique(se_1$tissue)){
    se.filt <- se_1[,se_1$tissue == tissue]
    if(length(se.filt$sample) > 1 ){
      sample_1 <- se.filt[,se.filt$sample == 190507]
      sample_2 <- se.filt[,se.filt$sample == 160421]
      diff_1 <- abs(assays(sample_1)$counts - assays(sample_2)$counts)
      tissue_summary_1 <- data.frame("dif" =  diff_1, "tissue" = as.factor(tissue), "method" = method)
      summary <- rbind(summary, tissue_summary_1)
    }
  }
  return(summary)
}
umis_diff <- calc_diff(umis,"umis")
mds_diff <- calc_diff(mds,"mds")
filter_diff <- calc_diff(filter,"filter")
identical(umis_diff,mds_diff)

summary <- data.frame("umis" =  umis_diff$dif,"filter" =  filter_diff$dif, "mds" =  mds_diff$dif, "tissue" = as.factor(umis_diff$tissue), "method" = umis_diff$method)
k = 200
mask <- summary$umis > k & summary$mds > k & summary$filter > k
summary_filt <- summary[mask,]
wilcox.test(summary_filt$umis, summary_filt$filter)
wilcox.test(summary_filt$umis, summary_filt$mds)
wilcox.test(summary_filt$mds, summary_filt$filter)

u <- data.frame("dif" = summary_filt$umis, "tissue" = summary_filt$tissue, "method" = "umis")
m <- data.frame("dif" = summary_filt$mds, "tissue" = summary_filt$tissue, "method" = "mds")
f <- data.frame("dif" = summary_filt$filter, "tissue" = summary_filt$tissue, "method" = "filter")
s <- rbind(u,m,f)
s

p<-ggplot(s, aes(x=as.numeric(s$dif), color=s$method)) +
  geom_density(alpha=0.4)
p
boxplot(dif~tissue, data = s, col = s$method)

p<-ggplot(s, aes(x=tissue, y=as.numeric(s$dif), color=method)) +
  geom_boxplot()
p




#### ZYAGEN vs BATCH only on UMIS
directory_batch = file.path(baseDir,"data/02_RNA-Seq/04_stringtie_counts/Batch01")
ctabs_batch<- iterate_files(directory_batch, 't_data.ctab')
ctabs_batch_umis <-ctabs_batch[sapply(ctabs_batch, function(str) str_detect(str,"umis"))]
# get Summarized experiment object
se_batch_umis <-extract_se(ctabs_batch_umis)

ctabs_zyagen<- iterate_files(directory_zyagen, 't_data.ctab')
ctabs_zyagen_umis <-ctabs_zyagen[sapply(ctabs_zyagen, function(str) str_detect(str,"umis"))]
# get Summarized experiment object
se_zyagen_umis <-extract_se(ctabs_zyagen_umis)
se_zyagen_umis<- se_zyagen_umis[,se_zyagen_umis$tissue == "Liver"]
colData(se_zyagen_umis)
se_zyagen_umis <- se_zyagen_umis[sapply(rownames(se_zyagen_umis), function(str) !str_detect(str,"MSTRG")),]
se_batch_umis <- se_batch_umis[sapply(rownames(se_batch_umis), function(str) !str_detect(str,"MSTRG")),]
se_zyagen_umis <- se_zyagen_umis[sapply(rowData(se_zyagen_umis)$gene_id, function(str) !str_detect(str,"MSTRG")),]
se_batch_umis <- se_batch_umis[sapply(rowData(se_batch_umis)$gene_id, function(str) !str_detect(str,"MSTRG")),]
se_zyagen_umis <- se_zyagen_umis[sapply(rowData(se_zyagen_umis)$t_name, function(str) !str_detect(str,"MSTRG")),]
se_batch_umis <- se_batch_umis[sapply(rowData(se_batch_umis)$t_name, function(str) !str_detect(str,"MSTRG")),]

intersection <- intersect(rownames(se_zyagen_umis), rownames(se_batch_umis))
se_zyagen_umis <- se_zyagen_umis[sapply(rownames(se_zyagen_umis), function(str) str %in% intersection),]
se_batch_umis <- se_batch_umis[sapply(rownames(se_batch_umis), function(str) str %in% intersection),]
length(intersection)
length(se_zyagen_umis)
se_all <-cbind(se_zyagen_umis,se_batch_umis)
dge_all<- DGEList(counts = assays(se_all)$counts, group = se_all$dataset, genes = as.data.frame(mcols(se_all)))

#dge_zyagen <- DGEList(counts = assays(se_zyagen_umis)$counts, group = se_zyagen_umis$dataset, genes = as.data.frame(mcols(se_zyagen_umis)))
#dge_batch <- DGEList(counts = assays(se_batch_umis)$counts, group = se_batch_umis$dataset, genes = as.data.frame(mcols(se_batch_umis)))


# Calc logCPM for easing the calculation TMM 
dge_all <- calcNormFactors(dge_all) 
assays(se_all)$logCPM <- cpm(dge_all, log = TRUE,normalized.lib.sizes=TRUE, prior.count=0.25)
#filter out lowly expressed genes
mask <- rowMeans(assays(se_all)$logCPM) > 1
se_all <- se_all[mask,]
dge_all <- dge_all[mask,]

plotMDS(dge_all, top = 500)




# Calc logCPM for easing the calculation TMM 
dge <- DGEList(counts = assays(se)$counts, group = se$dataset, genes = as.data.frame(mcols(se)))
dge <- calcNormFactors(dge) 
assays(se)$logCPM <- cpm(dge, log = TRUE,normalized.lib.sizes=TRUE, prior.count=0.25)
#filter out lowly expressed genes
mask <- rowMeans(assays(se)$logCPM) > 1
se <- se[mask,]
dge <- dge[mask,]
colData(se)
#MDS plot based on TMM values
plotMDS(dge)

# UMIS
## Create Summarized Experiment
ctabs<- iterate_files(directory_zyagen, 't_data.ctab')
coldata = extract_coldata(ctabs)
# get tx2gene table 
tx2gene <- read_tsv(ctabs[1])[, c("t_name","gene_id", "gene_name")]
rowdata = tx2gene
tx2gene = rowdata[,1:2]
# import tx 
txi <- tximport(ctabs, type = "stringtie", tx2gene = tx2gene, txOut = TRUE)
rownames(coldata) = coldata[["completeid"]]
rowdata = rowdata[match(rownames(txi[[1]]), as.character(rowdata[["t_name"]])),]
rownames(rowdata) = rowdata[["t_name"]]
#Create summarized experiment
se = SummarizedExperiment(assays = list(counts = txi[["counts"]],
                                        abundance = txi[["abundance"]],
                                        length = txi[["length"]]),
                          colData = coldata, 
                          rowData = rowdata)
# Create summarized experiment at gene level 
if (!is.null(tx2gene)){
  gi = summarizeToGene(txi, tx2gene = tx2gene)
  growdata = unique(rowdata[,2:3])
  growdata = growdata[match(rownames(gi[[1]]), growdata[["gene_id"]]),]
  rownames(growdata) = growdata[["tx"]]
  gse = SummarizedExperiment(assays = list(counts = gi[["counts"]],
                                           abundance = gi[["abundance"]],
                                           length = gi[["length"]]),
                             colData = DataFrame(coldata),
                             rowData = growdata)
}


se <- se[sapply(rownames(se), function(str) !str_detect(str,"MSTRG")),]
se <- se[sapply(rowData(se)$gene_id, function(str) !str_detect(str,"MSTRG")),]
se <- se[sapply(rowData(se)$t_name, function(str) !str_detect(str,"MSTRG")),]

dge <- DGEList(counts = assays(se)$counts, group = se$dataset, genes = as.data.frame(mcols(se)))

umis = se[,se$method == "umis"]
mds = se[,se$method == "mds"]
filter = se[,se$method == "filter"]

calc_diff <- function(se_1, method){
  summary <- data.frame()
  for(tissue in unique(se_1$tissue)){
    se.filt <- se_1[,se_1$tissue == tissue]
    if(length(se.filt$sample) > 1 ){
      sample_1 <- se.filt[,se.filt$sample == 190507]
      sample_2 <- se.filt[,se.filt$sample == 160421]
      diff_1 <- abs(assays(sample_1)$counts - assays(sample_2)$counts)
      tissue_summary_1 <- data.frame("dif" =  diff_1, "tissue" = as.factor(tissue), "method" = method)
      summary <- rbind(summary, tissue_summary_1)
    }
  }
  return(summary)
}
umis_diff <- calc_diff(umis,"umis")
mds_diff <- calc_diff(mds,"mds")
filter_diff <- calc_diff(filter,"filter")






# -------------------
# Calc logCPM for easing the calculation TMM 
dge <- DGEList(counts = assays(se)$counts, group = se$dataset, genes = as.data.frame(mcols(se)))
dge <- calcNormFactors(dge) 
assays(se)$logCPM <- cpm(dge, log = TRUE,normalized.lib.sizes=TRUE, prior.count=0.25)
#filter out lowly expressed genes
mask <- rowMeans(assays(se)$logCPM) > 1
se <- se[mask,]
dge <- dge[mask,]

#MDS plot based on TMM values
plotMDS(dge)







#-----------------------------------------------------------------------------
# remove Lowly expresses genes - se filt 
cpmcutoff <- round(10/min(dge$sample$lib.size/1e+06), digits = 1)
sprintf("Cutoff: %s", cpmcutoff)
nsamplescutoff <- min(table(se$method))
mask <- rowSums(cpm(dge) > cpmcutoff) >= nsamplescutoff
se.filt <- se[mask, ]
dge.filt <- dge[mask, ]

avgexp <- rowMeans(assays(se)$logCPM)
par(mfrow = c(1, 1), mar = c(4, 5, 1, 1))
par(mar = c(4, 5, 1, 1))
h <- hist(avgexp, xlab = expression("Expression level (" * log[2] * "CPM)"), main = "",
          las = 1, col = "grey", cex.axis = 1.2, cex.lab = 1.5)
x <- cut(rowMeans(assays(se.filt)$logCPM), breaks = h$breaks)
lines(h$mids, table(x), type = "h", lwd = 10, lend = 1, col = "darkred")
legend("topright", c("All genes", "Filtered genes"), fill = c("grey", "darkred"))

# How many did we filter out? 
print("Filtered out")
nrow(se) - nrow(se.filt)
print("Genes left")
nrow(se.filt)


boxplot_summary <- rbind(dif_umis,dif_filter, dif_mds)

p<-ggplot(boxplot_summary, aes(x=dif, color=method)) +
  geom_density()+
  xlim(1,100)
p+ geom_vline(aes(xintercept=mean(dif)),
              color="blue", linetype="dashed", size=1)
mu <- ddply(boxplot_summary, "method", summarise, grp.mean=mean(as.numeric(dif)))
mu
# ---------------------------------------------

# removing outliers
Ndesv <- 5
mean(dif_umis$dif)
sd(dif_umis$dif)
mask1 <- dif_umis$dif > (mean(dif_umis$dif) + Ndesv * sd(dif_umis$dif))
Conf <- (1 - 1 / Ndesv^2)
Conf


k = 1
mask <- dif_umis$dif > k
mask2 <-dif_filter$dif > k
mask3 <-dif_mds$dif > k
mask <-mask & mask2 &mask3

dif_umis.filt <- dif_umis[mask,]
dif_filter.filt <- dif_filter[mask,]
dif_mds.filt <- dif_mds[mask,]

boxplot_summary.filt <- rbind(dif_umis.filt,dif_filter.filt, dif_mds.filt)
group_by(boxplot_summary.filt, method) %>%
  summarise(
    count = n(),
    median = median(dif, na.rm = TRUE),
    IQR = IQR(dif, na.rm = TRUE)
  )

# after removing 
library(plyr)
mu <- ddply(boxplot_summary.filt, "method", summarise, grp.mean=mean(as.numeric(dif)))
p<-ggplot(boxplot_summary.filt, aes(x=as.numeric(dif), color=method)) +
  geom_density(alpha=0.4)+
  xlim(0,20)

p













d <- density(as.numeric(dif_umis.filt$dif))
plot(d)
nrow(dif_umis)
nrow(dif_filter)
unique(boxplot_summary$method)
ggboxplot(boxplot_summary, y = "dif", 
          color = "method",
          ylab = "Weight", xlab = "Groups")


ggplot(as.data.frame(boxplot_summary), y=dif, color= method) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("Tissue") + 
  xlab("Difference") + 
  ggtitle("UMI vs MD")+ ylim(0,5)




library("dplyr")
group_by(boxplot_summary, method) %>%
  summarise(
    count = n(),
    median = median(dif, na.rm = TRUE),
    IQR = IQR(dif, na.rm = TRUE)
  )

d <- density(as.numeric(dif_umis$dif))
dif_umis$dif
plot(d)
dif_umis

library(ggplot2)
# Basic density
p <- ggplot(dif_umis, aes(x=dif)) + 
  geom_density() + xlim(0,100)
p


# Add mean line
p+ geom_vline(aes(xintercept=mean(dif)),
              color="blue", linetype="dashed", size=1)
