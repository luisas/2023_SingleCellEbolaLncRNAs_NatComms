#source("/gpfs/projects/bsc83/Projects/Ebola/code/ebola/src/scripts/WGCNA_functions.R")
source("/home/luisas/Desktop/RUN-WGCNA-master/WGCNA_functions.R")
result_path <- file.path("/home/luisas/Desktop/cluster/data/RObjects/results")

library(reshape2)
library(ggplot2)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

immune.combined <- readRDS(file.path(result_path,"seurat_pbmc_rhemac10_merged_aftercellandgeneqc_afterScrublet.rds"))
count_matrix <- immune.combined@assays$RNA@counts
dat <- as.matrix(count_matrix)
determineSoftPowerWGCNA(data1=dat, outFile="powerPlots.png", propGenes=.1)
net <- runWGCNA(data1=dat, propGenes=.1, softPower=3, signedNetwork=TRUE)
plotModulesCut(referenceDataset=net, outFile="modClusterPlot.pdf")
e1 <- calculateModuleEigengenes(referenceDataset=net, split=2)

# Select the groups
groups <- gsub('-a.*dge.*$', '', x=rownames(e1))
groups <- t(data.frame(strsplit(groups, split="-", fixed=TRUE)))
groups <-groups[,c(2,5)]
groups[ ,2] <- gsub('H00', '', groups[ ,2])
rownames(groups) <- NULL
colnames(groups) <- c("Cond", "Hour")
e1 <- cbind(e1, groups)
# Melt
melted <- melt(data=e1, id.vars=c("Cond", "Hour"))
melted <- melted[order(as.numeric(melted$Hour)), ]
melted$variable <- gsub("ME", "", melted$variable)
# Plot
png(file="eigengenBoxplot.png")
q <- qplot(data=melted, y=value, x=Oxygen, facets=Media~variable, 
           geom=c("boxplot", "point"), ylab="Eignegene expression",
           colour=Oxygen)
q + scale_x_discrete(limits=c("1","5","20"))
dev.off()

# Find Hub Genes
hubs <- moduleHubGenes(referenceDataset=net, MEs=calculateModuleEigengenes(referenceDataset=net, split=2), nGenes=10, split=2)
hubs[ ,colnames(hubs) == "blue"]