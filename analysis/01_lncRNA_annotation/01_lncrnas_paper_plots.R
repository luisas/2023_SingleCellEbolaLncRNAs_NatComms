library(rtracklayer);  library(stringr); library(ggplot2); library(grid); library(gridExtra); library(RColorBrewer); library(readr); library(matrixStats)
library(GenomicRanges); library(dplyr); library(zeallot); library(ggpubr); library(plyr); library(Gviz)


# Define paths for data
source(file.path("../utils/00_datapaths.R"))
# Import Utils
source(file.path("../utils/01_lncrna_annotation_utils.R"))

# Reoccurring paths
datadir <- file.path(data_path, "01_bulk_RNA-Seq_lncRNAs_annotation/")
dir_counts_ref <- file.path(data_path, "01_bulk_RNA-Seq_lncRNAs_annotation/04_quantification/")

# Human reference for comparison
lncRNAs_ref_human <- import(file.path(datadir, "01_PreliminaryFiles_rheMac10/Homo_sapiens.GRCh38.100_known_lncrna.gtf"))
mrna_ref_human <- import(file.path(datadir, "01_PreliminaryFiles_rheMac10/Homo_sapiens.GRCh38.100_known_proteincoding.gtf"))

# Macaque reference
lncRNAs_ref <- import(file.path(datadir, "01_PreliminaryFiles_rheMac10/Macaca_mulatta.Mmul_10.100_known_lncrna.gtf"))
mRNAs_ref <- import(file.path(datadir, "01_PreliminaryFiles_rheMac10/Macaca_mulatta.Mmul_10.100_known_proteincoding.gtf"))
mRNAs_ref <- mRNAs_ref[!is.na(mRNAs_ref$gene_biotype)]
mRNAs_ref <- mRNAs_ref[mRNAs_ref$gene_biotype == "protein_coding"]


all <- import(file.path(data_path,"01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf"))
# How many novel lncRNAs do i identify (Genes)
all_novel_lnc <- all[substr(all$gene_id,1,4) %in% c( "MSTR"),]
orthologs <- readRDS(file.path(data_path, "/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/orthologs_geneid.rds"))

# Define palettes 
palette <- c("#f72631", "#fc6d70","#febcbd", "#3153a3","#6f73b4")
palette_border <- c("#ff9933","#F9DF74", rep("black",4))
palette_extensive <- c(rep(palette[1],2), palette[seq(2,length(palette))])

palette_expression <-palette[c(1,2,4)]
palette_expression_extensive <-palette_extensive[c(1,2,3,5)]
palette_expression_border <-palette_border[c(1,2,3,5)]

# ------------------------
#    1. 


iterate_files <- function(inpath, pattern_string){
  files <- list.files(path=inpath, pattern= pattern_string, full.names=TRUE, recursive=TRUE)
  return(files)
}
abundance_files <- iterate_files(dir_counts_ref, "*.tsv")

# Read all files
abundances <- list()
for(i in 1:length(abundance_files)) {
  file <- readr::read_tsv(abundance_files[i])
  # Sum up values for double entries (https://github.com/gpertea/stringtie/issues/192)
  file <- file %>% dplyr::group_by(`Gene ID`)%>% dplyr::summarise(TPM = sum(TPM))
  abundances[[i]] <- data.frame(file$TPM, row.names =file$`Gene ID` )
}

# Summarize all TPMs from all quantification files 
rn <- rownames(abundances[[1]])
dat <- abundances[[1]]
for(i in 2:length(abundances)) {
  dat <- merge(dat, abundances[[i]],  by= "row.names", all.x= F, all.y= F) [,-1]
  rownames(dat) <- rn
}

expression <- dat
expression <- as.matrix(expression)



# ----------------------------------------------------------------------
# MEDIAN
# ----------------------------------------------------------------------
median_expression <- data.frame(id=rownames(expression), expr=log(rowMedians(as.matrix(expression))))
median_expression['type'] <- "0"
median_expression <- add_type(median_expression, all_novel_lnc$gene_id, "Novel lncRNAs")
median_expression <- add_type(median_expression, lncRNAs_ref$gene_id, "Annotated lncRNAs")
median_expression <- add_type(median_expression, mRNAs_ref$gene_id, "mRNAs")
median_expression <- median_expression[!median_expression$type == "0",]
expr_plot <- plot_expression(median_expression, palette_expression)+theme(plot.title = element_text(size = 22), axis.text = element_text( angle =0, vjust = 0.9),axis.title.y  = element_text(size = 22))+labs(y = "median expression (logTPM)")+theme(panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                       panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                                       axis.line = element_line(colour = "black"),
                                                                                                                                                                                                                                                       axis.ticks.length=unit(.2, "cm"),
                                                                                                                                                                                                                                                       axis.text= element_text(size = 15, color = "black"),
                                                                                                                                                                                                                                                       axis.title = element_text(size = 15),
                                                                                                                                                                                                                                                       plot.title = element_text(hjust = 0.5, size = 17),
                                                                                                                                                                                                                                                       axis.line.x=element_blank(),
                                                                                                                                                                                                                                                       axis.ticks.x=element_blank(),axis.text.x = element_blank(),
                                                                                                                                                                                                                                                       panel.background = element_rect(fill = "white"))

pdf(file.path(plots, "01/median_expr.pdf"), width = 5, height = 7)
expr_plot
dev.off()

# Test for significance 
median_expression_novel <- median_expression[median_expression$type == "Novel lncRNAs",]$expr
median_expression_annot <- median_expression[median_expression$type == "Annotated lncRNAs",]$expr
median_expression_pc <- median_expression[median_expression$type == "mRNAs",]$expr
wilcox.test(median_expression_novel, median_expression_pc)
wilcox.test(median_expression_novel, median_expression_annot)
wilcox.test(median_expression_pc, median_expression_annot)


# SUPPL : sperate intergenic and antisense
median_expression <- data.frame(id=rownames(expression), expr=log(rowMedians(expression)))
median_expression['type'] <- "0"
median_expression <- add_type(median_expression, intergenic_lnc$gene_id, "Intergenic Novel lncRNAs")
median_expression <- add_type(median_expression, antisense_lnc$gene_id, "Antisense Novel lncRNAs")
median_expression <- add_type(median_expression, lncRNAs_ref$gene_id, "Annotated lncRNAs")
median_expression <- add_type(median_expression, mRNAs_ref$gene_id, "mRNAs")
median_expression <- median_expression[!median_expression$type == "0",]

levels <-  c("Intergenic Novel lncRNAs", "Antisense Novel lncRNAs","Annotated lncRNAs", "mRNAs")
pe_sep <- plot_expression(median_expression, palette_expression_extensive, level = levels, title = "Median Expression", palette_expression_border,sizes = c(1,1,0.5,0.5))+theme(plot.title = element_text(size = 22), axis.text = element_text( angle =0, vjust = 0.9),axis.title.y  = element_text(size = 15))+labs(y = "median expression (logTPM)")+theme(panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                             panel.grid.minor = element_blank(),                                                                                                                                                                                                                                                                                                                                                panel.background = element_rect(fill = "white"))

pdf(file.path(plots, "01/SUPPL_median_expr.pdf"), width = 5, height = 7)
pe_sep[[1]]
dev.off()





