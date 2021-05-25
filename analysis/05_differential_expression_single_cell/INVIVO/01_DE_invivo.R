knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
library(Seurat)
library(wesanderson)
library(Matrix)
library(SingleCellExperiment)
library(stringr)
library(rtracklayer)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(circlize)
library(reshape2)

theme_paper <- theme(legend.title = element_blank())+theme(panel.background = element_rect(fill = "white", colour = "white"))+theme(panel.background = element_rect(fill = "white", colour = "grey50"))+theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 18))


# Imports
source("../../utils/02_sc_utils.R")
source("../../utils/03_de_utils.R")
source("../../utils/04_utils_graph.R")

# Gene annotation
ref <- import(file.path("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf"))
ebola_ref <- import(file.path("/home/luisas/Desktop/cluster/data/00_RawData/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.gtf"))

# Output paths
robjectsdir <- file.path("/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/")

# LncRNAs
robjectsdir_stats <- "/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/"
all_lncrnas <- readRDS(file.path(robjectsdir_stats, "all_lncrnas.rds"))
annotated_mrnas <- readRDS(file.path(robjectsdir_stats,"annotated_mrnas.rds"))

# Seurat object
immune.combined <- readRDS("/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds")

# Ortholgs
orthologs <- readRDS("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/orthologs_geneid.rds")

# Check LncPedia
lncpedia_immune <- import("../../lncpedia_immune.gtf")
lncpedia_infection <- import("../../lncpedia_infection.gtf")

lncpedia <- unique(c(unique(unlist(lncpedia_immune@elementMetadata[, grepl("gene", colnames(lncpedia_immune@elementMetadata))])), unique(unlist(lncpedia_infection@elementMetadata[, grepl("gene", colnames(lncpedia_infection@elementMetadata))]))))

immlnc <- read.table("/home/luisas/Desktop/cluster/data/00_RawData/ImmLnc/Lnc_Pathways_Sig.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
reported_immlnc <- immlnc$lncRNA_symbol
# -------------------------------------------
#    Define threshold for the whole analysis
# -------------------------------------------

# Corrected q-value .05
FDRTHRESHOLD<- 0.05
# Fold-Change of 30%
FCTHRESHOLD <- log(1.23)

orthologs <- readRDS("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/orthologs_geneid.rds")
orthologs <- distinct(orthologs)
orthologs_red <- distinct(orthologs[,c("gene_id","lnc", "orthologGeneSymbol","alignScore")])
orthologs_max <- orthologs_red %>% dplyr::group_by(orthologGeneSymbol) %>% dplyr::summarise(alignScore=max(alignScore))
orthologs_max$key <- paste(orthologs_max$orthologGeneSymbol, orthologs_max$alignScore, sep ="_")
orthologs_red$key <- paste(orthologs_red$orthologGeneSymbol, orthologs_red$alignScore, sep ="_")
orthologs$key <- paste(orthologs$orthologGeneSymbol, orthologs$alignScore, sep ="_")
orthologs <- orthologs_max %>% left_join(orthologs, by = c("key"))
colnames(orthologs) <- gsub("orthologGeneSymbol.x","orthologGeneSymbol",colnames(orthologs))
orthologs <- distinct(orthologs)

files <- iterate_files(robjectsdir, "MAST_invivo_model_*")
get_celltype <- function(x){strsplit(basename(x), "_")[[1]][4]}
get_stage <- function(x){strsplit(strsplit(basename(x), "_")[[1]][5], '.', fixed = T )[[1]][1]}


# Read Rds files resulting from DE Analysis and assign to correct variable name
mast_objects_celltype <- unlist(lapply(files, function(name) {  appendix <- paste(get_celltype(name), get_stage(name),  sep="_")
                                varname=paste("mast", appendix,  sep="_") ;
                                assign(varname, readRDS(name), envir = .GlobalEnv);
                                print(varname);
                                # Assign right comparison depending on variable name
                                eval(parse(text=paste0(varname,"$comparison <- \"", appendix, "\"")), envir = .GlobalEnv);
                                eval(parse(text=paste0(varname,"$celltype <- \"", get_celltype(name), "\"")), envir = .GlobalEnv);
                                eval(parse(text=paste0(varname,"$stage <- \"", get_stage(name), "\"")), envir = .GlobalEnv);
                                return(varname)
                              }))

# Extract all celltypes for which some DE files were found
celltypes<- unique(unlist(lapply(mast_objects_celltype, function(x) str_split(x, "_")[[1]][2])))

# Create an object per celltype ( rbind )
mast_objects_celltype_all <- unlist(lapply(celltypes, function(celltype){
                                    varnames_selected <-mast_objects_celltype[unlist(lapply(mast_objects_celltype,function(x) (grepl( celltype,x))))];
                                    ex=paste("mast_",celltype, "_all <-rbind(",paste0(varnames_selected, collapse = ","), ")", sep="");
                                    eval(parse(text=ex), envir = .GlobalEnv);
                                    return(paste("mast_",celltype, "_all", sep = ""));
                                    }))

# Print names of the variables
# mast_objects_celltype
# mast_objects_celltype_all
# mast_B_early

# -------------------------------------------------------------------------
# Select genes
# -------------------------------------------------------------------------

get_de <- function(df_filtered,  subset = all_lncrnas, fdrthreshold = FDRTHRESHOLD, fc = FCTHRESHOLD){
  # Only select lncrnas
  df_filtered <- df_filtered[df_filtered$primerid %in% subset,]
  # Only select significant ones
  df_filtered <- df_filtered[fdr<fdrthreshold &abs(coef) > fc,]
  return(df_filtered)
}


# Extract significant lnc and mrna from MAST tables for each celltype and stage
lnc_de_objects <- unlist(lapply(mast_objects_celltype, function(x){
                                          varname <- paste("lnc_", gsub("mast_", "", x), sep = "");
                                          eval(parse(text=paste0(varname," <- get_de(", x, " ,subset = all_lncrnas )")), envir = .GlobalEnv);
                                          return(varname)
                                      }))

mrna_de_objects <- unlist(lapply(mast_objects_celltype, function(x){
                                          varname <- paste("mrna_", gsub("mast_", "", x), sep = "");
                                          eval(parse(text=paste0(varname," <- get_de(", x, " ,subset = annotated_mrnas )")), envir = .GlobalEnv);
                                          return(varname)
                                      }))


# Create one table of significant genes per celltype (Merging stages)
de_celltype_lnc <- unlist(lapply(celltypes, function(celltype){
                                    varnames_selected <- lnc_de_objects[unlist(lapply(lnc_de_objects,function(x) (grepl( celltype,x))))];
                                    ex=paste("lnc_",celltype, "_all <-rbind(",paste0(varnames_selected, collapse = ","), ")", sep="");
                                    eval(parse(text=ex), envir = .GlobalEnv);
                                    return(paste("lnc_",celltype, "_all", sep = ""));
                                    }))

de_celltype_mrna <- unlist(lapply(celltypes, function(celltype){
                                    varnames_selected <- mrna_de_objects[unlist(lapply(mrna_de_objects,function(x) (grepl( celltype,x))))];
                                    ex=paste("mrna_",celltype, "_all <-rbind(",paste0(varnames_selected, collapse = ","), ")", sep="");
                                    eval(parse(text=ex), envir = .GlobalEnv);
                                    return(paste("mrna_",celltype, "_all", sep = ""));
                                    }))


# Lists with gene names only - prepare for heatmap
# ---------------------------------------------------------------------------------------------------------------------------------------------
de_lists_lnc_my <- list(lnc_Monocyte_late$primerid, lnc_Monocyte_middle$primerid, lnc_Monocyte_early$primerid )
names(de_lists_lnc_my) <- c("lnc_Monocyte_late", "lnc_Monocyte_middle", "lnc_Monocyte_early")



#de_lists_lnc_neut <- list(lnc_Neutrophil_late$primerid, lnc_Neutrophil_middle$primerid, lnc_Neutrophil_early$primerid )
#names(de_lists_lnc_neut) <- c("lnc_Neutrophil_late", "lnc_Neutrophil_middle", "lnc_Neutrophil_early")

de_lists_lnc_T <- list(lnc_T_late$primerid, lnc_T_middle$primerid, lnc_T_early$primerid )
names(de_lists_lnc_T) <- c("lnc_T_late", "lnc_T_middle", "lnc_T_early")
de_lists_lnc_T <- de_lists_lnc_T[unlist(lapply(1:length(de_lists_lnc_T), function(index) length(de_lists_lnc_T[[index]]) != 0))]

de_lists_lnc_B <- list(lnc_B_late$primerid, lnc_B_middle$primerid, lnc_B_early$primerid )
names(de_lists_lnc_B) <- c("lnc_B_late", "lnc_B_middle", "lnc_B_early")
de_lists_lnc_B <- de_lists_lnc_B[unlist(lapply(1:length(de_lists_lnc_B), function(index) length(de_lists_lnc_B[[index]]) != 0))]

de_lnc_all <- c(de_lists_lnc_my, de_lists_lnc_T, de_lists_lnc_B)


de_lists_mrna_my <- list(mrna_Monocyte_late$primerid, mrna_Monocyte_middle$primerid, mrna_Monocyte_early$primerid )
names(de_lists_mrna_my) <- c("mrna_Monocyte_late", "mrna_Monocyte_middle", "mrna_Monocyte_early")

#de_lists_mrna_neu <- list(mrna_Neutrophil_late$primerid, mrna_Neutrophil_middle$primerid, mrna_Neutrophil_early$primerid )
#names(de_lists_mrna_my) <- c("mrna_Neutrophil_late", "mrna_Neutrophil_middle", "mrna_Neutrophil_early")

de_lists_mrna_T <- list(mrna_T_late$primerid, mrna_T_middle$primerid, mrna_T_early$primerid )
names(de_lists_mrna_T) <- c("mrna_T_late", "mrna_T_middle", "mrna_T_early")
de_lists_mrna_T <- de_lists_mrna_T[unlist(lapply(1:length(de_lists_mrna_T), function(index) length(de_lists_mrna_T[[index]]) != 0))]

de_lists_mrna_B <- list(mrna_B_late$primerid, mrna_B_middle$primerid, mrna_B_early$primerid )
names(de_lists_mrna_B) <- c("mrna_B_late", "mrna_B_middle", "mrna_B_early")
de_lists_mrna_B <- de_lists_mrna_B[unlist(lapply(1:length(de_lists_mrna_B), function(index) length(de_lists_mrna_B[[index]]) != 0))]

de_mrna_all <- c(de_lists_mrna_my, de_lists_mrna_T, de_lists_mrna_B)

# ---------------------------------------------------------------------------------------------------------------------------------------------


de_lists_lnc <- list(lnc_Monocyte_late$primerid,lnc_Monocyte_middle$primerid, lnc_Monocyte_early$primerid,lnc_T_late$primerid, lnc_T_middle$primerid, lnc_T_early$primerid, lnc_B_late$primerid, lnc_B_middle$primerid, lnc_B_early$primerid)
names(de_lists_lnc) <- c("lnc_Monocyte_late", "lnc_Monocyte_middle", "lnc_Monocyte_early","lnc_T_late", "lnc_T_middle", "lnc_T_early", "lnc_B_late", "lnc_B_middle", "lnc_B_early")


de_lists_mrna <- list(mrna_Monocyte_late$primerid,mrna_Monocyte_middle$primerid, mrna_Monocyte_early$primerid,mrna_T_late$primerid, mrna_T_middle$primerid, mrna_T_early$primerid, mrna_B_late$primerid, mrna_B_middle$primerid, mrna_B_early$primerid )
names(de_lists_mrna) <- c("mrna_Myeloid_late", "mrna_Myeloid_middle", "mrna_Myeloid_early","mrna_T_late", "mrna_T_middle", "mrna_T_early", "mrna_B_late", "mrna_B_middle", "mrna_B_early")

de_all_genes <- unique(unlist(c(unlist(de_lists_mrna), unique(de_lists_lnc))))
saveRDS(de_all_genes, file.path(robjectsdir, "de_all_genes.rds"))
saveRDS(de_lnc_all, file.path(robjectsdir, "de_lnc.rds"))


lnc_de_objects
mrna_de_objects




eval(parse(text= paste("lnc_de_all <- rbind(",paste0(de_celltype_lnc, collapse = ","),")")))
lnc_de_all <- lnc_de_all[lnc_de_all$celltype != "Neutrophil",]
lnc_de_all$type <- "lnc"
lnc_de_all$subtype <- ifelse(substr(lnc_de_all$primerid,1,5) == "MSTRG", "novel", "annotated")
eval(parse(text= paste("mrna_de_all <- rbind(",paste0(de_celltype_mrna, collapse = ","),")")))
mrna_de_all <- mrna_de_all[mrna_de_all$celltype != "Neutrophil",]
mrna_de_all$type <- "pc"
# Create df from stats fro seeing directinality
mrna_de_all$subtype <- "annotated"
de_all <- rbind(lnc_de_all, mrna_de_all)

# Remove DE in neutrophils
de_all <- de_all[de_all$celltype != "Neutrophil",]

de_all$direction <- ifelse(de_all$coef < 0 , "down", "up")
de_all$direction <- factor(de_all$direction, c("down", "up"))
de_all$stage <- factor(de_all$stage, c("early", "middle", "late"))


saveRDS(de_all, "/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_all_stages.rds")
saveRDS(unique(unlist(de_lists_lnc)), "/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/de_lnc.rds")


set.seed(123)

stats_complete <- rbind( mast_Monocyte_all,mast_B_all,mast_T_all  )

# Convert into matrix
m <-xtabs(coef~primerid+comparison, stats_complete)
attr(m, "class") <- NULL

m <- m[,c("Monocyte_early" ,"Monocyte_middle","Monocyte_late", "B_early" ,"B_middle","B_late", "T_early" ,"T_middle","T_late")]

m_significance<-xtabs(fdr~primerid+comparison, stats_complete)
attr(m_significance, "class") <- NULL
m_significance <- m_significance[,c("Monocyte_early" ,"Monocyte_middle","Monocyte_late", "B_early" ,"B_middle","B_late", "T_early" ,"T_middle","T_late")]

saveRDS(stats_complete, "/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/stats_complete.rds")
saveRDS(m, "/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/m.rds")
m <- readRDS("/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/04_DE/m.rds")

# Neutrophils nly appear later on
pal_celltypes <- brewer.pal(4, "Set2")
pal_celltypes <-wes_palette("GrandBudapest1", 4)
pal_celltypes[3] <- "#AE4E4E"
pal_celltypes <- pal_celltypes[c(2,1,3,4)]
immune.combined$group <- factor(immune.combined$group, levels = c("baseline", "early", "middle", "late"))
neut <- subset(immune.combined, ident = "Neutrophil")
DimPlot(immune.combined, split.by = "group", cols = pal_celltypes)+theme_sc
DimPlot(neut, pt.size = 1, group.by = "group", cols = rev(c(brewer.pal(3, "Set1"))))+ggtitle("Neutrophils colored by infection stage")+ylim(c(4,12))+theme_sc



# -------------------------------------------------------------------------
# Create Heatmap(s)
# -------------------------------------------------------------------------
de_lists_lnc <- de_lists_lnc[unlist(lapply(1:length(de_lists_lnc), function(index) length(de_lists_lnc[[index]]) != 0))]
h1_lnc<- hm(m= m,row_title = names(de_lists_lnc)[1], features  = de_lists_lnc[1][[1]], first = TRUE, row_names = FALSE)
hms <- lapply(2:length(de_lists_lnc), function(index) hm(m,row_title = names(de_lists_lnc)[index], features  = de_lists_lnc[index][[1]], row_names = FALSE))

length(hms)
heatmap_lnc <- h1_lnc %v% hms[[1]]%v%hms[[2]]%v%hms[[3]]%v%hms[[4]]%v%hms[[5]]%v%hms[[6]]
heatmap_lnc

# mrna ---------
de_lists_mrna <- de_lists_mrna[unlist(lapply(1:length(de_lists_mrna), function(index) length(de_lists_mrna[[index]]) != 0))]
h1 <- hm(m= m,row_title = names(de_lists_mrna)[1], features  = de_lists_mrna[1][[1]], first = TRUE, row_names = FALSE)
hms_mrna <- lapply(2:length(de_lists_mrna), function(index) hm(m,row_title = names(de_lists_mrna)[index], features  = de_lists_mrna[index][[1]], row_names = FALSE))

length(de_lists_mrna)

heatmap_mrna <- h1 %v% hms_mrna[[1]]%v%hms_mrna[[2]]%v%hms_mrna[[3]]%v%hms_mrna[[4]]%v%hms_mrna[[5]]%v%hms_mrna[[6]]%v%hms_mrna[[7]]%v%hms_mrna[[8]]
heatmap_mrna

de_stats_direction_plot <- de_all %>% group_by(primerid, celltype,stage,type) %>% dplyr::summarise(direction = toString((direction)))
de_stats_direction_plot <- de_all %>% group_by(primerid, celltype,stage,type,subtype) %>% dplyr::summarise(direction = toString((direction)))
de_stats_direction_plot$stage_short <- factor(toupper(substr(de_stats_direction_plot$stage,1,1)), levels=c("E", "M", "L"))

de_stats_direction_plot$direction <- factor(de_stats_direction_plot$direction, c("down", "up"))
# In general, agnostic of celltype and stage, how many are up and how many are down




ggplot(de_stats_direction_plot[de_stats_direction_plot$type == "pc",], aes(x = stage_short, fill = direction))+geom_bar(stat="count", position = "dodge")+theme_paper +facet_grid(. ~ celltype)+scale_fill_manual(values = rev(c("#FF6B6B", "#6BABFF")))+ylab("Number of DE mRNAs")+xlab("")



ggplot(de_stats_direction_plot[de_stats_direction_plot$type == "lnc",], aes(x = stage_short, fill = direction))+geom_bar(stat="count", position = "dodge")+theme_paper +facet_grid(. ~ celltype)+scale_fill_manual(values = rev(c("#FF6B6B", "#6BABFF")))+ylab("Number of DE lncRNAs")+xlab("")
ggplot(de_stats_direction_plot[de_stats_direction_plot$type == "lnc",], aes(x = stage_short, fill = direction))+geom_bar(stat="count", position = "dodge")+theme_paper +facet_grid(. ~ celltype+subtype)+scale_fill_manual(values = rev(c("#FF6B6B", "#6BABFF")))+ylab("Number of DE lncRNAs")+xlab("")+theme(legend.position = "")


de_stats_direction_plot <- de_all %>% group_by(primerid, celltype,type) %>% dplyr::summarise(direction = toString((direction)))
de_stats_direction_plot <- de_all %>% group_by(primerid, celltype,type,subtype) %>% dplyr::summarise(direction = toString((direction)))

ggplot(de_stats_direction_plot[de_stats_direction_plot$type == "lnc",], aes(x = celltype, col = direction, fill = direction))+geom_bar(stat="count", position= position_dodge())+theme_paper


ggplot(de_stats_direction_plot[de_stats_direction_plot$type == "pc",], aes(x = celltype, col = direction, fill = direction))+geom_bar(stat="count", position= position_dodge())+theme_paper

# Define which is the threshold to define when a gene is considered in CIS
cis_threshold <- 1000000
colocation <- readRDS("/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/07_colocation/colocation_lncvsDE_df.rds")
isg <- read.table("../isg.csv", stringsAsFactors = F)
distances <- data.frame(colocation, stringsAsFactors = F)
colnames(distances)  <-  c("lnc", "gene", "d")
distances$d <- as.integer(distances$d)

# Add gene names
genes <- ref[ref$type =="gene",]
correspondence <- data.frame(genes$gene_id, genename=genes$gene_name, stringsAsFactors = F)
rownames(correspondence) <- correspondence$genes.gene_id
distances$gene_name <- unlist(correspondence[distances$gene, ]$genename)

# Add lnc ortholog
orth_lnc <- data.frame(lnc = unique(distances$lnc), stringsAsFactors = F)
orth_lnc$orth <- unlist(lapply(orth_lnc$lnc, get_orthologname_))
rownames(orth_lnc) <- orth_lnc$lnc

# Check if ISG
distances$ISG <- ifelse(distances$gene_name %in% isg$V1, TRUE, FALSE)
distances$DEPC <- ifelse(distances$gene_name %in% gsub("-unknown", "", de_all[de_all$type == "pc",]$primerid), TRUE, FALSE)
distances[distances$gene  %in% gsub("-unknown", "", de_all[de_all$type == "pc",]$primerid),]$DEPC <- TRUE

comparisons_lnc_summary <- comparisons_lnc %>% dplyr::group_by(primerid) %>% dplyr::mutate(comparison_p = paste0(unique(celltype), collapse = ","))
comparisons_lnc_summary <- as.data.frame(distinct(comparisons_lnc_summary[,c("primerid", "comparison_p")]))
rownames(comparisons_lnc_summary) <- gsub("-unknown", "", comparisons_lnc_summary$primerid)
distances$our_module <- comparisons_lnc_summary[distances$lnc,]$comparison_p


# Explore
isg_dist <- distances[ distances$ISG == T & !is.na(distances$d) ,]
de_dist <- distances[distances$DEPC == T & !is.na(distances$d) ,]
close_isg <- isg_dist[isg_dist$d < cis_threshold ,]
close_isg_summary <- close_isg %>% dplyr::group_by(lnc) %>% dplyr::mutate(close = paste0(unique(gene_name), collapse = ","))
close_isg_summary <- distinct(close_isg_summary[,c("lnc","close")])
close_isg_summary <- as.data.frame(close_isg_summary)
rownames(close_isg_summary) <- close_isg_summary$lnc


close_depc <- de_dist[de_dist$d < cis_threshold ,]
close_depc_summary <- close_depc %>% dplyr::group_by(lnc) %>% dplyr::mutate(close = paste0(unique(gene_name), collapse = ","))
close_depc_summary <- distinct(close_depc_summary[,c("lnc","close")])
close_depc_summary <- as.data.frame(close_depc_summary)
rownames(close_depc_summary) <- close_depc_summary$lnc


# Add isg names
de_dist <- distances[ distances$DEPC == T & distances$d < cis_threshold ,]
#---------- HEATMAPS MAIN FIGURE
pal_celltypes <-wes_palette("GrandBudapest1", 4)
pal_celltypes[3] <- "#AE4E4E"
pal_celltypes <- pal_celltypes[c(2,1,3,4)]

png(filename = "/home/luisas/Desktop/cluster/data/plots/neut/T_lnc_FC.png", width = 1300, height = 1100)
FC_heatmap_subset_celltype(m, "T", unique(unlist(de_lists_lnc_T)), "T",  pal_celltypes[2],orthologs, lncpedia, reported_immlnc, legend = T)
dev.off()

png(filename = "/home/luisas/Desktop/cluster/data/plots/neut/B_lnc_FC.png", width = 1300, height = 1100)
FC_heatmap_subset_celltype(m, "B", unique(unlist(de_lists_lnc_B)), "B",  pal_celltypes[1],orthologs, lncpedia, reported_immlnc, legend = T )
dev.off()

png(filename = "/home/luisas/Desktop/cluster/data/plots/neut/Mono_lnc_FC.png", width = 1300, height = 1100)
FC_heatmap_subset_celltype(m, "Monocyte", unique(unlist(de_lists_lnc_my)), "Monocyte",  pal_celltypes[3],orthologs, lncpedia, reported_immlnc, legend = T )
dev.off()

