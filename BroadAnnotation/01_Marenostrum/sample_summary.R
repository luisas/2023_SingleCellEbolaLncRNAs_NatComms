
library(R.utils)
#pathstringtie <- "/home/luisas/Desktop/cluster/data2/00_RawData/BroadTranscriptomesComplete/stringtie2-gtfs/"
pathstringtie <- "/gpfs/projects/bsc83/Data/Ebola/00_RawData/BroadTranscriptomesComplete/stringtie2-gtfs"
samples_available <- read.table(file.path(pathstringtie,"available_samples.txt"))
complete_sample_table <- read.csv(file.path(pathstringtie,"/Master_Sample-Keys_TotalRNA_EBOV - Master_Sample-Keys_TotalRNA_EBOV.csv"))

# Add ID processing column
samples_available$ID.processing  <- ifelse(substr(samples_available$V1,1,1) == "A", 
                                                 substr(samples_available$V1,1,5), 
                                                 NA)


# Samples available - ID processing availble
samples_available_A <- samples_available [!is.na(samples_available$ID.processing ),]
samples_available_O <- samples_available [is.na(samples_available$ID.processing ),]

samples_available_A_info <- complete_sample_table[complete_sample_table$ID.Processing %in% samples_available_A$ID.processing,]
samples_available_A_info$X <- samples_available_A[samples_available_A$ID.processing %in% complete_sample_table$ID.Processing,]$V1


#write.csv(samples_available_O, file.path(pathstringtie,"/missing.csv"))

missing_labeled <- read.table(file.path(pathstringtie,"/missing_labeled.csv"), sep = ",")
# Check wether we have all the samples 
nrow(samples_available[samples_available$V1 != "folderdump", ]) == nrow(samples_available_A_info)+nrow(missing_labeled)

# prepare for merging
samples_available_A_info$X
names(missing_labeled) <- c("X", "Biosample", "Tissue", "ID.Cohort")
all <- rbind(samples_available_A_info[,c("X", "Biosample", "Tissue", "ID.Cohort")], missing_labeled)
total_number_samples <- nrow(all)
all[all[!is.na(all$Tissue),]$Tissue == "Whole blood",]$Tissue <- "Wholeblood"
all$infected <- NA
all[all$ID.Cohort %in% c("D000", "D-14", "D-30", "D-30?"),]$infected <- "Not infected"
all[!is.na(all$ID.Cohort) & is.na(all$infected),]$infected <- "infected"
all$X

#ggplot(all, aes(y = Tissue, fill = infected))+geom_bar()+theme_classic()+theme(legend.position = "top")

# Create subset for samples 
dir.create(file.path(pathstringtie,"/benchmark_subsets"), showWarnings = F)

# Test tissue effect 
# Select 3 tissue with enough samples

wholeblood <- all[all$Tissue== "Wholeblood",]
wholeblood <- wholeblood[!is.na(wholeblood$X),]

spleen <- all[all$Tissue == "Spleen",]
spleen <- spleen[!is.na(spleen$X),]

lymphnode <- all[all$Tissue == "Lymph node",]
lymphnode <- lymphnode[!is.na(lymphnode$X),]


n <- min(nrow(wholeblood), nrow(spleen), nrow(lymphnode))

make_link <- function(name,dir){
  prefix <- file.path(pathstringtie, "stringtie2-gtfs", name)
  file <- paste(prefix, "_stringtie2.gtf", sep="")
  print(file)
  link <- file.path(dir, paste(name, "_stringtie2.gtf", sep = "") )
  createLink(link=link, target = file )
}

wholeblood_names <- as.character(wholeblood[sample(nrow(wholeblood), n),]$X)
spleen_names <- as.character(spleen[sample(nrow(spleen), n),]$X)
lymphnode_names <- as.character(lymphnode[sample(nrow(lymphnode), n),]$X)

length(unique(wholeblood_names))


# ----- CREATE LINKS ------
dir.create(file.path(pathstringtie, "/benchmark_subsets/tissue_effect"), showWarnings = F)

dir_1 <- file.path(pathstringtie,"/benchmark_subsets/tissue_effect/00_spleen")
dir.create(dir_1, showWarnings = F)
# Create links for dir 1 
lapply(wholeblood_names, function(name) make_link(name, dir_1) ) 

dir_2 <- file.path(pathstringtie,"/benchmark_subsets/tissue_effect/01_spleen_wholeblood")
dir.create(dir_2, showWarnings = F)
lapply(c(wholeblood_names,spleen_names), function(name) make_link(name, dir_2) ) 

dir_3 <- file.path(pathstringtie,"/benchmark_subsets/tissue_effect/02_spleen_wholeblood_lymphnode") 
dir.create(dir_3, showWarnings = F)
lapply(c(wholeblood_names,spleen_names,lymphnode_names), function(name) make_link(name, dir_3) ) 



### Second benchmark 
# ----- CREATE LINKS ------
dir.create(file.path(pathstringtie, "/benchmark_subsets/inf-tissue-sample"), showWarnings = F)


wholeblood <- all[all$Tissue== "Wholeblood",]
wholeblood <- wholeblood[!is.na(wholeblood$X),]
wholeblood_ni <- wholeblood[wholeblood$infected == "Not infected",]
wholeblood_i <- wholeblood[wholeblood$infected == "infected",]

spleen_ni <- spleen[spleen$infected == "Not infected",]
spleen_i <- spleen[spleen$infected == "infected",]

lymphnode_ni <- lymphnode[lymphnode$infected == "Not infected",]
lymphnode_i <- lymphnode[lymphnode$infected == "infected",]

# -------------------------------
# SETTING 1 
# -------------------------------

# WB add 10 samples (uninfected) at the time
wb_ni_30 <- as.character(wholeblood_ni[sample(nrow(wholeblood_ni), 30),]$X)
wb_ni_20 <- wb_ni_30[11:30]
wb_ni_10 <- wb_ni_30[1:10]

dir <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/01_addsample")
dir.create(dir, showWarnings = F)

dir_1 <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/01_addsample/00_wbni_10")
dir.create(dir_1, showWarnings = F)
lapply(wb_ni_10, function(name) make_link(name, dir_1) ) 

dir_2 <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/01_addsample/01_wbni_20")
dir.create(dir_2, showWarnings = F)
lapply(wb_ni_20, function(name) make_link(name, dir_2) ) 

dir_3 <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/01_addsample/02_wbni_30")
dir.create(dir_3, showWarnings = F)
lapply(wb_ni_30, function(name) make_link(name, dir_3) ) 



# -------------------------------
# SETTING 1 
# -------------------------------

# Baseline 10 not infected 
# add 5 infected and 5 not infected 
# Again: add 5 infected and 5 not infected 
wb_ni_10
wb_10_inf <- as.character(wholeblood_i[sample(nrow(wholeblood_i), 10),]$X)
wb_5a_inf <- wb_10_inf[1:5]
wb_5b_inf <- wb_10_inf[6:10]
wb_5a_ni <- wb_ni_30[11:15]
wb_5b_ni <- wb_ni_30[16:20]
wb_infnotinf_10_10 <- c(wb_ni_10,wb_5a_inf, wb_5a_ni )
wb_infnotinf_10_20 <- c(wb_infnotinf_10_10,wb_5b_inf, wb_5b_ni )

dir <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/02_add_sampleMIXinf")
dir.create(dir, showWarnings = F)

dir_1 <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/02_add_sampleMIXinf/00_wbni_10")
dir.create(dir_1, showWarnings = F)
lapply(wb_ni_10, function(name) make_link(name, dir_1) ) 

dir_2 <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/02_add_sampleMIXinf/01_wbMIXinf_20")
dir.create(dir_2, showWarnings = F)
lapply(wb_infnotinf_10_10, function(name) make_link(name, dir_2) ) 

dir_3 <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/02_add_sampleMIXinf/02_wbMINinf_30")
dir.create(dir_3, showWarnings = F)
lapply(wb_infnotinf_10_20, function(name) make_link(name, dir_3) ) 



# Baseline 10 sample not infected whole blood 
# add 10 samples from a different tissue at the time.
wb_ni_10
wb_ni_spleen_10 <- as.character(spleen[sample(nrow(spleen), 10),]$X)
wb_ni_lymphnode_10 <- as.character(lymphnode[sample(nrow(lymphnode), 10),]$X)


dir <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/03_add_tissueNI")
dir.create(dir, showWarnings = F)

dir_1 <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/03_add_tissueNI/00_wbni_10")
dir.create(dir_1, showWarnings = F)
lapply(wb_ni_10, function(name) make_link(name, dir_1) ) 

dir_2 <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/03_add_tissueNI/01_wbspleen")
dir.create(dir_2, showWarnings = F)
lapply(c(wb_ni_10,wb_ni_spleen_10), function(name) make_link(name, dir_2) ) 

dir_3 <- file.path(pathstringtie,"/benchmark_subsets/inf-tissue-sample/03_add_tissueNI/02_wbspleenln")
dir.create(dir_3, showWarnings = F)
lapply(c(wb_ni_10,wb_ni_spleen_10,wb_ni_lymphnode_10), function(name) make_link(name, dir_3) ) 



### All samples (10 by 10 )

dir <- file.path(pathstringtie,"/benchmark_subsets/all/")
dir.create(dir, showWarnings = F)

all <- all[!is.na(all$X),]
all_tmp <- all 

# Initialize vars
dir_n <- 1
names <- c()
missing_samples <- length(unique(all_tmp$X))

while(missing_samples > 10){
  
  n_samples <- ifelse(length(unique(all_tmp$X))>=10, 10, length(unique(all_tmp$X)))
  a <- as.character(all_tmp[sample(nrow(all_tmp), n_samples),]$X)
  print("a")
  all_tmp <- all_tmp[!(all_tmp$X %in% a), ]
  
  names <- c(names, a)
  dir <- file.path(pathstringtie,"/benchmark_subsets/all/",dir_n)
  dir.create(dir, showWarnings = F)
  lapply(names, function(name) make_link(name, dir) ) 
  dir_n <- dir_n+1
  missing_samples <- length(unique(all_tmp$X))
  
}

# ---------------------------------------
#     TISSUE BENCHMARK (2nd design) 
# ---------------------------------------

# Add tissues
create_dir_subset <- function(names, n, name, subdir = "/benchmark_subsets/Tissue"){
  dir <- file.path(pathstringtie,subdir, paste(name, n, sep = "_"))
  dir.create(dir, showWarnings = F)
  names <- lapply(names, function(x) str_replace(x,"X",""))
  lapply(names, function(name) make_link(name, dir) )
}

# Add samples to each tissue 6 by 6 
wb_split_6 <- data.frame(split(wholeblood$X, ceiling(seq_along(wholeblood$X)/6)))
sp_split_6 <- data.frame(split(spleen$X, ceiling(seq_along(spleen$X)/6)))
ln_split_6 <- data.frame(split(lymphnode$X, ceiling(seq_along(lymphnode$X)/6)))


# Add wholeblood 6 by 6 
lapply(unlist(lapply(names(wb_split_6), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(wb_split_6[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*6), "00_wb" ) )
# Add spleen 6 by 6 
lapply(unlist(lapply(names(sp_split_6), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(sp_split_6[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*6), "01_sp" ) )
# Add ln 6 by 6 
lapply(unlist(lapply(names(ln_split_6), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(ln_split_6[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*6), "02_ln" ) )


# Add 3+3 wb sp
wb_split_3 <- data.frame(split(wholeblood$X, ceiling(seq_along(wholeblood$X)/3)))
sp_split_3 <- data.frame(split(spleen$X, ceiling(seq_along(spleen$X)/3)))
ln_split_3 <- data.frame(split(lymphnode$X, ceiling(seq_along(lymphnode$X)/3)))

# Add 2+2+2
wb_split_2 <- data.frame(split(wholeblood$X, ceiling(seq_along(wholeblood$X)/2)))
sp_split_2 <- data.frame(split(spleen$X, ceiling(seq_along(spleen$X)/2)))
ln_split_2 <- data.frame(split(lymphnode$X, ceiling(seq_along(lymphnode$X)/2)))

prep_df_3_3 <- function(set1, set2, set3 = NULL){
  if(is.null(set3)){
    set1 <- data.frame(set1)
    set2 <- data.frame(set2)
    intersection <- (intersect(names(set1), names(set2)))
    df <- rbind(set1[, intersection],set2[, intersection]) 
    list <- as.list(df)
    names(list) <- lapply(colnames(df), function(x) str_replace(x,"X",""))
  }else{
    intersection <- (intersect(intersect(names(set1), names(set2)), names(set3)))
    df <- rbind(set1[, intersection],set2[, intersection], set3[, intersection]) 
    list <- as.list(df)
    names(list) <- lapply(colnames(df), function(x) str_replace(x,"X",""))
  }
  return(data.frame(list))
}

wb_sp_3_3 <- prep_df_3_3(wb_split_3,sp_split_3 )
lapply(unlist(lapply(names(wb_sp_3_3), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(wb_sp_3_3[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*6), "03_wb_sp" ) )

wb_ln_3_3 <- prep_df_3_3(wb_split_3,ln_split_3 )
lapply(unlist(lapply(names(wb_ln_3_3), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(wb_ln_3_3[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*6), "04_wb_ln" ) )

ln_sp_3_3 <- prep_df_3_3(ln_split_3,sp_split_3 )
lapply(unlist(lapply(names(ln_sp_3_3), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(ln_sp_3_3[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*6), "05_ln_sp" ) )


wb_sp_ln_2_2_2 <- prep_df_3_3(wb_split_2,sp_split_2, ln_split_2)
lapply(unlist(lapply(names(wb_sp_ln_2_2_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(wb_sp_ln_2_2_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*6), "06_wb_sp_ln" ) )

# ---------------------------------------
#     Infection BENCHMARK (2nd design) 
# ---------------------------------------


wb_i_split_2 <- data.frame(split(wholeblood_i$X, ceiling(seq_along(wholeblood_i$X)/2)))
lapply(unlist(lapply(names(wb_i_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(wb_i_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*2),  "00_wb" ,"/benchmark_subsets/Infection" ) )



wb_ni_split_2 <- data.frame(split(wholeblood_ni$X, ceiling(seq_along(wholeblood_ni$X)/2)))
lapply(unlist(lapply(names(wb_ni_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(wb_ni_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*2),  "01_wb_notinf" ,"/benchmark_subsets/Infection" ) )



wb_i_split_1 <- data.frame(split(wholeblood_i$X, ceiling(seq_along(wholeblood_i$X)/1)))
wb_ni_split_1 <-  data.frame(split(wholeblood_ni$X, ceiling(seq_along(wholeblood_ni$X)/1)))

wb_i_ni_split_2 <- prep_df_3_3(wb_i_split_1,wb_ni_split_1 )
lapply(unlist(lapply(names(wb_i_ni_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(wb_i_ni_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*2),  "02_wb_inf_notinf" ,"/benchmark_subsets/Infection" ) )


# Spleen 
sp_i_split_2 <- data.frame(split(spleen_i$X, ceiling(seq_along(spleen_i$X)/2)))
lapply(unlist(lapply(names(sp_i_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(sp_i_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*2),  "03_sp" ,"/benchmark_subsets/Infection" ) )

sp_ni_split_2 <- data.frame(split(spleen_ni$X, ceiling(seq_along(spleen_ni$X)/2)))
lapply(unlist(lapply(names(sp_ni_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(sp_ni_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*2),  "04_sp_notinf" ,"/benchmark_subsets/Infection" ) )

sp_i_split_1 <- data.frame(split(spleen_i$X, ceiling(seq_along(spleen_i$X)/1)))
sp_ni_split_1 <-  data.frame(split(spleen_ni$X, ceiling(seq_along(spleen_ni$X)/1)))

sp_i_ni_split_2 <- prep_df_3_3(sp_i_split_1,sp_ni_split_1 )
lapply(unlist(lapply(names(sp_i_ni_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(sp_i_ni_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*2),  "05_sp_inf_notinf" ,"/benchmark_subsets/Infection" ) )

# Lymphnode 
ln_i_split_2 <- data.frame(split(lymphnode_i$X, ceiling(seq_along(lymphnode_i$X)/2)))
lapply(unlist(lapply(names(ln_i_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(ln_i_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*2),  "06_ln" ,"/benchmark_subsets/Infection" ) )

ln_ni_split_2 <- data.frame(split(lymphnode_ni$X, ceiling(seq_along(lymphnode_ni$X)/2)))
lapply(unlist(lapply(names(ln_ni_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(ln_ni_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*2),  "07_ln_notinf" ,"/benchmark_subsets/Infection" ) )

ln_i_split_1 <- data.frame(split(lymphnode_i$X, ceiling(seq_along(lymphnode_i$X)/1)))
ln_ni_split_1 <-  data.frame(split(lymphnode_ni$X, ceiling(seq_along(lymphnode_ni$X)/1)))

ln_i_ni_split_2 <- prep_df_3_3(ln_i_split_1,ln_ni_split_1 )
lapply(unlist(lapply(names(ln_i_ni_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(ln_i_ni_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*2),  "08_ln_inf_notinf" ,"/benchmark_subsets/Infection" ) )

# ---------------------------------------
#     Infection+Tissue BENCHMARK (2nd design) 
# ---------------------------------------
wb_ln_inf_split_2 <- (prep_df_3_3(wb_i_split_2,ln_i_split_2 ))
wb_ln_notinf_split_2 <- (prep_df_3_3(wb_ni_split_2,ln_ni_split_2 ))
lapply(unlist(lapply(names(wb_ln_inf_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(wb_ln_inf_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*4),  "00_wb_ln_inf" ,"/benchmark_subsets/InfectionANDTissue" ) )
lapply(unlist(lapply(names(wb_ln_notinf_split_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(wb_ln_notinf_split_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*4),  "01_wb_ln_notinf" ,"/benchmark_subsets/InfectionANDTissue" ) )



wb_inf_ni_2 <- (prep_df_3_3(wb_i_split_1,wb_ni_split_1 ))
ln_inf_ni_2 <- (prep_df_3_3(ln_i_split_1,ln_ni_split_1 ))
wb_ln_inf_ni_2 <- prep_df_3_3(wb_inf_ni_2, ln_inf_ni_2)
lapply(unlist(lapply(names(wb_ln_inf_ni_2), function(x) str_replace(x,"X",""))), function(subset) create_dir_subset(unlist(wb_ln_inf_ni_2[,1:as.numeric(subset), drop = F]), as.character(as.numeric(subset)*4),  "02_wb_ln_mix" ,"/benchmark_subsets/InfectionANDTissue" ) )


