
# converts a bed file after geneToPred to a bed12 file.
library(readr)
# Read command line
args = commandArgs(trailingOnly=TRUE)



get_difference <- function(row){
  difference_list <- as.double(strsplit(toString(row$V12), ",")[[1]])-as.double(strsplit(toString(row$V11), ",")[[1]])
  difference_list_string <- paste0(paste(difference_list, collapse = ","), ",")
  return(difference_list_string)
}

# Apply modifications
file <- read.table(args[1])
file <- as.data.frame(file)
file$V13 <- unlist(lapply(1:nrow(file), function(r_num) { get_difference(file[r_num, ]) }))
file$V14 <- file$V11
file$V11 <-NULL
file$V12 <-NULL

#save
outpath <- args[2]
write_tsv(file, outpath, col_names = FALSE)