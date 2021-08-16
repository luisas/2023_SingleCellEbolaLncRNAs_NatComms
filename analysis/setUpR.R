
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c( "rtracklayer", "stringr","ggplot2", "grid", "gridExtra", "RColorBrewer", "readr", "matrixStats",
                            "GenomicRanges", "dplyr", "zeallot", "tximport", "ggpubr", "plyr",
                            "formattable" ))
