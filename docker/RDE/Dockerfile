FROM r-base:3.6.2
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>

RUN R -e "install.packages(c('readr'), repos='https://cloud.r-project.org/')"
RUN apt-get update && \
  apt-get install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev
RUN  R -e "install.packages(c('devtools'), repos='https://cloud.r-project.org/')"


RUN  R -e "install.packages(c('MAST', 'SingleCellExperiment', 'Seurat'), repos='https://cloud.r-project.org/')"

RUN  R -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org/')"
RUN  R -e "BiocManager::install(c('MAST','SingleCellExperiment', 'Seurat'), update=FALSE, ask=FALSE)"
RUN  R -e "BiocManager::install(c('data.table', 'stringr'), update=FALSE, ask=FALSE)"
RUN  R -e "BiocManager::install(c('DESeq2', 'zinbwave'), update=FALSE, ask=FALSE)"
RUN  R -e "BiocManager::install(c('scater', 'purrr', 'rtracklayer', 'stringr'), update=FALSE, ask=FALSE)"
RUN  R -e "BiocManager::install(c('reshape2'), update=FALSE, ask=FALSE)"
RUN  R -e "BiocManager::install(c('Matrix'), update=FALSE, ask=FALSE)"