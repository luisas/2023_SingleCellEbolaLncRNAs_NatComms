FROM r-base:3.6.3
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>

RUN R -e "install.packages(c('readr'), repos='https://cloud.r-project.org/')"
RUN apt-get update && \
  apt-get install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev

RUN  R -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org/')"
RUN  R -e "BiocManager::install(c('tximport','rtracklayer', 'dplyr', 'matrixStats', 'forcats'), update=FALSE, ask=FALSE)"
