FROM r-base:3.6.2
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>

RUN R -e "install.packages(c('readr'), repos='https://cloud.r-project.org/')"
RUN apt-get update && \
  apt-get install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev
RUN  R -e "install.packages(c('devtools'), repos='https://cloud.r-project.org/')"


RUN R -e "library(devtools)"

RUN  R -e "install.packages(c('caTools'), repos='https://cloud.r-project.org/')"

RUN R -e "library(devtools) ; install_github('rajewsky-lab/dropbead')"
