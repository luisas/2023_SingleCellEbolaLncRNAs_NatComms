FROM ubuntu:16.04
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>



RUN apt-get update \
  && apt-get install -y python3-pip python3-dev wget \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN wget https://sourceforge.net/projects/rna-cpat/files/v3.0.0/CPAT-3.0.0.tar.gz \
    && tar -zxvf CPAT-3.0.0.tar.gz

RUN apt-get install -y r-base 
RUN pip3 install CPAT
