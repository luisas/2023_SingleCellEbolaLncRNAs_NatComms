# Docker image for RNASeq analysis
FROM ubuntu:16.04
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>

# ---- Contained ------
# FastQC v.0.11.5
# Picard v.2.2.0
# Hisat2 v.2.2.1
# Samtools v.1.9
# Stringtie v.1.3.6
# Umi_tools v1.1.0
# HTSeq v0.11.1

RUN apt-get update && apt-get install -y software-properties-common && 	apt-get install -y openjdk-8-jre &&\
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        openjdk-8-jdk \
        python3 \
				libfindbin-libs-perl \
        python3-pip \
        unzip \
        vim-common \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Java
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN apt-get -qq update && apt-get -y upgrade && \
	apt install -y wget libfindbin-libs-perl software-properties-common unzip

RUN apt-get install -y  --allow-unauthenticated bison byacc

RUN apt-get install -y make

#STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz &&\
    tar -xzvf 2.7.3a.tar.gz &&\
    rm ./2.7.3a.tar.gz


#Bioawk
RUN wget https://github.com/lh3/bioawk/archive/master.zip && unzip master.zip && cd bioawk-master && make && cp bioawk /usr/local/bin/


ENV PATH="/STAR-2.7.3a/bin/Linux_x86_64/:${PATH}"





RUN mkdir -p /scratch
