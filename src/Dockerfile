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


# FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -O /tmp/fastqc.zip && \
    unzip /tmp/fastqc.zip -d /opt/ && \
    rm /tmp/fastqc.zip && \
    chmod 777 /opt/FastQC/fastqc

ENV PATH="/opt/FastQC/:${PATH}"


# Picard
RUN mkdir /apps && mkdir /apps/PICARD && mkdir /apps/PICARD/2.20.0
RUN wget https://github.com/broadinstitute/picard/releases/download/2.20.0/picard.jar &&\
    mv picard.jar /apps/PICARD/2.20.0/


#STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz &&\
    tar -xzvf 2.7.3a.tar.gz &&\
    cp ./2.7.3a.tar.gz /usr/bin/


## ---- OLD
# Hisat2
#RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
#    unzip hisat2-2.1.0-Linux_x86_64.zip
#RUN cp -p hisat2-2.1.0/hisat2 hisat2-2.1.0/hisat2-* /usr/bin
#RUN rm -rf hisat2-2.1.0-Linux_x86_64.zip

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -jxf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    make && \
    make install && \
    cp samtools /usr/bin/
RUN rm -rf samtools-1.9.tar.bz2


# Picard
RUN apt-get update && apt-get install -y picard-tools



# Stringtie
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.6.Linux_x86_64.tar.gz && \
    tar zxf stringtie-1.3.6.Linux_x86_64.tar.gz && \
    cp ./stringtie-1.3.6.Linux_x86_64/stringtie /usr/bin/
RUN rm -rf stringtie-1.3.6.Linux_x86_64.tar.gz

# Umi-tools
RUN apt-get update && \
    apt-get install --yes python3-pip && \
    pip3 install --upgrade pip

RUN pip3 install umi_tools==1.0.0 && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*


# HtSeq counts
RUN pip3 install HTSeq


# Gff compare
RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.10.6.Linux_x86_64.tar.gz &&\
    tar -xzvf gffcompare-0.10.6.Linux_x86_64.tar.gz &&\
    rm gffcompare-0.10.6.Linux_x86_64.tar.gz

ENV PATH="/gffcompare-0.10.6.Linux_x86_64/:${PATH}"




RUN mkdir -p /scratch
