FROM ubuntu:16.04

MAINTAINER Luisa Santus <luisa.santus95@gmail.com>
# -----------------------------
# -----------------------------
# Install JAVA and FASTQC
# commands taken from
RUN apt-get update && apt-get install -y software-properties-common

RUN apt-get update && \
	apt-get install -y openjdk-8-jre && \
	rm -rf /var/lib/apt/lists/*

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN apt-get -qq update && apt-get -y upgrade && \
	apt install -y wget libfindbin-libs-perl software-properties-common unzip


# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -O /tmp/fastqc.zip && \
    unzip /tmp/fastqc.zip -d /opt/ && \
    rm /tmp/fastqc.zip && \
    chmod 777 /opt/FastQC/fastqc

ENV PATH="/opt/FastQC/:${PATH}"

# -----------------------------
# -----------------------------
# Install PICARD
# Commands taken from broadinstitute/picard

RUN mkdir /apps && mkdir /apps/PICARD && mkdir /apps/PICARD/2.20.0
RUN wget https://github.com/broadinstitute/picard/releases/download/2.20.0/picard.jar &&\
    mv picard.jar /apps/PICARD/2.20.0/


# -----------------------------
# -----------------------------
# part of the commands taken from kathrinklee/rna-seq-pipeline-hisat2
# Install hisat2
RUN apt-get update  && apt-get install -y \
    build-essential \
    gcc-multilib \
    zlib1g-dev \
    curl \
    wget \
    cmake \
    python \
    python-pip \
    python-dev \
    python2.7-dev \
    python-numpy \
    python-matplotlib \
    hdf5-tools \
    libhdf5-dev \
    hdf5-helpers \
    libhdf5-serial-dev \
    libssh2-1-dev \
    libcurl4-openssl-dev \
    icu-devtools \
    libssl-dev \
    libxml2-dev \
    r-bioc-biobase \
    git \
    apt-utils \
    pigz
# Install Python
RUN add-apt-repository ppa:jonathonf/python-3.6 -y &&\
    apt-get update &&\
    apt-get install python3.6 -y &&\
    cd /usr/bin && unlink python && ln -s /usr/bin/python3.6 python && cd


RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
    unzip hisat2-2.1.0-Linux_x86_64.zip
RUN cp -p hisat2-2.1.0/hisat2 hisat2-2.1.0/hisat2-* /usr/bin

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -jxf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    make && \
    make install && \
    cp samtools /usr/bin/

# Install stringtie
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.6.Linux_x86_64.tar.gz && \
    tar zxf stringtie-1.3.6.Linux_x86_64.tar.gz && \
    cp ./stringtie-1.3.6.Linux_x86_64/stringtie /usr/bin/

# Remove useless things
RUN rm -rf hisat2-2.1.0
RUN rm -rf samtools-1.9
RUN rm -rf stringtie-1.3.6.Linux_x86_64

RUN apt-get update &&\
	  apt-get -y install python3-pip

RUN pip install umi_tools



#wget http://sourceforge.net/projects/rseqc/files/RSeQC-2.6.4.tar.gz &&\
# tar zxf RSeQC-2.6.4.tar.gz
