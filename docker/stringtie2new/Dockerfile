FROM ubuntu:16.04
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>


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
       git \
       python3 \
       libfindbin-libs-perl \
       python3-pip \
       unzip \
       vim-common \
       wget \
       zlib1g-dev \
   && rm -rf /var/lib/apt/lists/*



# Stringtie
RUN git clone https://github.com/mpertea/stringtie2 && \
    cd stringtie2 && \
    make release


# Gff compare
RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.10.6.Linux_x86_64.tar.gz &&\
   tar -xzvf gffcompare-0.10.6.Linux_x86_64.tar.gz &&\
   rm gffcompare-0.10.6.Linux_x86_64.tar.gz

ENV PATH="/gffcompare-0.10.6.Linux_x86_64/:${PATH}"


RUN mkdir -p /scratch
