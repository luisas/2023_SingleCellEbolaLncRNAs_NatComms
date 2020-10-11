FROM ubuntu:16.04
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>



RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
    build-essential \
    ca-certificates \
    gcc \
    git \
    libpq-dev \
    make \
    python-pip \
    python2.7 \
    python2.7-dev \
    ssh \
    && apt-get autoremove \
    && apt-get clean

RUN pip2 install "numpy==1.16.6"
RUN pip2 install "biopython==1.76"

RUN wget http://cpc2.gao-lab.org/data/CPC2-beta.tar.gz \
    && gzip -dc CPC2-beta.tar.gz | tar xf - \
    && cd CPC2-beta/libs/libsvm \
    && gzip -dc libsvm-3.18.tar.gz | tar xf - \
    && cd libsvm-3.18 \
    && make clean && make

ENV CPC_HOME="/CPC2-beta/bin"

RUN cd $CPC_HOME \
    && sed -i "s/if not os.path.isfile(options.fasta):/if (not os.path.isfile(options.fasta) and os.path.islink(options.fasta)):/" CPC2.py
