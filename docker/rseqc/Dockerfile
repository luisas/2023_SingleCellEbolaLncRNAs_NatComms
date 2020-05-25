FROM ubuntu:16.04
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>



RUN apt-get update \
  && apt-get install -y wget \
	&& apt-get install -y gcc \
	&& apt-get install -y python-dev \
	&& apt-get install -y zlib1g-dev libssl-dev \
	&& apt-get install -y samtools r-base\
	&& wget https://bootstrap.pypa.io/get-pip.py \
	&& python get-pip.py \
	&& pip install RSeQC \
	&& rm get-pip.py \
	&& rm -rf /var/lib/apt/lists/*
