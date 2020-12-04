FROM ubuntu:16.04
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>



RUN apt-get update && apt-get install -y zlib1g-dev build-essential wget

RUN apt-get update && apt-get -y upgrade

RUN apt-get install -y python3-pip

RUN pip3 install --upgrade pip

RUN pip3 install pandas llvmlite==0.31.0 scrublet

