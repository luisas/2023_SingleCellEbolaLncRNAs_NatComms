
FROM ubuntu:16.04
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>

RUN apt-get update && apt-get install -y wget &&\
                      apt-get install -y software-properties-common &&\
                      apt-get install -y unzip &&\
                      apt-get install -y openjdk-8-jdk

RUN wget https://github.com/broadinstitute/Drop-seq/releases/download/v2.3.0/Drop-seq_tools-2.3.0.zip

RUN unzip Drop-seq_tools-2.3.0.zip && rm Drop-seq_tools-2.3.0.zip

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

ENV PATH "/Drop-seq_tools-2.3.0:${PATH}"
