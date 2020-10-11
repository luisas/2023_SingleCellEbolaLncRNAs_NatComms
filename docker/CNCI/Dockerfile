FROM ubuntu:16.04
MAINTAINER Luisa Santus <luisa.santus95@gmail.com>

RUN apt-get update \
  && apt-get install -y python-pip python-dev wget r-base \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python2.7.17 python \
  && pip install --upgrade pip

RUN pip install numpy sklearn xgboost

RUN wget http://cnit.noncode.org/CNIT/data/CNIT.tar.gz \
    && tar -zxvf CNIT.tar.gz

RUN pip install joblib

ENV CNIT_HOME="/CNCI2"

RUN cd CNCI2 \
    && grep -v "print \"" CNCI2.py > temp.py \
    && mv temp.py CNCI2.py

RUN sed -i "s/self.n_jobs/self.nthread/g" /usr/local/lib/python2.7/dist-packages/xgboost/sklearn.py
