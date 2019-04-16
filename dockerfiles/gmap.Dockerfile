FROM ubuntu:16.04

LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"

#Adapted from genomicpariscentre/gmap with thanks to author: Aurélien Birer

ARG VERSION="2019-03-15"
WORKDIR /tmp

RUN apt update \
  && apt install --yes \
    make \
    g++ \
    wget \
    perl \
  && wget http://research-pub.gene.com/gmap/src/gmap-gsnap-${VERSION}.tar.gz \
  && tar xzvf gmap-gsnap-${VERSION}.tar.gz \
  && cd gmap-${VERSION} \
  && ./configure --prefix=/usr/local \
  && make \
  && make check \
  && make install \
  && rm -rf /tmp/* \
  && apt remove --purge --yes \
    make \
    g++ \
    wget \
  && apt autoremove --purge --yes