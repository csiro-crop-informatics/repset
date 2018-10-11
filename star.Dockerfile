FROM ubuntu:16.04
LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y wget libgomp1

#overwrite at build time with e.g. --build-arg VERSION=1.3.1
ARG VERSION=2.6.1b

#Install dart
WORKDIR /star
RUN wget https://github.com/alexdobin/STAR/archive/${VERSION}.tar.gz \
  && tar xzvf ${VERSION}.tar.gz
ENV PATH "$PATH:/star/STAR-${VERSION}/bin/Linux_x86_64/"

