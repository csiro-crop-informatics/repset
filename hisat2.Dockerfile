FROM ubuntu:16.04
LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y wget unzip python2.7

# Install hisat2
ARG HISAT_VERSION=2.0.5
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-${HISAT_VERSION}-Linux_x86_64.zip \
  && unzip hisat2-${HISAT_VERSION}-Linux_x86_64.zip \
  && rm hisat2-${HISAT_VERSION}-Linux_x86_64.zip
RUN cp -p hisat2-${HISAT_VERSION}/hisat2 hisat2-${HISAT_VERSION}/hisat2-* /usr/bin
