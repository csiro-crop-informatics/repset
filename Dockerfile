FROM ubuntu:16.04
LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y wget build-essential automake make unzip pigz ruby-full

# Install hisat2
ARG HISAT_VERSION=2.0.5
WORKDIR /docker
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-${HISAT_VERSION}-Linux_x86_64.zip \
    && unzip hisat2-${HISAT_VERSION}-Linux_x86_64.zip \
    && rm hisat2-${HISAT_VERSION}-Linux_x86_64.zip
#RUN cp -p hisat2-${HISAT_VERSION}/hisat2 hisat2-${HISAT_VERSION}/hisat2-* /usr/bin
#ENV PATH "$PATH:/docker/hisat2-2.0.5/"

#Install biokanga
#overwrite at build time with e.g. --build-arg KANGA_VERSION=4.3.11
ARG KANGA_VERSION=4.3.10
WORKDIR /docker
RUN wget https://github.com/csiro-crop-informatics/biokanga/archive/v${KANGA_VERSION}.tar.gz \
  && tar xzf v${KANGA_VERSION}.tar.gz \
  && cd biokanga-${KANGA_VERSION} \
  && autoreconf -f -i \
  && ./configure \
  && make
#ENV PATH "$PATH:/docker/biokanga-${KANGA_VERSION}/biokanga"

# Setup evaluation framework
ARG EVALU_VERSION=0.3
WORKDIR /project/itmatlab/aligner_benchmark
RUN wget https://github.com/rsuchecki/biokanga_benchmark/archive/v${EVALU_VERSION}.tar.gz \
  && tar xzvf v${EVALU_VERSION}.tar.gz \
  && mv biokanga_benchmark-${EVALU_VERSION}/* . \
  && rmdir biokanga_benchmark-${EVALU_VERSION} \
  && gem install erubis \
  && mkdir -p /project/itmatlab/aligner_benchmark/dataset/human/genome /project/itmatlab/aligner_benchmark/tool_results /project/itmatlab/aligner_benchmark/statistics/
