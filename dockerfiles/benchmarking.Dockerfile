FROM ubuntu:16.04
LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y gawk wget ruby-full ruby-bundler && gem install erubis

# Setup evaluation framework
#overwrite at build time with e.g. --build-arg EVALU_VERSION=0.2
ARG EVALU_VERSION=0.9
WORKDIR /benchmarking
RUN wget https://github.com/rsuchecki/biokanga_benchmark/archive/v${EVALU_VERSION}.tar.gz \
  && tar xzvf v${EVALU_VERSION}.tar.gz \
  && mv biokanga_benchmark-${EVALU_VERSION}/* . \
  && rmdir biokanga_benchmark-${EVALU_VERSION}

ENV PATH "$PATH:/benchmarking"

# ARG SAMTOOLS_VERSION=1.9
# WORKDIR /tmp
# ADD https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}}/samtools-${SAMTOOLS_VERSION}.tar.bz2 .

# RUN apt-get install -y ncurses ncurses-dev musl-dev g++ make zlib-dev \
#     && tar xjvf samtools-${TOOL_VERSION}.tar.bz2 \
#     && cd /tmp/samtools-${TOOL_VERSION} && make \
#     && mv /tmp/samtools-${TOOL_VERSION}/samtools /usr/bin \
#     && rm -rf /tmp/*
RUN apt install --yes libcurl3-gnutls \
  && wget http://mirrors.kernel.org/ubuntu/pool/universe/s/samtools/samtools_1.4.1-1build1_amd64.deb \
  && wget http://mirrors.kernel.org/ubuntu/pool/universe/h/htslib/libhts2_1.5-1_amd64.deb \
  && dpkg -i samtools_*.deb libhts2_*.deb \
  && rm *.deb \
  && apt clean