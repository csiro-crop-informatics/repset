FROM ubuntu:16.04
LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y gawk wget ruby-full ruby-bundler && gem install erubis

# Setup evaluation framework
#overwrite at build time with e.g. --build-arg EVALU_VERSION=0.2
ARG EVALU_VERSION=0.7
WORKDIR /benchmarking
RUN wget https://github.com/rsuchecki/biokanga_benchmark/archive/v${EVALU_VERSION}.tar.gz \
  && tar xzvf v${EVALU_VERSION}.tar.gz \
  && mv biokanga_benchmark-${EVALU_VERSION}/* . \
  && rmdir biokanga_benchmark-${EVALU_VERSION}

ENV PATH "$PATH:/benchmarking"