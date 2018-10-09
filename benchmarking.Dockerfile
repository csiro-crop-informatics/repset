FROM ubuntu:16.04
LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y wget build-essential automake make unzip pigz ruby-full python2.7 samtools vim less

# Setup evaluation framework
#overwrite at build time with e.g. --build-arg EVALU_VERSION=0.2
ARG EVALU_VERSION=0.4
WORKDIR /project/itmatlab/aligner_benchmark
RUN wget https://github.com/rsuchecki/biokanga_benchmark/archive/v${EVALU_VERSION}.tar.gz \
  && tar xzvf v${EVALU_VERSION}.tar.gz \
  && mv biokanga_benchmark-${EVALU_VERSION}/* . \
  && rmdir biokanga_benchmark-${EVALU_VERSION} \
  && gem install erubis \
  && mkdir -p \
      /project/itmatlab/aligner_benchmark/dataset/human/{genome,dataset_t{1,2,3}r{1,2,3}} \
      /project/itmatlab/aligner_benchmark/tool_results \
      /project/itmatlab/aligner_benchmark/statistics/
