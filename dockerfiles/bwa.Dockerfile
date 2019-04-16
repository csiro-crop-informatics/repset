FROM ubuntu:16.04

LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"

RUN apt-get update && apt-get install --yes build-essential gcc-multilib apt-utils zlib1g-dev wget gawk

# Get source code from git
RUN apt-get install --yes git
WORKDIR /tmp
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /tmp/bwa
RUN git checkout v0.7.17

# Compile
RUN make
RUN cp -p bwa /usr/local/bin

# Cleanup
RUN rm -rf /tmp/bwa \
  && apt-get clean \
  && apt-get remove --yes --purge build-essential gcc-multilib apt-utils zlib1g-dev wget