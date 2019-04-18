FROM ubuntu:18.04
LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN apt-get update \
  && apt-get install -y wget libgomp1 build-essential gcc-multilib apt-utils zlib1g-dev


#overwrite at build time with e.g. --build-arg VERSION=2.6.1
ARG VERSION=2.7.0f

#Install star
WORKDIR /star
RUN wget https://github.com/alexdobin/STAR/archive/${VERSION}.tar.gz \
  && tar -xzf ${VERSION}.tar.gz \
  && cd STAR-${VERSION}/source \
  && make

# Cleanup
RUN rm /star/${VERSION}.tar.gz \
  && apt-get clean \
  && apt-get remove --yes --purge wget build-essential gcc-multilib apt-utils zlib1g-dev

ENV PATH "$PATH:/star/STAR-${VERSION}/bin/Linux_x86_64/"
