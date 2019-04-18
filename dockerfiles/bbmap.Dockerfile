FROM java:openjdk-8u111-jre

LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
#SHELL ["/bin/bash", "-c"]

RUN apt-get update

#overwrite at build time with e.g. --build-arg VERSION=38.26
ARG VERSION=38.44

#Install dart
WORKDIR /
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_${VERSION}.tar.gz \
  && tar xzvf BBMap_${VERSION}.tar.gz && rm BBMap_${VERSION}.tar.gz
ENV PATH "${PATH}:/bbmap"