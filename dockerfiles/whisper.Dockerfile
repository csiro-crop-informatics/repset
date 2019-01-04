FROM ubuntu:16.04

RUN apt-get update --fix-missing \
    && apt-get install -y wget \
    && wget -O /bin/whisper https://github.com/refresh-bio/Whisper/releases/download/v1.1/whisper \
    && wget -O /bin/whisper-index https://github.com/refresh-bio/Whisper/releases/download/v1.1/whisper-index \
    && chmod +x /bin/whisper /bin/whisper-index \
    && apt-get remove -y wget \
    && apt-get autoremove -y \
    && apt-get clean