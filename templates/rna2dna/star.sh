#!/usr/bin/env bash

STAR \
  --runThreadN ${task.cpus} \
  --genomeDir genomeDir \
  --readFilesIn ${reads} \
  --readFilesCommand zcat \
  ${ALIGN_PARAMS}



