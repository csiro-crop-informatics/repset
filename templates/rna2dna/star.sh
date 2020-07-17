#!/usr/bin/env bash

STAR \
  --runThreadN ${task.cpus} \
  --genomeDir genomeDir \
  --readFilesIn ${reads} \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  ${ALIGN_PARAMS}