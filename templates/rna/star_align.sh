#!/usr/bin/env bash

STAR \
  --runThreadN ${task.cpus} \
  --genomeDir genomeDir \
  --readFilesIn ${r1} ${r2} \
  ${ALIGN_PARAMS}



