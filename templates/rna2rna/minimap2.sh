#!/usr/bin/env bash

minimap2 \
  -ax sr \
  -t ${task.cpus} \
  ${ALIGN_PARAMS} \
  ${ref}.mmi \
  ${reads} \
  | samtools view -b \
  > out.bam