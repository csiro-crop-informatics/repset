#!/usr/bin/env bash

bwa mem \
  -t ${task.cpus} \
  ${ALIGN_PARAMS} \
  ${ref} \
  ${reads} \
  | samtools view -b \
  > out.bam
