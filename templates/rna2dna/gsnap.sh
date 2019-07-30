#!/usr/bin/env bash

gsnap \
  -D genomeDir \
  -d GENOME \
  -A sam \
  --nthreads ${task.cpus} \
  ${ALIGN_PARAMS} \
  <(zcat ${reads[0]}) \
  <(zcat ${reads[1]}) \
  > out.sam