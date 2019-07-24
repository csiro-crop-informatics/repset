#!/usr/bin/env bash

rapmap quasimap \
  --selAln \
  -t ${task.cpus} \
  ${ALIGN_PARAMS} \
  --index ${ref}.idx \
  -1 ${reads[0]} \
  -2 ${reads[1]} \
  --writeUnmapped \
  --output out.sam

