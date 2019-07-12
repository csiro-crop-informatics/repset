#!/usr/bin/env bash

gsnap \
  -D genomeDir \
  -d ${idxmeta.target} \
  -A sam \
  --nthreads ${task.cpus} \
  ${ALIGN_PARAMS} \
	${r1} ${r2} \
  > sam