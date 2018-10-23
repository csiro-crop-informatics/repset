#!/usr/bin/env bash

gsnap \
  -D genomeDir \
  -d ${idxmeta.target} \
  -A sam \
  --merge-distant-samechr \
	--novelsplicing 1 \
  --nthreads ${task.cpus} \
  --batch 5 \
	${r1} ${r2} \
  > sam