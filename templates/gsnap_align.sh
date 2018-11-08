#!/usr/bin/env bash

gsnap \
  -D genomeDir \
  -d ${idxmeta.target} \
  -A sam \
  --merge-distant-samechr \
	--novelsplicing 1 \
  --max-mismatches 0.05 \
  --pairmax-rna 100000 \
  --localsplicedist 100000 \
  --adapter-strip paired \
  --nthreads ${task.cpus} \
  --batch 5 \
	${r1} ${r2} \
  > sam