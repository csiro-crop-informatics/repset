#!/usr/bin/env bash

bwa mem \
  -t ${task.cpus} \
  ${idxmeta.target} \
  1.fq.gz 2.fq.gz \
  > out.sam
