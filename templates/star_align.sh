#!/usr/bin/env bash

STAR \
  --runThreadN ${task.cpus} \
  --genomeDir ./ \
  --readFilesIn ${r1} ${r2}

