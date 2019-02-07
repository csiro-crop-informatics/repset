#!/usr/bin/env bash

STAR \
  --runThreadN ${task.cpus} \
  --genomeDir genomeDir \
  --readFilesIn ${r1} ${r2} \
  --outFilterMismatchNoverReadLmax 0.05



