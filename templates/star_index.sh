#!/usr/bin/env bash

STAR --runThreadN ${task.cpus} \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles ${ref}

