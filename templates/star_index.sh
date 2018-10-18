#!/usr/bin/env bash

mkdir -p genomeDir
STAR --runThreadN ${task.cpus} \
  --runMode genomeGenerate \
  --genomeDir genomeDir \
  --genomeFastaFiles ${ref}

