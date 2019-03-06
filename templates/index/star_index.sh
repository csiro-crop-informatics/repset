#!/usr/bin/env bash

mkdir -p genomeDir
STAR --runThreadN ${task.cpus} \
  --runMode genomeGenerate \
  --genomeDir genomeDir \
  --genomeFastaFiles ${ref}
mv Log.out indexing_Log.out #otherwise file overwritten by align process which in turn triggers re-runs of align due to cache mismatch
