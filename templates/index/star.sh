#!/usr/bin/env bash

mkdir -p genomeDir
STAR --runThreadN ${task.cpus} \
  --runMode genomeGenerate \
  --genomeDir genomeDir \
  --genomeFastaFiles ${ref}

#fixed and no longer needed in STAR approx > 2.7.0f <= 2.7.5a . Previously file overwritten by align process which in turn triggers re-runs of align process due to cache mismatch
[ -f Log.out ] && mv Log.out indexing_Log.out || echo "no Log.out to cleanup"