#!/usr/bin/env bash

whisper \
  -t ${task.cpus} \
  -dist_paired 100000 \
  -out sam \
  -temp \${TMPDIR}/whisper_tmp_
  ${idxmeta.target} \
  ${r1} \
  ${r2}