#!/usr/bin/env bash

bwa mem \
  -t ${task.cpus} \
  ${ALIGN_PARAMS} \
  ${ref} \
  ${reads} \
  > out.sam
