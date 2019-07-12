#!/usr/bin/env bash

minimap2 \
  -x sr \
  -I 150G \
  -t ${task.cpus} \
  -d ${ref}.mmi \
  ${ref}