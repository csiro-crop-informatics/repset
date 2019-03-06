#!/usr/bin/env bash

bwa mem \
  -t ${task.cpus} \
  -L 1 \
  ${idxmeta.target} \
  ${r1} <(sed 's/b\$/a/' ${r2}) \
  | gawk -vOFS="\\t" '\$1 !~ /^@/ && and(\$2,128) {sub(/a\$/,"b",\$1)};{print}' \
  > sam
