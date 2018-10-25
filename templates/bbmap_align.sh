#!/usr/bin/env bash

bbmap.sh \
  in=${r1} \
  in2=${r2} \
  threads=${task.cpus} \
  maxindel=100000 \
  ambiguous=best \
  intronlen=20 \
  keepnames=t \
  sam=1.3 \
  minid=95 \
  local=t \
  out=sam