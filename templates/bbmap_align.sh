#!/usr/bin/env bash

bbmap.sh \
  in=${r1} \
  in2=${r2} \
  threads=${task.cpus} \
  maxindel=200000 \
  ambiguous=best \
  intronlen=20 \
  keepnames=t \
  sam=1.3 \
  out=sam