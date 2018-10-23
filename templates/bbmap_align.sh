#!/usr/bin/env bash

bbmap.sh \
  in=${r1} \
  in2=${r2} \
  threads=${task.cpus} \
  maxindel=200000 \
  ambig=random \
  intronlen=20 \
  keepnames=t \
  out=sam