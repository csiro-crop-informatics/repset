#!/usr/bin/env bash

#stripping FASTA description line otherwise gsnap falls over at mapping
sed '/^>/s/ .*//' ${ref} > _REF_ \
&& gmap_build -D genomeDir -d GENOME _REF_
rm _REF_
