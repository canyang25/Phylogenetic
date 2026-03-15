#!/usr/bin/env bash
set -e

mkdir -p results/trees/iqtree

iqtree \
  -s results/msa/M.aln.fa \
  -m MFP \
  -B 1000 \
  --bnni \
  -alrt 1000 \
  -nt AUTO \
  --prefix results/trees/iqtree/M
