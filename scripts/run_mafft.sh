#!/usr/bin/env bash
set -e

mkdir -p results/msa

mafft --globalpair --maxiterate 1000 --reorder --thread 8 results/msa/input/PB2.fa > results/msa/PB2.aln.fa
mafft --globalpair --maxiterate 1000 --reorder --thread 8 results/msa/input/PB1.fa > results/msa/PB1.aln.fa
mafft --globalpair --maxiterate 1000 --reorder --thread 8 results/msa/input/PA.fa  > results/msa/PA.aln.fa
mafft --globalpair --maxiterate 1000 --reorder --thread 8 results/msa/input/HA.fa  > results/msa/HA.aln.fa
mafft --globalpair --maxiterate 1000 --reorder --thread 8 results/msa/input/NP.fa  > results/msa/NP.aln.fa
mafft --globalpair --maxiterate 1000 --reorder --thread 8 results/msa/input/NA.fa  > results/msa/NA.aln.fa
mafft --globalpair --maxiterate 1000 --reorder --thread 8 results/msa/input/M.fa   > results/msa/M.aln.fa
mafft --globalpair --maxiterate 1000 --reorder --thread 8 results/msa/input/NS.fa  > results/msa/NS.aln.fa
