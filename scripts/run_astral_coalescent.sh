#!/usr/bin/env bash

# Coalescent method: ASTRAL-IV
# Builds one ML gene tree per influenza segment with IQ-TREE, then summarizes
# them into a species tree using quartet-based coalescent inference.
#
# Assumptions: segments are independent loci, alignment is homologous,
# gene-tree discordance reflects ILS.
# Limitations: reassortment violates the coalescent model; only 8 loci.

mkdir -p results/trees/coalescent/gene_trees

# One ML tree per segment
for seg in PB2 PB1 PA HA NP NA M NS; do
    iqtree -s results/msa/${seg}.aln.fa -m MFP -nt 4 \
        --prefix results/trees/coalescent/gene_trees/${seg}
done

# Combine gene trees
cat results/trees/coalescent/gene_trees/{PB2,PB1,PA,HA,NP,NA,M,NS}.treefile \
    > results/trees/coalescent/influenza_segments.gene_trees.nwk

# Run ASTRAL-IV
astral4 -R \
    -i results/trees/coalescent/influenza_segments.gene_trees.nwk \
    -o results/trees/coalescent/influenza_segments.astral4.tree \
    2> results/trees/coalescent/influenza_segments.astral4.log
