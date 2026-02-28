#Distance and Parsimony Methods

# NJ  with TN93 distances:
# Builds a tree by joining the closest pair of taxa iteratively.
# Assumes distances are additive and the TN93 model fits the data.
# Limitation: heuristic, not guaranteed optimal and sensitive to distance errors.

# Parsimony:
# Finds the tree with the fewest total character changes.
# Assumes fewer changes = better explanation, no substitution model needed.
# Limitation: prone to long-branch attraction; heuristic search can get stuck.

# install.packages("adegenet", dep=TRUE)
# install.packages("phangorn", dep=TRUE)

library(ape)
library(adegenet)
library(phangorn)

dna <- fasta2DNAbin(file="results/msa/M.aln.fa")
dna2 <- as.phyDat(dna)

dir.create("results/trees", showWarnings=FALSE, recursive=TRUE)

# Distance: NJ tree
D <- dist.dna(dna, model="TN93", pairwise.deletion=TRUE)
tre <- njs(D)
tre <- ladderize(tre)

pdf("results/trees/nj_tree.pdf")
plot(tre, cex=.6)
title("NJ tree")
dev.off()

# Parsimony: MP tree
tre.ini <- njs(dist.dna(dna, model="raw", pairwise.deletion=TRUE))
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2)

pdf("results/trees/parsimony_tree.pdf")
plot(tre.pars, cex=0.6)
title("Most Parsimonious Tree")
dev.off()
