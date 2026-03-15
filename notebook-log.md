# Notebook Log

## Dataset(Feb.2)
I’m using a **public viral sequencing dataset from NCBI SRA** for this semester. These are **raw reads** (SRA runs → FASTQ). The 10 SRR IDs are in `data/metadata/srr_ids.txt`.

## Feb 7, 2026 — Sequencing HW: QC with FastQC
For this homework I downloaded the reads with `fastq-dump` and checked read quality using **FastQC**. I saved FASTQs to `data/raw/` and FastQC reports to `data/qc/fastqc/`.

Commands I used:

```bash
fastq-dump --split-files -O data/raw <SRR_ID>
fastqc -o data/qc/fastqc -t 4 data/raw/*.fastq
```


# HW MSA Notes(Feb.17)

## Chosen MSA Method
- Method: **MAFFT G-INS-i**
- Command mode: `--globalpair --maxiterate 1000 --reorder`

## Why This Method
- My dataset is Influenza A, so I analyzed it segment by segment (`PB2 PB1 PA HA NP NA M NS`), which makes it reasonable to align the same segment across samples.
- I used MAFFT in a higher-accuracy mode because I only had a small number of samples.
- Before running MAFFT, I built one consensus sequence per sample per segment (using reference-guided mapping and masking low-depth bases), then aligned those consensus sequences.

## Algorithm Description
MAFFT does multiple sequence alignment by first building a progressive alignment and then improving it with iterative refinement steps.

## Assumptions
- The same segment (for example `PB2`) is homologous across all samples.
- There are no major rearrangements within each segment in this dataset.

## Limitations
- Because consensus was reference-guided, some bias toward the reference sequence is possible.
- Low-depth positions are masked as `N`, which means some information is missing.
- Influenza can reassort segments, so interpretation across segments can be more complicated.
- Any within-sample variation (or mixed infection) is simplified into one consensus sequence.

## Feb 27, 2026 — Distance and Parsimony HW

Used the M segment alignment (`results/msa/M.aln.fa`). Built two trees in R using `ape` and `phangorn`. Script: `scripts/distance_parsimony.R`. Output trees in `results/trees/`.

**Neighbor Joining:** computes pairwise distances with TN93 model, then builds a tree by joining the closest pair iteratively. Assumes distances are additive. Limitation: heuristic, not guaranteed optimal.

**Parsimony:** finds the tree requiring the fewest mutations. No substitution model needed. Limitation: prone to long-branch attraction; heuristic search can get stuck locally.

## Mar 14, 2026 — Maximum Likelihood HW with IQ-TREE

I chose **IQ-TREE** as my maximum likelihood method. I used the same M segment alignment (`results/msa/M.aln.fa`) so I can compare the ML tree to my earlier neighbor joining and parsimony trees.

## Why This Method
- IQ-TREE is a maximum likelihood phylogeny program.
- It can test substitution models automatically with `-m MFP`.
- It is commonly used and has a simple command line workflow.

## Algorithm Description
IQ-TREE searches for the tree that gives the highest likelihood for the aligned sequences under a chosen substitution model. It uses heuristic tree search, and it can also estimate branch support with ultrafast bootstrap and SH-aLRT.

## Assumptions
- The sequences in the alignment are homologous.
- The M segment alignment is correct enough for tree inference.
- The chosen substitution model is a reasonable fit for the data.
- One tree can represent the history of this segment.

## Limitations
- Maximum likelihood is still model-based, so a poor model can affect the tree.
- My sequences are consensus sequences from reference-guided mapping, so some reference bias is possible.
- Low-depth sites were masked with `N`, so some information is missing.
- Influenza can reassort between segments, so one segment tree does not always represent the whole virus history.

Commands I used:

```bash
conda env update -f envs/phylo-msa.yml
bash scripts/run_iqtree.sh
```

Main IQ-TREE command:

```bash
iqtree2 \
  -s results/msa/M.aln.fa \
  -m MFP \
  -B 1000 \
  --bnni \
  -alrt 1000 \
  -nt AUTO \
  --prefix results/trees/iqtree/M
```
