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
