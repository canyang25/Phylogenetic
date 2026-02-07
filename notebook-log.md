# Notebook Log

## Dataset
I’m using a **public viral sequencing dataset from NCBI SRA** for this semester. These are **raw reads** (SRA runs → FASTQ). The 10 SRR IDs are in `data/metadata/srr_ids.txt`.

## Feb 7, 2026 — Sequencing HW: QC with FastQC
For this homework I downloaded the reads with `fastq-dump` and checked read quality using **FastQC**. I saved FASTQs to `data/raw/` and FastQC reports to `data/qc/fastqc/`.

Commands I used:

```bash
fastq-dump --split-files -O data/raw <SRR_ID>
fastqc -o data/qc/fastqc -t 4 data/raw/*.fastq
```
