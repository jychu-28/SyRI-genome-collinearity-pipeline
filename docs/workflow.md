# Workflow Details

## Overview

This pipeline is designed for large chromosome-level genomes. Direct whole-chromosome alignment can be slow and memory-intensive, so each chromosome is extracted and split into 100 Mb chunks. Neighboring chunks are aligned, merged, corrected back to chromosome-level coordinates, and then used for SyRI structural variation detection.

## Step 1: Chromosome Extraction and Chunk Splitting

Input whole-genome FASTA files contain chromosomes 1-5. Each target chromosome is extracted into an individual FASTA file, then split into chunks:

```text
Chr1.fa
Chr1_chunk0000.fa
Chr1_chunk0001.fa
...
```

Default settings:

```text
chunk size: 100 Mb
overlap:    10 kb
```

## Step 2: Chunk Alignment

Adjacent chunk pairs are aligned using:

```bash
minimap2 -x asm20 --cs --eqx
```

The pipeline only aligns query chunks around the expected corresponding reference chunk, reducing unnecessary all-vs-all comparisons.

## Step 3: PAF Coordinate Correction

minimap2 reports coordinates relative to each chunk. `fix_paf_coords.py` restores the original chromosome coordinates:

```text
chromosome_coordinate = chunk_coordinate + chunk_index * chunk_size
```

It also replaces chunk lengths with real chromosome lengths from `.fai` files.

## Step 4: CIGAR Tag Conversion

SyRI requires a `cg:Z:` CIGAR tag when using PAF input. minimap2 can produce a `cs:Z:` tag with `--cs`, so `add_cigar_to_paf.py` converts:

```text
cs:Z::100*ag:50+actg:20
```

into:

```text
cg:Z:100=1X50=4I20=
```

## Step 5: SyRI Structural Variation Detection

SyRI is run with PAF input:

```bash
syri -c merged.fixed.cg.paf -r ref.fa -q query.fa -F P --nc 1 --nosnp
```

The main output is:

```text
syri.out
syri.summary
syri.vcf
```

## Step 6: plotsr Visualization

SyRI outputs from all five chromosome pairs are merged and visualized with plotsr:

```bash
plotsr --sr merged.syri.out --genomes genomes.txt -o comparison.png
```

The final PNG shows genome collinearity and rearrangement patterns.

