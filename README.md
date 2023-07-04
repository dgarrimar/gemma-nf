# gemma-nf

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.1-blue.svg)](http://nextflow.io)

A pipeline for multi-trait genome-wide association studies (GWAS) using [GEMMA](https://github.com/genetics-statistics/GEMMA).

The pipeline performs the following analysis steps:

* Split genotype file 
* Preprocess genotype and phenotype data
* Compute kinship
* Test for association between phenotypes and genetic variants
* Collect summary statistics

The pipeline uses [Nextflow](http://www.nextflow.io) as the execution backend. Please check [Nextflow documentation](http://www.nextflow.io/docs/latest/index.html) for more information.

## Requirements

- Unix-like operating system (Linux, MacOS, etc.)
- Java 8 or later 
- [Docker](https://www.docker.com/) (v1.10.0 or later) or [Singularity](http://singularity.lbl.gov) (v2.5.0 or later)

## Quickstart (~2 min)

1. Install Nextflow:
    ```
    curl -fsSL get.nextflow.io | bash
    ```

2. Make a test run:
    ```
    nextflow run dgarrimar/gemma-nf -with-docker
    ```

**Notes**: move the `nextflow` executable to a directory in your `$PATH`. Set `-with-singularity` to use Singularity instead of Docker. 

(*) Alternatively you can clone this repository:
```
git clone https://github.com/dgarrimar/gemma-nf
cd gemma-nf
nextflow run gemma.nf -with-docker
```

## Pipeline usage

Launching the pipeline with the `--help` parameter shows the help message:

```
nextflow run gemma.nf --help
```

```
N E X T F L O W  ~  version 22.04.0
Launching `gemma.nf` [modest_ritchie] DSL1 - revision: a8f7d96105

G E M M A - N F
============================================================================================
Performs multi-trait GWAS using GEMMA (https://github.com/genetics-statistics/GEMMA)

Usage:
    nextflow run gemma.nf [options]

Parameters:
 --geno GENOTYPES            genotype data in VCF format, indexed
 --pheno PHENOTYPES          (covariate-adjusted) phenotype data
 --maf MAF                   MAF filter (default: 0.01
 --l VARIANTS/CHUNK          number of variants tested per chunk (default: 500)
 --t THREADS                 number of threads (deafault: 1)
 --dir DIRECTORY             output directory (default: result)
 --out OUTPUT                output file prefix (default: gemma.tsv)
```
