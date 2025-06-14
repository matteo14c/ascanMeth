# aScanMeth: Allele-Specific Methylation Analysis from Oxford Nanopore Sequencing Data

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](link-to-your-build-status)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![GitHub Stars](https://img.shields.io/github/stars/matteo14c/aScanMeth.svg?style=social)](https://github.com/matteo14c/aScanMeth)
[![Galaxy Tool](https://img.shields.io/badge/Galaxy%20Tool-Available-purple)](link-to-galaxy-tool-shed)

---

## Background

**DNA methylation** is a key epigenetic mechanism that modulates gene expression, with profound implications for human health and disease. **Allele-specific methylation (ASM)**, where methylation levels differ between maternal and paternal haplotypes, has been implicated in a range of conditions, including cancer and developmental disorders. Recent advances in sequencing technologies, particularly single-molecule sequencing, now enable haplotype-resolved genome assemblies and the simultaneous detection of multiple epigenetic modifications directly on DNA and RNA.

Notwithstanding the potential implications and applications of allele-specific methylation in precision medicine and/or for advancing our molecular knowledge of biological systems, bioinformatics methods for an integrated and streamlined characterization of ASM from single-molecule sequencing data are currently lacking.

---

## aScanMeth

**aScanMeth**, is a fully automated and customizable computational workflow for the identification of ASM from haplotype-resolved raw Oxford Nanopore sequencing data. 
aScanMeth automates the identification and phasing of genetic variants into discrete haplotype blocks and enables the reliable detection of differentially methylated alleles. 
The method is highly implemented in the form of:

1. A standalone command line tool and workflow.
2. A fully automated Galaxy workflow*. **⚠️ Note:** This feature is not currently implemented, but will be implemented soon **⚠️**.


Please see below for additional details.


---

## Dependencies

aScanMeth depends on the following software:

1. Clair3 for variant calling. see: https://github.com/HKU-BAL/Clair3 .
2. modkit for handling modBam files, see: https://github.com/nanoporetech/modkit.
3. whatshap for haplotype phasing, see:  https://github.com/whatshap/whatshap
4. R and Rscript: https://cran.r-project.org/.

All the dependencies need to be installed and properly configured 

---

## Installation
Install the dependencies, and clone the main git repository of **aScanMeth**

---

## Usage

### Galaxy workflow

A test environment will be made available soon. **⚠️Note:** This feature is not currently implemented, but will be implemented soon **⚠️**.

### Command line

Obtain a phased modBam file by using Clair3, whatshap and modkit. See the file XXX for an example worfklow. **⚠️Note:** This feature is not currently implemented, but will be implemented soon **⚠️**. As a result you should obtain 2 distinct bedmethyl files, one for each haplotype.
Please see  https://github.com/nanoporetech/modkit for a more detailed explanation of the bedmethyl format.

To execute a differential methylation analysis run:
```R
Rscript aScanMeth.R --h1 1.bed --h2 2.bed --n 20
```
Where 1.bed and 2.bed are haplotype resolved bedmethyl files as obtained from modkit and n is the number of consecutive CpGs tested.   

### Galaxy workflow

A standalone galaxy instance will be made available to test the galaxy workflow soon.
This feature is not currently implemented.

---

## Output

The current version of **aScaMeth** produces 3 files as its main output:

1. results_D.CpG.csv: results for differential density of CpG. This one indicates/reports whether the paternal/maternal haplotypes has an increased number of CpGs w.r.t the other haplotype;
2. results_H.CpG.csv: results for differential *5hmC* methylation;
3. results_M.CpG.csv: results for differential *5mC* methylation. 

An example can be found in the folder [exampleResults](https://github.com/matteo14c/ascanMeth/tree/main/exampleResults).
The format is a csv file with 10 columns, and example is provided below:
| chr | start | end | h1 | h2 | h1stat | h2stat | aScanstat | pvalue | FDR |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1000 | 2000 | mCpG H1 | mCpG H1 | 0.85 | 0.92 | 0.78 | 0.001 | 0.015 |

The columns:

* chr, start, end indicate genomics coordinates;
* h1, and h2 report the number of CpG (results_D.CpG.csv) or the number of modified CpGs (results_H.CpG.csv or results_M.CpG.csv) for the two haplotypes;
* FDR reports the FDR-corrected p-value for differential methylation/CpG density.





