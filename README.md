# aScanMeth: Allele-Specific Methylation Analysis from Oxford Nanopore Sequencing Data

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](link-to-your-build-status)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![GitHub Stars](https://img.shields.io/github/stars/your_username/aScanMeth.svg?style=social)](https://github.com/your_username/aScanMeth)
[![Galaxy Tool](https://img.shields.io/badge/Galaxy%20Tool-Available-purple)](link-to-galaxy-tool-shed)

---

## Background

**DNA methylation** is a key epigenetic mechanism that modulates gene expression, with profound implications for human health and disease. **Allele-specific methylation (ASM)**, where methylation levels differ between maternal and paternal haplotypes, has been implicated in a range of conditions, including cancer and developmental disorders. Recent advances in sequencing technologies, particularly single-molecule sequencing, now enable haplotype-resolved genome assemblies and the simultaneous detection of multiple epigenetic modifications directly on DNA and RNA.

Notwithstanding the potential implications and applications of allele-specific methylation in precision medicine and/or for advancing our molecular knowledge of biological systems, bioinformatics methods for an integrated and streamlined characterization of ASM from single-molecule sequencing data are currently lacking.

---

## Main Findings

Here we present **aScanMeth**, a fully automated and customizable computational workflow for the identification of ASM from haplotype-resolved raw Oxford Nanopore sequencing data. aScanMeth automates the identification and phasing of genetic variants into discrete haplotype blocks and enables the reliable detection of differentially methylated alleles. The method is highly effective, achieving performance comparable to—or exceeding—that of existing approaches. The workflow is implemented as a standalone tool and also provided in the form of a modular and user-friendly Galaxy pipeline, ensuring reproducibility and seamless integration into any instance of the popular bioinformatics workflow manager Galaxy.

---

## Biological Relevance

By enabling the efficient and reproducible inference of haplotype-specific methylation patterns, aScanMeth offers a valuable resource for epigenetic research. Its potential applications range from identifying imprinted regions in the genome to advancing our molecular understanding of the implications of methylation in the modulation of gene expression, and in health and disease.

---

## Significance

aScanMeth enables accurate, haplotype-resolved detection of allele-specific DNA methylation from Oxford Nanopore sequencing data, providing a powerful and user-friendly tool for advancing epigenomic research and precision medicine.
