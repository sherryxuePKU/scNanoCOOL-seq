scNanoCOOL-seq
=========

Source code for manuscript entitles ***scNanoCOOL-seq: a long-read single-cell sequencing method for multi-omics profiling within individual cells***.

## Overview
--------
scNanoCOOL-seq technology mainly uses bisulfite conversion and single-barcode PCR amplification strategies and finally obtained DNA fragments of approximately 900 bp in length. 
As shown below, the experimental process of scNanoCOOL-seq mainly includes: 
1) cell lysis in a gentle way and in vitro methylation by DNA transferase M.CviPI enzyme, followed by nucleocytoplasmic separation of single cells;
2) the nucleus of single cells is processed with single-linker amplification after bisulfite treatment and PBAT-like library construction and is used for single-molecule sequencing on the Nanopore platform;
3) the cytoplasmic part of the single cell is used for scRNA-seq library construction (optimized STRT-seq). 

<p align="center"> 
<img src="insert_pic2.png" style="width: 100%; height: 100%"/>​
</p>

## Content
-------
* `bin/`: functions used in each section from pre-processing to downstream analysis
* `data/`: example data provided to reproduce the result
* `reference/`: necessary files as the reference
* `vignettes/`: brief introductions for each section

## Data
-------
The raw data is accessible at SRA BioProject [PRJNA905717](https://www.ncbi.nlm.nih.gov/sra/PRJNA905717). Parsed data is provided upon request.

## Citation
--------
Please cite [PMID: 37700167](https://www.nature.com/articles/s41422-023-00873-5)

## Contact
-------
* Computational analysis: Xiaohui Xue (xxhui@stu.pku.edu.cn) 

