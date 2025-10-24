
## Overview

scNanoCOOL-seq is a sequencing method based on nanopore sequencing platform, which involves **bisulfite conversion** and adopt a **single-barcode PCR amplification** strategy. Up to now, its DNA fragments is approximately **900 bp** in length. 
As shown below, the experimental process of scNanoCOOL-seq mainly includes: 
1) **Preprocess**: cell lysis in a gentle way and in vitro methylation by DNA transferase M.CviPI enzyme, followed by nucleocytoplasmic separation of single cells;
2) **DNA methylome**: the nucleus of single cells is processed with single-linker amplification after bisulfite treatment and PBAT-like library construction and is used for single-molecule sequencing on the nanopore platform;
3) **RNA transcriptome**: the cytoplasmic part of the single cell is used for scRNA-seq library construction (optimized STRT-seq). 

<p align="center"> 
<img src="protocol.png" style="width: 90%; height: 90%"/>​
</p>

## Content

* `bin/`: functions used in each section from pre-processing to downstream analysis
* `reference/`: necessary files as the reference
* `vignettes/`: brief introductions for each section

## Data

The raw data is accessible at SRA BioProject [PRJNA905717](https://www.ncbi.nlm.nih.gov/sra/PRJNA905717).  
Processed data is provided upon request.

## Citation

If you use scripts or code in your research, please cite:

```
Lin J, Xue X, Wang Y, et al. scNanoCOOL-seq: a long-read single-cell sequencing method for multi-omics profiling within individual cells. Cell Res. 2023;33(11):879-882. doi:10.1038/s41422-023-00873-5
```

## Contact

Xiaohui Xue (xxhui@stu.pku.edu.cn) 

