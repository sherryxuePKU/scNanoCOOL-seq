## Table of contents

- [Quick start](#quick-start)
- [What's included](#whats-included)
- [Usage](#Usage)
- [Last modification](#last-modification)
- [Contact](#contact)

## Quick start

the folder /bin/CNV/ is used to call copy number variations (CNVs) by [FREEC-11.5](https://github.com/BoevaLab/FREEC) with following steps

- Step1: prepare reference as input of FREEC -> $chr_len & $freec_db
- Step2: calc GC-ratio normalized signal in each bin and infer the copy number -> .bam_ratio.txt
- Step3: plot the result ( heatmap & dotplot )


## What's included

```text
CNV/
├── README.md
├── callCNV.sh
├── plot_dotplot.R
└── plot_heatmap.R

```

## Usage

Step1: prepare reference

Follow the manual of [FREEC](http://boevalab.inf.ethz.ch/FREEC/)

Step2: calc the signal -> .bam_ratio.txt

    ## for single cell
    BASEDIR=/path/to/this folder
    FUNCTION=${BASEDIR}/callCNV.sh
    CONFIG=scNanoCool_config

    source $FUNCTION $work_dir $sc_sample $CONFIG
    do_CNV 1
    do_CNV 5
    do_CNV 10

Step3: plot the result

Run plot_dotplot.R and plot_heatmap.R in the Rstudio

## Last modification

2022/10/15

## Last modification

xxhui@stu.pku.edu.cn
