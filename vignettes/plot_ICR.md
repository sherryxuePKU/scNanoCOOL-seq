## Table of contents

- [Quick start](#quick-start)
- [What's included](#whats-included)
- [Usage](#Usage)
- [Last modification](#last-modification)
- [Contact](#contact)

## Quick start

the folder /bin/ICR/ is used to plot read-level methylation of imprinting controled regions (ICRs) in the following steps. It is worth mentioning that the main script is modified from the package CluBCpG (https://github.com/waterlandlab/CluBCpG) and adapted to NOMe-seq.

- Step1: merge the .bam files of single cells to pseudo-bulk
- Step2: generate the methylation matrix of WCG sites and plot the matrix
- Step3: de novo find SNV in the given regions and plot 

## What's included

```text
ICR/
├── README.md
├── plotReadLevel.sh
└── script/
    ├── PlotTanghulu.py
    ├── PlotTanghulu.SNV.py
    └── PlotTanghulu.R
```

## Usage
Step1: merge the .bam files
Align the .bam files as shown in the NDR folder

Common preset

    BASEDIR=/path/to/thisfolder
    FUNCTION=${BASEDIR}/plotReadLevel.sh
    CONFIG=scNanoCool_config

Step2: generate the WCG matrix and plot the result

    source $FUNCTION $work_dir $psbulk_sample $CONFIG
    do_plot_tanghulu.ICR

Step2: generate the WCG matrix and plot the result

    source $FUNCTION $work_dir $psbulk_sample $CONFIG
    do_plot_tanghulu.SNV $region ICR

## Last modification

2022/10/16

## Last modification

xxhui@stu.pku.edu.cn
