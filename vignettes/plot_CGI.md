## Table of contents

- [Quick start](#quick-start)
- [What's included](#whats-included)
- [Usage](#Usage)
- [Last modification](#last-modification)
- [Contact](#contact)

## Quick start

the folder /bin/CGI/ is used to plot read-level methylation of CpG islands (CGIs) and their flanking regions in the following steps. It is worth mentioning that the main script is modified from the package [CluBCpG](https://github.com/waterlandlab/CluBCpG) and adapted to NOMe-seq.

- Step1: merge the .bam files of single cells to pseudo-bulk
- Step2: generate the methylation matrix of WCG sites and plot the matrix marked the boundary of CGIs

## What's included

```text
CGI/
├── README.md
├── plotReadLevel_CGI.sh
└── script/
    ├── PlotTanghulu.py
    └── PlotTanghulu_CGI.R
```

## Usage
Step1: merge the .bam files
Align the .bam files as shown in the NDR folder

Step2: generate the WCG matrix and plot the result

    BASEDIR=/path/to/thisfolder
    FUNCTION=${BASEDIR}/plotReadLevel_CGI.sh
    CONFIG=scNanoCool_config

    source $FUNCTION $work_dir $psbulk_sample $CONFIG
    do_plot_tanghulu.CGI

## Last modification

2022/10/15

## Contact

xxhui@stu.pku.edu.cn
