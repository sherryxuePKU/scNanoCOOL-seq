## Table of contents

- [Quick start](#quick-start)
- [What's included](#whats-included)
- [Usage](#usage)
- [Last modification](#last-modification)
- [Contact](#contact)

## Quick start

the folder /bin/PCA/ is used to generate matrix of given region for WCG/GCH methylation and contains following steps

- Step1: calc methylation level in single-cell level -> .bed
- Step2: merge all the .bed to matrix -> .mtx
- Step3: perform PCA analysis on the matrix and plot the result


## What's included

```text
reduce_dimension/
├── README.md
├── generateMatirx.sh
└── reduce_dimension.R
```

## Usage
Common preset

    BASEDIR=/path/to/this folder
    FUNCTION=${BASEDIR}/generateMatrix.sh
    CONFIG=scNanoCool_config

Step1: calc methylation level -> .bed

    ## for single cells
    source $FUNCTION $work_dir $sc_sample $CONFIG
    do_element_meth 3 GCH TSS_0.75k
    do_element_meth 3 WCG TSS_0.75k

Step2: merge to matrix -> .mtx

    ## for batch
    source $FUNCTION $work_dir $batch_sample $CONFIG
    do_generate_matrix GCH TSS_0.75k 0.7 nanoCool_5aza
    do_generate_matrix WCG TSS_0.75k 0.7 nanoCool_5aza

Step2: reduce dimension & plot

Run reduce_dimension.R in the Rstudio

## Last modification

2022/10/15

## Contact

xxhui@stu.pku.edu.cn
