## Table of contents

- [Quick start](#quick-start)
- [What's included](#whats-included)
- [Usage](#Usage)
- [Last modification](#last-modification)
- [Contact](#contact)

## Quick start

the folder /bin/Preprocess/ is used to preprocess the raw data of scNanoCool-seq to the methylation information of WCG and GCH site with the following steps

- Step1: demultiplex the raw data to single-cell .fastq files
- Step2: trim the barcode and random primer 
- Step3: align to assembly hg38 by bismark-MM2
- Step4: filter the aligned reads by the read-level methylation level of CH sites
- Step5: calc the methylation level of WCG and GCH site
- Step6: generate the meta-statistics for quality control   

_To be mentioned, we used the development version of [bismark (v0.23.0devmm2 - MM2 dev)](https://github.com/FelixKrueger/Bismark/tree/MM2)_

## What's included

```text
Preprocess/
├── README.md
├── preprocess.sh
└── script/
    ├── bulk_test2_v3p.py
    ├── cut_NDR_support_dep.py
    └── NDR_tabix-5GCH.pl
```

## Usage
Common preset
```
BASEDIR=/path/to/thisfolder
FUNCTION=${BASEDIR}/preprocess.sh
CONFIG=scNanoCool_config
```

Step1: demultiplex
```
s01.demultiplex.sh $work_dir $batch $work_dir/barcode/Barcode.$batch.fa
```

Step2-5: trim, align, filter, extract meth

_Noteworthy that the function do_align requires a lot of resources, at least 14 cores and 42 G memory per cell_
```
## for single cell
source $FUNCTION $work_dir $sc_sample $CONFIG
do_trim
do_align
do_filter_read_mCH
do_methylpy_meth_extractor
```

Step6: meta-statistics for quality control
```
## for single cell
source $FUNCTION $work_dir $sc_sample $CONFIG
do_methylpy_CT_CA_efficiency_v2
do_meth_profile_qc GCH TSS
do_computeMatrix WCG GENE
do_basic_stat
do_read_length
```

## Last modification

2022/10/15

## Last modification

xxhui@stu.pku.edu.cn
