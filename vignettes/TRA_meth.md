## Table of contents

- [Quick start](#quick-start)
- [What's included](#whats-included)
- [Usage](#Usage)
- [Last modification](#last-modification)
- [Contact](#contact)

## Quick start

the folder /bin/TRA/ is used to evaluate the methylation level of WCG/GCH site flanking the breakpoint of BCR-ABL1 translocation with following steps

- Step1: align the reads to assembly hg38 and contig, respectively by single cell
- Step2: extract the wildtype and translocated reads aligned across the breakpoint 
- Step3: merge the wildtype and translocated reads by pseudo-bulk 
- Step4: calc the methylation level of WCG/GCH site   


## What's included

```text
TRA/
├── README.md
├── TRA_meth.sh
└── reference/
    ├── contig.BCR_ABL1.fusion_point.bed
    └── hg38.BCR_ABL1.fusion_point.bed
```

## Usage
Common preset

    BASEDIR=/path/to/thisfolder
    FUNCTION=${BASEDIR}/TRA_meth.sh
    CONFIG=scNanoCool_config
    
    region_bcr_contig=reference/contig.BCR_ABL1.fusion_point.bed
    region_bcr_hg38=reference/hg38.BCR_ABL1.fusion_point.bed

Step1: align to hg38 or contig

Align the read to hg38 as shown in the folder of preprocess

    ## for single cells to contig
    source $FUNCTION $work_dir $sc_sample $CONFIG
    do_align_to_contig

Step2: extract targeted reads

    ## for single cells
    source $FUNCTION $work_dir $psbulk_sample $CONFIG
    do_extract_fusion_reads.contig $region_bcr_contig BCR_ABL1
    do_extract_fusion_reads.hg38 $region_bcr_hg38 BCR_ABL1

Step3: merge the reads

    ## for pseudo-bulk
    source $FUNCTION $work_dir $psbulk_sample $CONFIG
    do_merge_bam.contig $region_bcr_contig BCR_ABL1
    do_merge_bam.hg38 $region_bcr_hg38 BCR_ABL1

Step4: calc the methylation level

    ## for pseudo-bulk
    source $FUNCTION $work_dir $psbulk_sample $CONFIG
    do_methylpy_TRA_hg38.v2 BCR_ABL1
    do_methylpy_TRA_contig BCR_ABL1

## Last modification

2022/10/15

## Last modification

xxhui@stu.pku.edu.cn
