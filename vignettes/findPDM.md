## Table of contents

- [Quick start](#quick-start)
- [What's included](#whats-included)
- [Usage](#Usage)
- [Last modification](#last-modification)
- [Contact](#contact)

## Quick start

the folder /bin/PDM/ is used to plot read-level methylation of passively demethylated regions (PDMs) in the following steps. It is worth mentioning that the main script is modified from the package [CluBCpG](https://github.com/waterlandlab/CluBCpG) and adapted to NOMe-seq.

- Step1: merge the .bam files of single cells to pseudo-bulk
- Step2: find the region with high read-level WCG methylation level in the control sample
- Step3: generate the methylation matrix of WCG sites and plot the matrix
- Step4: de novo find SNV in the given regions and plot 

## What's included

```text
PDM/
├── README.md
├── findPDM.sh
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

Step2: find the PDM regions by chromosome

    source $FUNCTION $work_dir $psbulk_sample $CONFIG
    for chr in {1..22} X Y
    do
        do_cluBCpG-pdm chr$chr nanoCool_5aza.1uM_345.CHlt80per
    done

Step3: generate the WCG matrix and plot the result

1) merge the chr-level result 
2) filter the merged result
    ```
    rm -f ../../meta/nanoCool_5aza.1002.pdm.*.bed
    cat CompleteBins.nanoCool_5aza.DMSO_all.CHlt80per.nanoCool_5aza.1uM_345.CHlt80per.chr*.csv | grep -v "Not" \
    | awk -v OFS="\t" min_cpg=5 min_read=3 '{split($1,a,"_");split($2,b,",");if($3>=min_read && $4>=min_cpg && $5>=min_read && $6>=min_cpg){if(b[1]==0 && b[2]==0)print a[1],a[2]-500,a[2],"fullzero_"NR,$2>>"../../meta/nanoCool_5aza.1002.pdm.fullzero.bed";if(b[2]==0 && b[3]==0)print a[1],a[2]-500,a[2],"fullmeth_"NR,$2>>"../../meta/nanoCool_5aza.1002.pdm.fullmeth.bed";if(b[1]==0 && b[3]==0)print a[1],a[2]-500,a[2],"mosaic_"NR,$2>>"../../meta/nanoCool_5aza.1002.pdm.mosaic.bed"}}'
    ```
3) plot the regions
    ```
    region_fullzero=meta/nanoCool_5aza.1002.pdm.fullzero.bed
    region_fullmeth=meta/nanoCool_5aza.1002.pdm.fullmeth.bed
    region_mosaic=meta/nanoCool_5aza.1002.pdm.mosaic.bed

    source $FUNCTION $work_dir $psbulk_sample $CONFIG

    do_plot_tanghulu.region $region_fullzero PDM_fullzero
    do_plot_tanghulu.region $region_fullmeth PDM_fullmeth
    do_plot_tanghulu.region $region_mosaic PDM_mosaic

    do_plot_tanghulu.SNV $region_fullmeth PDM_fullmeth
    do_plot_tanghulu.SNV $region_fullzero PDM_fullzero
    ```

## Last modification

2022/10/16

## Contact

xxhui@stu.pku.edu.cn
