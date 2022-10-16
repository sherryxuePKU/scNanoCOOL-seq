## Table of contents

- [Quick start](#quick-start)
- [What's included](#whats-included)
- [Usage](#Usage)
- [Last modification](#last-modification)
- [Contact](#contact)

## Quick start

the folder /bin/NDR/ is used to call nucleosome-depleted regions (NDRs) by home-made script with following steps

- Step1: merge .bam files of single cells to pseudo-bulk .bam file and calc the GCH methylation level
- Step2: calc NDRs with the home-made script by chromosome 
- Step3: filter the NDRs by the number of GCH sites located within the NDRs and covered reads 
- Step4: merge the results of chromosomes' for pseudo-bulk
- Step5: validate the NDRs with benchmarking data from ENCODE
- Step6: calc the methylation level of WCG/GCH and the signals of other omics (e.g. ChIP-seq) around the NDRs   


## What's included

```text
NDR/
├── README.md
├── callNDR.sh
└── script/
    ├── bulk_test2_v3p.py
    ├── cut_NDR_support_dep.py
    ├── NDR_tabix-5GCH.pl
    ├── AbsoluteDistance_NOMe_NDR.pl
    ├── enrichheatmap_NDR.R
    ├── NDR_benchmark.R
    └── plot_meth_signal.R
```

## Usage
Common preset
```
BASEDIR=/path/to/thisfolder
FUNCTION_NDR=${BASEDIR}/callNDR.sh
FUNCTION_preprocess=${BASEDIR}/preprocess.sh
CONFIG=scNanoCool_config
```

Step1: merge .bam files and calc GCH
```
## for pseudo-bulk
source $FUNCTION_preprocess $work_dir $psbulk_sample $CONFIG
source $FUNCTION_NDR $work_dir $psbulk_sample $CONFIG
do_merge_bam
do_methylpy_meth_extractor
```

Step2 & Step3: calc NDRs & filter them
```
## for chromosome
source $FUNCTION_NDR $work_dir $psbulk_sample $CONFIG
for chr in {1..22} X Y
do
    do_call_NDR_v2 chr${chr}
done
```

Step4: merge the result and classify the NDRs into proximal NDRs and distal NDRs
```
## for pseudo-bulk
source $FUNCTION_NDR $work_dir $psbulk_sample $CONFIG
do_merge_chrNDR
```

Step5: benchmarking

1) calc the precision-recall matrix by [intervene](https://github.com/asntech/intervene)
    ```
    intervene pairwise -i GRCh38_cCRE_K562.bed.gz ENCFF057UYP.bed.gz ENCFF185XRG.bed.gz $scNanoCool_NDR $scCOOL_NDR
    ``` 
2) plot the result of intervene

    run the NDR_benchmark.R in the Rstudio

Step6: calc the signal of other omics

1) calc the methylation level of WCG/GCH
    ```
    ## GCH
    perl AbsoluteDistance_NOMe_NDR.pl $proximal_NDR nanoCool_5aza.0.NOMe.genome.CHlt80.GCH_report.txt > nanoCool_5aza.0.CHlt80per.3Depth_5GCH.proximal_NDR.GCH_methlevel.tsv

    ## WCG
    perl AbsoluteDistance_NOMe_NDR.pl $proximal_NDR nanoCool_5aza.0.NOMe.genome.CHlt80.WCG_report.txt > nanoCool_5aza.0.CHlt80per.3Depth_5GCH.proximal_NDR.WCG_methlevel.tsv
    ```
2) calc the signal of ChIP-seq
    ```
    ## CTCF
    Rscript enrichheatmap_NDR.R ENCFF041ODC.bed.gz \
    $proximal_NDR nanoCool_5aza.0.CHlt80per.3Depth_5GCH.proximal_NDR.CTCF_signal.rds

    ## POLR2A
    Rscript enrichheatmap_NDR.R ENCFF355MNE.bed.gz  \
    nanoCool_5aza.0.CHlt80per.3Depth_5GCH.proximal_NDR.POLR2A_signal.rds
    ```
3) plot the methylation levels and signals around proximal NDRs
    run the plot_meth_signal.R in the Rstudio

### Download links of reference data
| Antibody | Link of ENCODE |
| ------ | ----------- |
| CTCF  | ENCFF041ODC (https://www.encodeproject.org/files/ENCFF041ODC/) |
| POLR2A | ENCFF355MNE (https://www.encodeproject.org/files/ENCFF355MNE/) |

## Last modification

2022/10/15

## Last modification

xxhui@stu.pku.edu.cn
