#!/bin/sh
############################################
### nanoCool Pipeline                    ###
### Author: Xiaohui Xue                  ###
### Last Modification: 2022-10-12        ###
### Email: xxh@stu.pku.edu.cn            ###
############################################

dir=$1
sp=$2

echo "workdir=$dir"
echo "sample=$sp"

CONFIG=nanoCool_config

while read line
do
    eval "$line"
done < $CONFIG

function do_cluBCpG-pdm(){
    echo "=================================="
    echo "Find passively demethylated regions begin at: "`date +%Y-%m-%d,%H:%M:%S`
    chr=$1
    sample_paired=$2
    suffix=clean_bismark_mm2.sort.rmdup.bam
    bam_a=$bam_dir/$sp/$sp.$suffix
    bam_b=$bam_dir/$sample_paired/$sample_paired.$suffix

    source $activate cluBCpG
    $clubcpg_pdm \
    --bin_size 500 \
    -a $bam_a -b $bam_b \
    -r $ref -o $bam_dir/$sp \
    -n 20 -chr $chr

    mv $bam_dir/$sp/CompleteBins.$sp.$suffix.$sample_paired.$suffix.$chr.csv \
    $bam_dir/$sp/CompleteBins.$sp.$sample_paired.$chr.csv

    echo "Find passively demethylated regions end at: "`date +%Y-%m-%d,%H:%M:%S` 
}

