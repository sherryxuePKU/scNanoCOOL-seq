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

CONFIG=$3

while read line
do
    eval "$line"
done < $CONFIG

function do_CNV(){
	echo "=================================="
	echo "Calling CNV begin at: "`date +%Y-%m-%d,%H:%M:%S`
	window=$1
	mkdir -p $cnv_dir/$sp/${window}M
	
	bam=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.CHlt80per.bam
	config=$cnv_dir/$sp/${window}M/$sp.${window}M_config.tmp.txt
	echo "
[general]

chrLenFile = $chr_len
ploidy = 3
breakPointThreshold = .6
chrFiles = $freec_db
window = ${window}000000
maxThreads = 2
outputDir = $cnv_dir/$sp/${window}M
BedGraphOutput = TRUE
samtools = $samtools_exe
minExpectedGC = 0.10
maxExpectedGC = 0.30
sex=XY
[sample]

mateFile = $bam
inputFormat = BAM
mateOrientation = 0

[control]

mateFile = /path/to/control.bam
inputFormat = BAM
mateOrientation = 0

	" > $config
	
	$freec_exe -conf $config
	echo "Calling CNV end at: "`date +%Y-%m-%d,%H:%M:%S`
}
