#!/bin/sh
dir=$1
sp=$2
#barcode=$dir/barcode/Barcode.M5.fa
barcode=$3
nanoplexer_exe=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/xuexiaohui/software/miniconda3/bin/nanoplexer

mkdir -p $dir/00.raw_data/$sp

$nanoplexer_exe -t 8 \
-b $barcode \
-p $dir/00.raw_data/$sp \
$dir/00.raw_data/$sp.fastq.gz

#inner=`echo $sp | awk '{split($1,a,"_");print a[2]}'`
#rename .fastq ${inner}.fastq $dir/00.raw_data/$sp/*
