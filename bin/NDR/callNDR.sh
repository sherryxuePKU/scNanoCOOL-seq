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

## Step1:
function do_merge_bam(){
    echo "=================================="
    echo "Merge bam files begin at: "`date +%Y-%m-%d,%H:%M:%S`
    list=$dir/meta/$sp.sample.list
    bam_list=$dir/meta/$sp.bam.list
    outdir=$bam_dir/$sp
    bam_suffix=clean_bismark_mm2.sort.rmdup.CHlt80.bam

    mkdir -p $outdir
    awk '{print "'$bam_dir'""/"$1"/"$1".""'$bam_suffix'"}' $list > $bam_list
    
    $bamtools_exe merge -list $bam_list \
    -out $outdir/$sp.$bam_suffix &&\
    $samtools_exe sort -o $outdir/$sp.$bam_suffix \
    $outdir/$sp.$bam_suffix &&\
    $samtools_exe index $outdir/$sp.$bam_suffix
    echo "Merge bam files end at: "`date +%Y-%m-%d,%H:%M:%S`    
}

## Step2
function do_call_NDR_v2(){
    echo "=================================="
    echo "Call NDR by MethGC begin at: "`date +%Y-%m-%d,%H:%M:%S`
    chr=$1
    depth=3
    input=$meth_dir/$sp/${sp}.NOMe.genome.GCH_report.txt
    chr_GCH=$meth_dir/$sp/$chr.${sp}.GCH_report.txt.gz
    chr_CT=$meth_dir/$sp/$chr.${sp}.CT_total.xls
    output=$NDR_dir/$sp/$sp.NOMe.$chr.NDR_MethGC.${depth}Depth_5GCH.tsv

    mkdir -p $NDR_dir/$sp
    awk '{OFS="\t";if($1=="'$chr'")print $1,$3,$7,$6}' $input |\
    $bgzip_exe -c > $chr_GCH && $tabix_exe -s 1 -b 2 -e 2 $chr_GCH 
    awk '{if($1=="'$chr'"){un+=$7;meth+=$6}}END{print un"\t"meth}' $input > $chr_CT

    $python2_exe $findNDR $chr_GCH $chr_CT $depth > $output.tmp &&\
    $perl_exe $NDR_5GCH $output.tmp $chr_GCH > $output

    rm $chr_GCH $chr_CT $output.tmp $chr_GCH.tbi

    echo "Call NDR by MethGC end at: "`date +%Y-%m-%d,%H:%M:%S`    
}

## Step2: merge NDRs called by chr
function do_merge_chrNDR(){
    echo "=================================="
    echo "Call NDR by MethGC begin at: "`date +%Y-%m-%d,%H:%M:%S`
    depth=3
    merged_NDR=$NDR_dir/$sp/$sp.NOMe.genome.NDR_MethGC.${depth}Depth_5GCH.tsv

    cat $NDR_dir/$sp/$sp.NOMe.chr*.NDR_MethGC.${depth}Depth_5GCH.tsv |\
    sort -k1,1 -k2,2n > $merged_NDR &&\
    rm $NDR_dir/$sp/$sp.NOMe.chr*.NDR_MethGC.${depth}Depth_5GCH.tsv

    awk -v OFS="\t" '{len=$3-$2;if(len>=240)print $1,$2,$3}' $merged_NDR > $merged_NDR.len_ge240bp.bed

    $bedtools_exe intersect -wa -a $merged_NDR.len_ge240bp.bed -b $TSS_2k_bed > $merged_NDR.proximal_NDR.bed
    $bedtools_exe intersect -wa -v -a $merged_NDR.len_ge240bp.bed -b $TSS_2k_bed > $merged_NDR.distal_NDR.bed 
    $bedtools_exe intersect -wa -a $merged_NDR.len_ge240bp.bed -b $promoter_bed > $merged_NDR.promoter_NDR.bed
    echo "Call NDR by MethGC end at: "`date +%Y-%m-%d,%H:%M:%S`
}

## Step3: calc GCH_meth in NDR (optional)
function do_computeMatrix_NDR(){
	echo "=================================="
	echo "ComputeMatrix begin at: "`date +%Y-%m-%d,%H:%M:%S`
	region=$1
    C_type=GCH
	input=$meth_dir/$sp/${sp}.NOMe.genome.${C_type}_report.txt
    outdir=$info_dir/meth_profile_data/${C_type}_NDR
	abs_mtx=$outdir/${sp}.${C_type}_NDR.absoluteDistance.txt
    abs_df=$outdir/${sp}.${C_type}_NDR.absoluteDistance_merge.txt
	
    mkdir -p $outdir
	
	$perl_exe $AbD_NDR_pl $region $input > $abs_mtx &&\
    $rscript_exe $merge_mtx_r $abs_mtx $abs_df 

    rm -f $abs_mtx
	echo "ComputeMatrix end at: "`date +%Y-%m-%d,%H:%M:%S`
}

