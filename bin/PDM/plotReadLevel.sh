#!/bin/sh
############################################
### nanoCool Pipeline                    ###
### Author: Xiaohui Xue                  ###
### Last Modification: 2022-10-16        ###
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

## plot the regions without SNVs
function do_plot_tanghulu.region(){
    echo "=================================="
    echo "draw tanghulu plot of given regions begin at: "`date +%Y-%m-%d,%H:%M:%S`
    # region=$1
    region=$1
    name=$2
    ctype=WCG
    
    in_bam=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.bam
    outdir=$info_dir/Tanghulu/$sp/${name}_noSNV

    mkdir -p $outdir

    source $activate cluBCpG
    cat $region | while read line
    do
        rg_name=`echo $line | awk '{print $4}'`
        full_pos=`echo $line | awk '{print $1":"$2"-"$3}'`
        tmp_bam=$bam_dir/$sp/$sp.$rg_name.bam
        out_mtx=$outdir/$sp.$rg_name.${ctype}_addCB.txt
        out_pdf=$outdir/$sp.$rg_name.${ctype}_addCB.pdf

        $samtools_exe view -bSh -q 20 $in_bam $full_pos > $tmp_bam &&\
        $samtools_exe index $tmp_bam &&\
        $python $tanghulu_py -b $tmp_bam -c $full_pos -o $out_mtx -t $ctype
        $rscript_exe $tanghulu_r $out_mtx $out_pdf $rg_name":"$title
        rm -f $tmp_bam $tmp_bam.bai $out_mtx
    done
   
    $pdfunite $outdir/$sp.*.pdf $outdir/../$sp.${name}_noSNV.merged.pdf

    echo "draw tanghulu plot of given regions end at: "`date +%Y-%m-%d,%H:%M:%S`        
}

## plot the regions with SNVs
function do_plot_tanghulu.SNV(){
    echo "=================================="
    echo "draw tanghulu plot of given regions begin at: "`date +%Y-%m-%d,%H:%M:%S`
    # region=$1 ## format for region.bed: chr start end name
    # region=$ICR
    region=$1
    name=$2
    ctype=WCG

    in_bam=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.bam
    snv_bed=$bam_dir/$sp/$sp.longshot_filter.${name}_SNV.bed

    $samtools_exe view -bSh $in_bam -L $region -o $in_bam.$name.bam
    $samtools_exe index $in_bam.$name.bam

    $longshot_exe -F -n --min_cov 2 --bam $in_bam.$name.bam \
    --ref $ref --out $bam_dir/$sp/$sp.longshot_raw.vcf

    grep -v "#" $bam_dir/$sp/$sp.longshot_raw.vcf \
     | awk -v OFS="\t" '{
        split($8,a,";");sub(/DP=/,"",a[1]);sub(/AC=/,"",a[2]);split(a[2],b,",");
        if(!($4=="C" && $5=="T") && !($4=="G" && $5=="A") && b[1]>=2 && b[2]>=2)
        print $1,$2-1,$2,$4,$5,a[1],a[2]}' \
    > $bam_dir/$sp/$sp.longshot_filter.SNV.bed

    $bedtools_exe intersect -a $region \
    -b $bam_dir/$sp/$sp.longshot_filter.SNV.bed -wa -wb \
     | $bedtools_exe groupby -i - -g 1,2,3,4 -c 4 -o count \
     | awk '{if($5==1)print}' > $bam_dir/$sp/$sp.longshot_filter.single_SNV.bed

    $bedtools_exe intersect -wa -wb \
    -a $bam_dir/$sp/$sp.longshot_filter.single_SNV.bed \
    -b $bam_dir/$sp/$sp.longshot_filter.SNV.bed > $snv_bed

    # rm -f  $bam_dir/$sp/$sp.longshot_filter.single_SNV.bed

    outdir=$info_dir/Tanghulu/$sp/${name}_SNV

    mkdir -p $outdir

    source $activate cluBCpG
    cat $snv_bed | while read line
    do
        rg_name=`echo $line | awk '{print $4}'`
        plot_pos=`echo $line | awk '{print $1":"$2"-"$3}'`
        rg_pos=`echo $line | awk '{print $1":"$2"-"$3}'`
        snv_pos=`echo $line | awk '{print $7}'`
        tmp_bam=$bam_dir/$sp/$sp.$rg_name.bam
        out_mtx=$outdir/$sp.$rg_name.${ctype}_addCB.txt
        out_pdf=$outdir/$sp.$rg_name.${ctype}_addCB.pdf

        $samtools_exe view -bSh $in_bam $rg_pos  > $tmp_bam &&\
        $samtools_exe index $tmp_bam &&\
        $python $tanghulu_py -b $tmp_bam -c $plot_pos -o $out_mtx -t $ctype -s $snv_pos
        $rscript_exe $tanghulu_snv_r $out_mtx $out_pdf $rg_name":"$rg_pos $snv_pos

        rm -f $tmp_bam $tmp_bam.bai
    done

    for gene in `cat $snv_bed | awk '{split($4,a,"_");print a[1]}' | sort -u`
    do
       $pdfunite $outdir/$sp.${gene}_*.pdf $outdir/$sp.$gene.pdf
       rm -f $outdir/$sp.${gene}_*.pdf
    done

    rm -f $bam_dir/$sp/$sp.longshot_filter.single_SNV.bed \
    $in_bam.$name.bam $in_bam.$name.bam.bai

    echo "draw tanghulu plot of given regions end at: "`date +%Y-%m-%d,%H:%M:%S`        
}
