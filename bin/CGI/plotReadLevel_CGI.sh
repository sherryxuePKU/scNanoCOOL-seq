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


function do_plot_tanghulu.CGI(){
    echo "=================================="
    echo "draw tanghulu plot of given regions begin at: "`date +%Y-%m-%d,%H:%M:%S`
    # region=$1
    # region=$ICR
    ctype=WCG
    
    in_bam=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.bam
    outdir=$info_dir/Tanghulu/$sp/CGI_noSNV

    mkdir -p $outdir
    # mkdir -p $info_dir/Tanghulu/$sp/plot
    # mkdir -p $info_dir/Tanghulu/$sp/matrix

    region=$CGI_ref.len_filtered.bed
    awk '{len=$3-$2;if(len>=300 && len <=500)print}' $CGI_ref > $region

    source $activate cluBCpG
    cat $region | while read line
    do
        rg_name=`echo $line | awk '{print $4}'`
        full_pos=`echo $line | awk '{print $1":"$2-500"-"$3+500}'`
        title=`echo $line | awk '{print $1":"$2"-"$3}'`
        tmp_bam=$bam_dir/$sp/$sp.$rg_name.bam
        out_mtx=$outdir/$sp.$rg_name.${ctype}_addCB.txt
        out_pdf=$outdir/$sp.$rg_name.${ctype}_addCB.pdf

        $samtools_exe view -bSh $in_bam $title  > $tmp_bam &&\
        $samtools_exe index $tmp_bam &&\
        $python $tanghulu_py -b $tmp_bam -c $full_pos -o $out_mtx -t $ctype
        $rscript_exe $tanghulu_CGI_r $out_mtx $out_pdf $rg_name":"$title
	rm -f $tmp_bam $tmp_bam.bai
    done
   
    $pdfunite $outdir/$sp.CGI*.pdf $outdir/$sp.merged.CGI.pdf 
    rm -f $region
    echo "draw tanghulu plot of given regions end at: "`date +%Y-%m-%d,%H:%M:%S`        
}
