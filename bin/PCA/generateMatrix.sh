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

function do_element_meth(){
    echo "=================================="
    echo "Calc regional WCG or GCH meth begin at: "`date +%Y-%m-%d,%H:%M:%S`
    minC=$1
    Ctype=$2
    element=$3

    case $element in
    "TSS_1k") ref=$TSS_1k_bed;;
    "TSS_0.75k") ref=$TSS_750_bed;;
    "genebody") ref=$genebody_bed;;
    "tile") 
        ref=None
        window=$4
        depth=$5
    ;;
    esac

    if [ "$element" != "tile" ];then
        outdir=$info_dir/PCA/${Ctype}_${element}
        mkdir -p $outdir

        input=$meth_dir/$sp/${sp}.NOMe.genome.${Ctype}_report.txt
        output=$outdir/$sp.${element}_${Ctype}.bed

        $bedtools_exe intersect -wa -wb \
        -a $ref -b $input \
        | awk '{OFS="\t";print $1,$2,$3,$4,$11/($11+$12)}' \
        | $bedtools_exe groupby -i - -g 1,2,3,4 -c 5,5 -o count,mean \
        | awk '{OFS="\t";if($5>="'$minC'")print $1,$2,$3,$6}' > $output
    else
        outdir=$info_dir/PCA/${Ctype}_${element}.${window}bp_${depth}X
        mkdir -p $outdir

        input=$meth_dir/$sp/$sp.NOMe.genome.${Ctype}_report.txt
        output=$outdir/${sp}.${window}bp_${depth}X_${Ctype}.bed

        $perl_exe $tile_pl $input $window $depth $output
    fi
    echo "Calc regional WCG or GCH meth end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_generate_matrix(){
    echo "=================================="
    echo "Generate matrix for reduce dimension begin at: "`date +%Y-%m-%d,%H:%M:%S`
    Ctype=$1
    element=$2
    na_per=$3
    prefix=$4

    indir=$info_dir/PCA/${Ctype}_${element}
    output=$info_dir/PCA/$prefix.${Ctype}_${element}_NA_Per${na_per}.mtx

    cd $indir &&\
    $bedtools_exe unionbedg -i ./* \
    -header -filler NA \
    -names `ls | cut -d "." -f 1 | xargs` \
    | awk -v per=$na_per \
    '{na=0;for(i=4;i<=NF;i++){if($i=="NA")na+=1}
    if(na<=int(per*(NF-3)))print $0}' > $output
    
    echo "Generate matrix for reduce dimension end at: "`date +%Y-%m-%d,%H:%M:%S`      
}

