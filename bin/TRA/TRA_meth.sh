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

## bismark_mm2 branch align to K562_TRA_contig.fa
function do_align_to_contig(){
    echo "=================================="
    echo "Align by bismark_minimap2 to contig begin at: "`date +%Y-%m-%d,%H:%M:%S`
    mkdir -p $bam_dir/$sp.contig
    trim=$trim_dir/$sp.clean.fastq
    bam=$bam_dir/$sp.contig/$sp.clean_bismark_mm2
    $bismark_mm2 --mm2 --parallel 1 --non_directional \
    --fastq --mm2_maximum_length 100000 \
    --path_to_minimap2 $minimap2_dir \
    --samtools_path $bin_dir/samtools-1.9/bin \
    --output_dir $bam_dir/$sp.contig \
    --temp_dir $bam_dir/$sp.contig \
    $bismark_mm2_contig_db $trim &&\
    $samtools_exe view -uSb -@ 3\
    -t $bismark_mm2_contig_db/K562_translocation_contig.fa $bam.bam |\
    $samtools_exe sort -m 900M -@ 3 -T $sp -o $bam.sort.bam &&\
    $samtools_exe rmdup -s $bam.sort.bam $bam.sort.rmdup.bam &&\
    $samtools_exe index $bam.sort.rmdup.bam &&\
    rm -f $bam.sort.bam
    echo "Align by bismark_minimap2 to contig end at: "`date +%Y-%m-%d,%H:%M:%S`    
}

function do_extract_fusion_reads.contig(){
    echo "=================================="
    echo "add tag to reads of some regions begin at: "`date +%Y-%m-%d,%H:%M:%S`
    region=$1
    name=$2
    bam_prefix=$bam_dir/$sp.contig/$sp.clean_bismark_mm2.sort.rmdup

    # $methylpy_exe bam-quality-filter \
    # --input-file $bam_prefix.bam \
    # --output-file $bam_prefix.CHlt80per \
    # --ref-fasta $bismark_mm2_contig_db/K562_translocation_contig.fa \
    # --min-num-ch 1 --max-mch-level 0.8 &&\
    # $samtools_exe index $bam_prefix.CHlt80per.bam

    $samtools_exe view -bSh -L $region $bam_prefix.CHlt80per.bam > $bam_prefix.$name.bam &&\
    $samtools_exe view $bam_prefix.$name.bam -h |\
    awk '{if($1~/@/)print;else print $0"\t""CB:Z:""'$sp'"}' |\
    $samtools_exe view -bSh - > $bam_prefix.${name}_Mut_CBtag.bam &&\
    $samtools_exe index $bam_prefix.${name}_Mut_CBtag.bam

    echo "add tag to reads of some regions end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_merge_bam.contig(){
    echo "=================================="
    echo "Merge bam files begin at: "`date +%Y-%m-%d,%H:%M:%S`
    region=$1
    name=$2
    # upstream=$2
    # downstream=$3
    
    list=$dir/meta/$sp.TRA.sample.list
    bam_list=$dir/meta/$sp.bam.list
    outdir=$bam_dir/$sp
    bam_suffix=clean_bismark_mm2.sort.rmdup.${name}_Mut_CBtag

    ## merge bam files
    mkdir -p $outdir.contig
    awk '{print "'$bam_dir'""/"$1".contig/"$1".""'$bam_suffix'"".bam"}' $list > $bam_list
    
    $bamtools_exe merge -list $bam_list -out $outdir.contig/$sp.$bam_suffix.bam &&\
    $samtools_exe sort -o $outdir.contig/$sp.$bam_suffix.sort.bam $outdir.contig/$sp.$bam_suffix.bam

    ## intersect with upstream region
    upstream=`awk '{print $1":"$2-150"-"$2-140}' $region`

    $samtools_exe index $outdir.contig/$sp.$bam_suffix.sort.bam &&\
    $samtools_exe view -bSh $outdir.contig/$sp.$bam_suffix.sort.bam \
    $upstream -o $outdir.contig/$sp.$bam_suffix.upstream.bam
    $samtools_exe index $outdir.contig/$sp.$bam_suffix.upstream.bam

    ## intersect with downstream region
    downstream=`awk '{print $1":"$3+150"-"$2+160}' $region`

    $samtools_exe view -bSh $outdir.contig/$sp.$bam_suffix.upstream.bam \
    $downstream -o $outdir.contig/$sp.$bam_suffix.fusion.bam &&\
    $samtools_exe index $outdir.contig/$sp.$bam_suffix.fusion.bam

    rm -f $outdir.contig/$sp.$bam_suffix.bam \
    $outdir.contig/$sp.$bam_suffix.upstream.bam \
    $outdir.contig/$sp.$bam_suffix.upstream.bam.bai \
    $outdir.contig/$sp.$bam_suffix.sort.bam \
    $outdir.contig/$sp.$bam_suffix.sort.bam.bai 
    echo "Merge bam files end at: "`date +%Y-%m-%d,%H:%M:%S`    
}

function do_methylpy_TRA_contig(){
    echo "=================================="
    echo "Calc methylation level of translocated reads begin at: "`date +%Y-%m-%d,%H:%M:%S`
    # region=$1
    name=$1

    bam_mut_suffix=clean_bismark_mm2.sort.rmdup.${name}_Mut_CBtag
    input=$bam_dir/$sp.contig/$sp.$bam_mut_suffix.fusion
    output_prefix=$meth_dir/$sp.contig/$sp.${name}.NOMe.contig

    # $methylpy_exe bam-quality-filter \
    # --input-file $input.bam \
    # --output-file $input.CHlt80per.bam \
    # --ref-fasta $ref \
    # --min-num-ch 1 \
    # --max-mch-level 0.8 &&\
    # $samtools_exe index $input.CHlt80per.bam
    mkdir -p $meth_dir/$sp.contig

    $samtools_exe sort -o $input.sort.bam $input.bam

    $methylpy_exe call-methylation-state \
    --input-file $input.sort.bam --sample $sp.${name} \
    --ref-fasta $bismark_mm2_contig_db/K562_translocation_contig.fa \
    --paired-end False --num-upstream-bases 1 --num-downstream-bases 1 \
    --num-procs 4 --path-to-output $meth_dir/$sp.contig

    zcat $meth_dir/$sp.contig/allc_${sp}.${name}.tsv.gz \
    | awk -v OFS="\t" '{if($4 ~ /[AT]CG/ && $1!="lambda")print "chr"$1,$2-1,$2,$3,$4,$5,$6-$5,"WCG",$4}' \
    > $output_prefix.WCG_report.txt
    
    zcat $meth_dir/$sp.contig/allc_${sp}.${name}.tsv.gz \
    | awk -v OFS="\t" '{if($4 ~ /GC[ATC]/ && $1!="lambda")print "chr"$1,$2-1,$2,$3,$4,$5,$6-$5,"GCH",$4}' \
    > $output_prefix.GCH_report.txt
    
    awk -v OFS="\t" '{if($1!="chrLambda")print $1,$2,$3,$6/($6+$7)}' $output_prefix.WCG_report.txt \
    | sort -k1,1 -k2,2n > $output_prefix.WCG_meth.tmp.bedgraph &&\
    $bedGraphToBigWig $output_prefix.WCG_meth.tmp.bedgraph \
    $chrom_size_contig $output_prefix.WCG_meth.bw

    awk -v OFS="\t" '{if($1!="chrLambda")print $1,$2,$3,$6/($6+$7)}' $output_prefix.GCH_report.txt \
    | sort -k1,1 -k2,2n > $output_prefix.GCH_meth.tmp.bedgraph &&\
    $bedGraphToBigWig $output_prefix.GCH_meth.tmp.bedgraph \
    $chrom_size_contig $output_prefix.genome.GCH_meth.bw

    rm -f $output_prefix.WCG_meth.tmp.bedgraph \
    $output_prefix.GCH_meth.tmp.bedgraph \
    $meth_dir/$sp.contig/allc_${sp}.${name}.tsv.gz \
    $meth_dir/$sp.contig/allc_${sp}.${name}.tsv.gz.idx

    echo "Calc methylation level of translocated reads end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_extract_fusion_reads.hg38(){
    echo "=================================="
    echo "add tag to reads of some regions begin at: "`date +%Y-%m-%d,%H:%M:%S`
    region=$1
    name=$2
    bam_prefix=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup

    $samtools_exe view -bSh -L $region $bam_prefix.CHlt80per.bam > $bam_prefix.$name.bam &&\
    $samtools_exe view $bam_prefix.$name.bam -h |\
    awk '{if($1~/@/)print;else print $0"\t""CB:Z:""'$sp'"}' |\
    $samtools_exe view -bSh - > $bam_prefix.${name}_WT_CBtag.bam &&\
    $samtools_exe index $bam_prefix.${name}_WT_CBtag.bam

    echo "add tag to reads of some regions end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_merge_bam.hg38(){
    echo "=================================="
    echo "Merge bam files begin at: "`date +%Y-%m-%d,%H:%M:%S`
    region=$1
    name=$2
    
    list=$dir/meta/$sp.TRA.sample.list
    bam_list=$dir/meta/$sp.bam.list
    outdir=$bam_dir/$sp.CHlt80per
    bam_wt_suffix=clean_bismark_mm2.sort.rmdup.${name}_WT_CBtag
    bam_mut_suffix=clean_bismark_mm2.sort.rmdup.${name}_Mut_CBtag.fusion

    ## merge bam files
    mkdir -p $outdir
    awk '{print "'$bam_dir'""/"$1"/"$1".""'$bam_wt_suffix'"".bam"}' $list > $bam_list
    
    $bamtools_exe merge -list $bam_list -out $outdir/$sp.$bam_wt_suffix.bam &&\
    $samtools_exe sort -o $outdir/$sp.$bam_wt_suffix.sort.bam $outdir/$sp.$bam_wt_suffix.bam

    ## filter fusion read
    cat <($samtools_exe view $bam_dir/$sp.contig/$sp.$bam_mut_suffix.bam |cut -f 1 | uniq) \
    <($samtools_exe view $outdir/$sp.$bam_wt_suffix.sort.bam |cut -f 1 | uniq) |\
    sort | uniq -u > $outdir/$sp.${name}_WT.reads.list

    $python3_exe $extract_bam \
    -b $outdir/$sp.$bam_wt_suffix.sort.bam \
    -n $outdir/$sp.${name}_WT.reads.list \
    -o $outdir/$sp.${name}_WT_CBtag.filtered.bam &&\
    $samtools_exe sort -o $outdir/$sp.${name}_WT_CBtag.filtered.sort.bam \
    $outdir/$sp.${name}_WT_CBtag.filtered.bam

    mv $outdir/$sp.${name}_WT_CBtag.filtered.sort.bam $outdir/$sp.$bam_wt_suffix.bam
    $samtools_exe index $outdir/$sp.$bam_wt_suffix.bam 

    rm -f $outdir/$sp.${name}_WT_CBtag.filtered.bam \
    $outdir/$sp.${name}_WT.reads.list \
    $outdir/$sp.$bam_wt_suffix.sort.bam \

    ## intersect with upstream and downstream regions
    ## intersect with upstream region
    upstream=`awk '{print $1":"$2-150"-"$2-140}' $region`

    $samtools_exe view -bSh $outdir/$sp.$bam_wt_suffix.bam \
    $upstream -o $outdir/$sp.$bam_wt_suffix.upstream.bam
    $samtools_exe index $outdir/$sp.$bam_wt_suffix.upstream.bam

    ## intersect with downstream region
    downstream=`awk '{print $1":"$3+150"-"$2+160}' $region`

    $samtools_exe view -bSh $outdir/$sp.$bam_wt_suffix.upstream.bam \
    $downstream -o $outdir/$sp.$bam_wt_suffix.fusion.bam &&\
    $samtools_exe index $outdir/$sp.$bam_wt_suffix.fusion.bam

    rm -f $outdir/$sp.$bam_wt_suffix.bam \
    $outdir/$sp.$bam_wt_suffix.bam.bai \
    $outdir/$sp.$bam_wt_suffix.upstream.bam \
    $outdir/$sp.$bam_wt_suffix.upstream.bam.bai

    echo "Merge bam files end at: "`date +%Y-%m-%d,%H:%M:%S`    
}

function do_methylpy_TRA_hg38.v2(){
    echo "=================================="
    echo "Calc methylation level of translocated reads begin at: "`date +%Y-%m-%d,%H:%M:%S`
    # region=$1
    name=$1

    bam_wt_suffix=clean_bismark_mm2.sort.rmdup.${name}_WT_CBtag
    input=$bam_dir/$sp.CHlt80per/$sp.$bam_wt_suffix.fusion
    output_prefix=$meth_dir/$sp.CHlt80per/$sp.${name}.NOMe.genome

    # $methylpy_exe bam-quality-filter \
    # --input-file $input.bam \
    # --output-file $input.CHlt80per.bam \
    # --ref-fasta $ref \
    # --min-num-ch 1 \
    # --max-mch-level 0.8 &&\
    # $samtools_exe index $input.CHlt80per.bam
    $samtools_exe sort -o $input.sort.bam $input.bam

    $methylpy_exe call-methylation-state \
    --input-file $input.sort.bam --sample $sp.${name} \
    --ref-fasta $ref --paired-end False \
    --num-upstream-bases 1 --num-downstream-bases 1 \
    --num-procs 4 --path-to-output $meth_dir/$sp.CHlt80per

    zcat $meth_dir/$sp.CHlt80per/allc_${sp}.${name}.tsv.gz \
    | awk -v OFS="\t" '{if($4 ~ /[AT]CG/ && $1!="lambda")print "chr"$1,$2-1,$2,$3,$4,$5,$6-$5,"WCG",$4}' \
    > $output_prefix.WCG_report.txt
    
    zcat $meth_dir/$sp.CHlt80per/allc_${sp}.${name}.tsv.gz \
    | awk -v OFS="\t" '{if($4 ~ /GC[ATC]/ && $1!="lambda")print "chr"$1,$2-1,$2,$3,$4,$5,$6-$5,"GCH",$4}' \
    > $output_prefix.GCH_report.txt
    
    awk -v OFS="\t" '{if($1!="chrLambda")print $1,$2,$3,$6/($6+$7)}' $output_prefix.WCG_report.txt \
    | sort -k1,1 -k2,2n > $output_prefix.WCG_meth.tmp.bedgraph &&\
    $bedGraphToBigWig $output_prefix.WCG_meth.tmp.bedgraph \
    $chrom_size $output_prefix.WCG_meth.bw

    awk -v OFS="\t" '{if($1!="chrLambda")print $1,$2,$3,$6/($6+$7)}' $output_prefix.GCH_report.txt \
    | sort -k1,1 -k2,2n > $output_prefix.GCH_meth.tmp.bedgraph &&\
    $bedGraphToBigWig $output_prefix.GCH_meth.tmp.bedgraph \
    $chrom_size $output_prefix.GCH_meth.bw

    rm -f $output_prefix.WCG_meth.tmp.bedgraph \
    $output_prefix.GCH_meth.tmp.bedgraph \
    $meth_dir/$sp.CHlt80per/allc_${sp}.${name}.tsv.gz \
    $meth_dir/$sp.CHlt80per/allc_${sp}.${name}.tsv.gz.idx

    echo "Calc methylation level of translocated reads end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_merge_bam.hg38.v2(){
    echo "=================================="
    echo "Merge bam files begin at: "`date +%Y-%m-%d,%H:%M:%S`
    # region=$1
    list=$dir/meta/$sp.sample.list
    bam_list=$dir/meta/$sp.bam.list
    bam_prefix=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup
    bam_chlt80_prefix=$bam_dir/$sp.CHlt80per/$sp.CHlt80per.clean_bismark_mm2.sort.rmdup
    bam_wt_prefix=$bam_dir/$sp/$sp.BCR_ABL1_all_CBtag

    $samtools_exe view -bSh -L $region $bam_chlt80_prefix.bam > $bam_wt_prefix.bam &&\
    $samtools_exe sort -o $bam_wt_prefix.sort.bam $bam_wt_prefix.bam &&\
    $samtools_exe index $bam_wt_prefix.sort.bam

    cat <($samtools_exe view $bam_prefix.BCR_ABL1_Mut_CBtag.fusion.bam |cut -f 1 | uniq) \
    <($samtools_exe view $bam_dir/$sp/$sp.BCR_ABL1_all_CBtag.bam |cut -f 1 | uniq) |\
    sort | uniq > $bam_dir/$sp/$sp.BCR_ABL1_WT.reads.list

    $python3_exe $extract_bam \
    -b $bam_dir/$sp/$sp.BCR_ABL1_all_CBtag.bam \
    -n $bam_dir/$sp/$sp.BCR_ABL1_WT.reads.list \
    -o $bam_dir/$sp/$sp.BCR_ABL1_WT_CBtag.bam &&\
    $samtools_exe sort -o $bam_dir/$sp/$sp.BCR_ABL1_WT_CBtag.sort.bam
    $samtools_exe index $bam_dir/$sp/$sp.BCR_ABL1_WT_CBtag.bam

    rm -f $bam_wt_prefix.bam $bam_dir/$sp/$sp.BCR_ABL1_WT_CBtag.bam \
    $bam_wt_prefix.BCR_downstream.bam $bam_wt_prefix.BCR_downstream.bam.bai \
    $bam_wt_prefix.ABL1_upstream.bam $bam_wt_prefix.ABL1_upstream.bam.bai

    echo "Merge bam files end at: "`date +%Y-%m-%d,%H:%M:%S`    
}
