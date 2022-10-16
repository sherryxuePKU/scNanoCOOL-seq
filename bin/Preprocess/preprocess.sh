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

### Function
## Step1: trim barcode and filter
function do_trim() {
    echo "=================================="
    echo "Trim barcode begin at: "`date +%Y-%m-%d,%H:%M:%S`
    mkdir -p $trim_dir
    mkdir -p $indir/QC
    input=$indir/$sp.fastq
    tmp=$trim_dir/$sp.clean.tmp.fastq
    output=$trim_dir/$sp.clean.fastq
    NanoStat_fastq=$indir/QC/$sp.raw.nanostat_report.txt
    bc=`echo $sp | awk -F "[_.]" '{print $3}'`
    bc_seq=`cat $barcode_list | awk -v a=$bc '$1==a {print $2}'`
    #bc_seq_rev=`echo $bc_seq | tr a-z A-Z|tr ATCG TAGC | rev`
    ## trim primer
    $cutadapt_exe -j 4 -g $bc_seq --revcomp -e 0.2 -m 100 -o $tmp $input
    ## filter length < 100
    $seqtk_exe trimfq -b 6 -e 9 $tmp > $output &&\
    rm -f $tmp
    
    $nanostat_exe --fastq $indir/$sp.fastq -t 4 >$NanoStat_fastq

    # gzip $input $output
    echo "Trim barcode end at: "`date +%Y-%m-%d,%H:%M:%S`
}

## Step2: bismark_mm2 branch (6 core, 100G, 8h)
function do_align(){
    echo "=================================="
    echo "Align by bismark_minimap2 begin at: "`date +%Y-%m-%d,%H:%M:%S`
    mkdir -p $bam_dir/$sp
    trim=$trim_dir/$sp.clean.fastq.gz
    bam=$bam_dir/$sp/$sp.clean_bismark_mm2
    $bismark_mm2 --mm2 --parallel 1 --non_directional \
    --fastq --mm2_maximum_length 100000 \
    --path_to_minimap2 $minimap2_dir \
    --samtools_path $bin_dir/samtools-1.9/bin \
    --output_dir $bam_dir/$sp \
    --temp_dir $bam_dir/$sp \
    $bismark_mm2_db $trim &&\
    $samtools_exe view -uSb -@ 3\
    -t $bismark_mm2_db/hg38.genome_lambda.fa $bam.bam |\
    $samtools_exe sort -m 900M -@ 3 -T $sp -o $bam.sort.bam &&\
    $samtools_exe rmdup -s $bam.sort.bam $bam.sort.rmdup.bam &&\
    $samtools_exe index $bam.sort.rmdup.bam &&\
    rm -f $bam.sort.bam
    echo "Align by bismark_minimap2 end at: "`date +%Y-%m-%d,%H:%M:%S`
}

## Step3: read-level quality control
function do_filter_read_mCH(){
    echo "=================================="
    echo "Filter read-level CH methylation begin at: "`date +%Y-%m-%d,%H:%M:%S`
    bam_prefix=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup
    filter_report=$bam_dir/$sp/$sp.mCH_filter.report.txt

    $methylpy_exe bam-quality-filter \
    --input-file $bam_prefix.bam \
    --output-file $bam_prefix.CHlt80per \
    --ref-fasta $ref --min-num-ch 1 --max-mch-level 0.8 &&\
    $samtools_exe index $bam_prefix.CHlt80per.bam 

    ## extract reads mapped to lambda
    $samtools_exe view -bSh $bam_prefix.bam Lambda -o $bam_prefix.lambda.bam &&\
    $samtools_exe index $bam_prefix.lambda.bam

    echo -e -n "Sample\tBefore_filter\tAfter_filter\n" > $filter_report
    before=`$samtools_exe view -F 2308 -c $bam_prefix.bam`
    after=`$samtools_exe view -F 2308 -c $bam_prefix.CHlt80per.bam`
    echo -e -n "$sp\t$before\t$after\n" >> $filter_report
    echo "Filter read-level CH methylation end at: "`date +%Y-%m-%d,%H:%M:%S`  
}

## Step4: call methylation by methylpy (1-based)
function do_methylpy_meth_extractor(){
    echo "=================================="
    echo "Extract DNA methylation site begin at: "`date +%Y-%m-%d,%H:%M:%S`
    bam=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup
    mkdir -p $meth_dir/$sp

    $methylpy_exe call-methylation-state \
    --input-file $bam.bam --sample $sp.CHlt80per \
    --ref-fasta $ref --paired-end False \
    --num-upstream-bases 1 --num-downstream-bases 1 \
    --num-procs 4 --path-to-output $meth_dir/$sp

    zcat $meth_dir/$sp/allc_${sp}.CHlt80per.tsv.gz |\
    awk -v OFS="\t" '{if($4 ~ /[AT]CG/ && $1!="lambda")print "chr"$1,$2-1,$2,$3,$4,$5,$6-$5,"WCG",$4}' > $meth_dir/$sp/$sp.NOMe.genome.WCG_report.txt
    zcat $meth_dir/$sp/allc_${sp}.CHlt80per.tsv.gz |\
    awk -v OFS="\t" '{if($4 ~ /[AT]CG/ && $1=="lambda")print "chr"$1,$2-1,$2,$3,$4,$5,$6-$5,"WCG",$4}' > $meth_dir/$sp/$sp.NOMe.lambda.WCG_report.txt
    
    zcat $meth_dir/$sp/allc_${sp}.CHlt80per.tsv.gz |\
    awk -v OFS="\t" '{if($4 ~ /GC[ATC]/ && $1!="lambda")print "chr"$1,$2-1,$2,$3,$4,$5,$6-$5,"GCH",$4}' > $meth_dir/$sp/$sp.NOMe.genome.GCH_report.txt
    zcat $meth_dir/$sp/allc_${sp}.CHlt80per.tsv.gz |\
    awk -v OFS="\t" '{if($4 ~ /GC[ATC]/ && $1=="lambda")print "chr"$1,$2-1,$2,$3,$4,$5,$6-$5,"GCH",$4}' > $meth_dir/$sp/$sp.NOMe.lambda.GCH_report.txt
    
    awk -v OFS="\t" '{if($1!="chrLambda")print $1,$2,$3,$6/($6+$7)}' \
    $meth_dir/$sp/$sp.NOMe.genome.WCG_report.txt |\
    sort -k1,1 -k2,2n > $meth_dir/$sp/$sp.NOMe.genome.WCG_meth.tmp.bedgraph &&\
    $bedGraphToBigWig $meth_dir/$sp/$sp.NOMe.genome.WCG_meth.tmp.bedgraph \
    $chrom_size $meth_dir/$sp/$sp.NOMe.genome.WCG_meth.bw

    awk -v OFS="\t" '{if($1!="chrLambda")print $1,$2,$3,$6/($6+$7)}' \
    $meth_dir/$sp/$sp.NOMe.genome.GCH_report.txt |\
    sort -k1,1 -k2,2n > $meth_dir/$sp/$sp.NOMe.genome.GCH_meth.tmp.bedgraph &&\
    $bedGraphToBigWig $meth_dir/$sp/$sp.NOMe.genome.GCH_meth.tmp.bedgraph \
    $chrom_size $meth_dir/$sp/$sp.NOMe.genome.GCH_meth.bw

    rm -f $meth_dir/$sp/$sp.NOMe.genome.WCG_meth.tmp.bedgraph \
    $meth_dir/$sp/$sp.NOMe.genome.GCH_meth.tmp.bedgraph
    
    echo "Extract DNA methylation site end at: "`date +%Y-%m-%d,%H:%M:%S`
}

## Step5-1: stat quality control-conversion ratio
function do_methylpy_CT_CA_efficiency_v2(){
    echo "=================================="
    echo "Calculate CT_conversion_ratio and CA_effiency begin at: "`date +%Y-%m-%d,%H:%M:%S`
    bam=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.lambda

    $methylpy_exe call-methylation-state \
    --input-file $bam.bam --sample $sp.lambda \
    --ref-fasta $ref --paired-end False \
    --num-upstream-bases 1 --num-downstream-bases 1 \
    --num-procs 4 --path-to-output $meth_dir/$sp

    meth_summary=$meth_dir/$sp/$sp.methylpy_meth_summary.v3.txt
    mkdir -p $meth_dir/$sp
    echo -e -n "Sample\tCT_conversion_WCG\tCT_conversion_WCH\tCA_conversion\t" > $meth_summary
    CT_conv_WCG=`zcat $meth_dir/$sp/allc_${sp}.lambda.tsv.gz | awk '{if($1=="Lambda" && $4 ~ /[AT]CG/){ratio+=1-$5/$6;num+=1}}END{print ratio/num}'`
    CT_conv_WCH=`zcat $meth_dir/$sp/allc_${sp}.lambda.tsv.gz | awk '{if($1=="Lambda" && $4 ~ /[AT]C[ATC]/){ratio+=1-$5/$6;num+=1}}END{print ratio/num}'`
    CA_conv_GCH=`zcat $meth_dir/$sp/allc_${sp}.lambda.tsv.gz | awk '{if($1=="Lambda" && $4 ~ /GC[ATC]/){ratio+=$5/$6;num+=1}}END{print ratio/num}'`
    
    echo -e -n "CT_conversion_WCG.CHlt80\tCT_conversion_WCH.CHlt80\tCA_conversion.CHlt80\n" >> $meth_summary
    CT_conv_WCG_filtered=`zcat $meth_dir/$sp/allc_${sp}.CHlt80per.tsv.gz | awk '{if($1=="Lambda" && $4 ~ /[AT]CG/){ratio+=1-$5/$6;num+=1}}END{print ratio/num}'`
    CT_conv_WCH_filtered=`zcat $meth_dir/$sp/allc_${sp}.CHlt80per.tsv.gz | awk '{if($1=="Lambda" && $4 ~ /[AT]C[ATC]/){ratio+=1-$5/$6;num+=1}}END{print ratio/num}'`
    CA_conv_GCH_filtered=`zcat $meth_dir/$sp/allc_${sp}.CHlt80per.tsv.gz | awk '{if($1=="Lambda" && $4 ~ /GC[ATC]/){ratio+=$5/$6;num+=1}}END{print ratio/num}'`
    echo -e -n "$sp\t$CT_conv_WCG\t$CT_conv_WCH\t$CA_conv_GCH\t" >> $meth_summary
    echo -e -n "$CT_conv_WCG_filtered\t$CT_conv_WCH_filtered\t$CA_conv_GCH_filtered\n" >> $meth_summary

    echo "Calculate CT_conversion_ratio and CA_effiency end at: "`date +%Y-%m-%d,%H:%M:%S`
}

## Step5-2 stat quality control-compute matrix for GCH_TSS profile
function do_meth_profile_qc(){
    echo "=================================="
    echo "Collect TSS or Genebody Meth begin at: "`date +%Y-%m-%d,%H:%M:%S`
    query=$1
    region=$2

    case $region in
	"GENE") ref=$genebody_v2_ref;;
    "TSS") ref=$TSS_ref;;
	esac

    input=$meth_dir/$sp/$sp.NOMe.genome.${query}_report.txt
    outdir=$info_dir/meth_profile_data/${query}_${region}
    output=$outdir/${sp}.${query}_${region}_avemethyrate.bed
    # p=$profile_dir/QC.$region/${sp}_${query}_in${region}_methylation.pdf

    mkdir -p $info_dir/meth_profile_data/${query}_${region}
    $bedtools_exe intersect -a $ref -b $input -wa -wb |\
    awk '{OFS="\t";print $1,$2,$3,$4,$5,$12/($12+$13)}' |\
    $bedtools_exe groupby -i - -g 1,2,3,4,5 -c 6 -o mean > $output

    #mkdir -p $profile_dir/QC.$region
    #$rscript_exe $cal_genebody $mean_bed $p
    echo "Collect TSS or Genebody Meth begin at: "`date +%Y-%m-%d,%H:%M:%S`
}

## Step5-3: stat quality control-compute matrix for WCG_genebody profile
function do_computeMatrix(){
	echo "=================================="
	echo "ComputeMatrix begin at: "`date +%Y-%m-%d,%H:%M:%S`
	region=$2
    C_type=$1
	input=$meth_dir/$sp/${sp}.NOMe.genome.${C_type}_report.txt
    outdir=$info_dir/meth_profile_data/${C_type}_${region}
	abs_mtx=$outdir/${sp}.${C_type}_${region}.absoluteDistance.txt
    abs_df=$outdir/${sp}.${C_type}_${region}.absoluteDistance_merge.txt
	
    mkdir -p $outdir
	case $region in
	"genebody") ref=$genebody_ref
	;;
	"CGI") ref=$CGI_ref
	;;
	esac
	
	$perl_exe $AbD_pl $ref $input > $abs_mtx &&\
    $rscript_exe $merge_mtx_r $abs_mtx $abs_df 

    rm -f $abs_mtx
	echo "ComputeMatrix end at: "`date +%Y-%m-%d,%H:%M:%S`
}

## Step5-4: stat quality control-basic info
function do_basic_stat(){
    echo "=================================="
    echo "Basic stat begin at: "`date +%Y-%m-%d,%H:%M:%S`
    mkdir -p $info_dir/basic_stat
    stat_report=$info_dir/basic_stat/$sp.basic_stat.txt
    bismark_mm2_report=$bam_dir/$sp/$sp.clean_bismark_mm2_SE_report.txt
    NanoStat_fastq=$indir/QC/$sp.raw.nanostat_report.txt
    NanoStat_bam=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.nanostat_report.txt
    meth=$meth_dir/$sp/allc_${sp}.CHlt80per.tsv.gz
    mCHfilter_report=$bam_dir/$sp/$sp.mCH_filter.report.txt
    bam=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.bam
    bam_mCHlt80=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.CHlt80per.bam
    NanoStat_mCHlt80_bam=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.mCHlt80.nanostat_report.txt

    echo -e -n "Sample\tTotal_reads\tTotal_reads(M)\tTotal_bases\tTotal_bases(Gb)\tRaw_median_length\tRaw_median_quality\t" >$stat_report
    echo -e -n "Clean_reads\tAligned_reads\tMapping_ratio\t" >>$stat_report
    echo -e -n "Dedupped_reads\tDedupped_reads(M)\tDuplication" >>$stat_report
    echo -e -n "mCHlt80_read\tcmCHlt80_filter_ratio\t" >>$stat_report
    echo -e -n "Final_bases\tFinal_bases(Gb)\tYield\tCoverage\t" >>$stat_report
    echo -e -n "WCG_site\tWCG_coverage\tWCG_meth\tGCH_site\tGCH_coverage\tGCH_meth\t" >>$stat_report
    echo -e -n "Conversion_ratio\tCA_conv_efficiency\tConversion_ratio_mCHlt80\tCA_conv_efficiency_mCHlt80\t" >>$stat_report
    echo -e -n "Median_length\tMedian_quality\tLongest_reads_length\t" >>$stat_report
    echo -e -n "CGI_coverage_fulllen\tPromoter_coverage_fulllen\tCGI_coverage_90len\tPromoter_coverage_90len\t" >> $stat_report
    echo -e -n "Promoter_coverage_3CpG\tCGI_coverage_3CpG\n" >> $stat_report

    ## sample name
    echo -e -n "$sp\t" >> $stat_report

    ## total_reads, total_bases, raw_median_length, raw_median_quality
    #$nanostat_exe --fastq $indir/$sp.fastq.gz -t 2 >$NanoStat_fastq
    total_read=`cat $NanoStat_fastq | grep "Number of reads:" | awk -F ":" '{print $2}' | sed 's/[ \t,]//g'`
    total_read_M=$(printf "%.2f" `echo "scale=2;$total_read/1000000"|bc`)
    total_base=`cat $NanoStat_fastq  | grep "Total bases:" | awk -F ":" '{print $2}' | sed 's/[ \t,]//g'`
    total_base_G=$(printf "%.2f" `echo "scale=2;$total_base/1000000000"|bc`)
    raw_median_length=`cat $NanoStat_fastq  | grep "Median read length:" | awk -F ":" '{print $2}' | sed 's/[ \t,]//g'`
    raw_median_quality=`cat $NanoStat_fastq  | grep "Median read quality:" | awk -F ":" '{print $2}' | sed 's/[ \t,]//g'`
    echo -e -n "$total_read\t$total_read_M\t$total_base\t$total_base_G\t$raw_median_length\t$raw_median_quality\t" >>$stat_report

    ## clean_reads, aligned_reads, mapping_ratio from bismark report
    clean_read=`cat $bismark_mm2_report | grep "Sequences analysed in total:" | awk -F ":" '{print $2}' | sed 's/[ \t]//g'`
    aligned_reads=`cat $bismark_mm2_report | grep "Number of alignments with a unique best hit from the different alignments:" | awk -F ":" '{print $2}' | sed 's/[ \t]//g'`
    mapping_ratio=`cat $bismark_mm2_report | grep "Mapping efficiency:" | awk -F ":" '{print $2}' | sed 's/[% \t]//g' | awk '{print $1/100}'` 
    echo -e -n "$clean_read\t$aligned_reads\t$mapping_ratio\t" >>$stat_report

    ## Dedupped_reads, Final_bases, Depth, Coverage, Yield
    $nanostat_exe -t 2 --bam $bam --no_supplementary -n $NanoStat_bam
    dedupped_read=`cat $NanoStat_bam | grep "Number of reads:" | awk -F ":" '{print $2}' | sed 's/[ \t,]//g'`
    dedupped_read_M=$(printf "%.2f" `echo "scale=2;$dedupped_read/1000000"|bc`)
    duplication=$(printf "%.2f" `echo "scale=2;1-${dedupped_read}/${aligned_reads}"|bc`)
    echo -n -e "$dedupped_read\t$dedupped_read_M\t${duplication}\t" >>$stat_report

    ## mCH80 filter ratio 
    mCHlt80_read=`sed 1d $mCHfilter_report | awk '{print $3}'`
    mCH_filter_ratio=$(printf "%.2f" `echo "scale=2;1-${mCHlt80_read}/${dedupped_read}"|bc`)
    echo -e -n "$mCHlt80_read\t$mCH_filter_ratio\t" >>$stat_report

    ## Final used bases
    $nanostat_exe -t 2 --bam $bam_mCHlt80 --no_supplementary -n $NanoStat_mCHlt80_bam
    final_base=`cat $NanoStat_mCHlt80_bam  | grep "Total bases:" | awk -F ":" '{print $2}' | sed 's/[ \t,]//g'`
    final_base_G=$(printf "%.2f" `echo "scale=2;$final_base/1000000000"|bc`)
    yield=$(printf "%.2f" `echo "scale=2;${final_base}/${total_base}"|bc`)
    echo -n -e "$final_base\t$final_base_G\t$yield\t" >>$stat_report
    
    ## Genomic coverage
    $samtools_exe depth $bam_mCHlt80 | awk -v size=$hg38_genome_size '{OFS="\t";if($1!="Lambda"){pos+=1}}END{printf pos/size"\t"}' >>$stat_report

    ## WCG_site, WCG_meth, GCH_site, GCH_meth
    WCG=`zcat $meth | awk -v size=$WCG_site_num '{if($4 ~ /[AT]CG/){ratio+=$5/$6;num+=1}}END{printf num"\t"num/size"\t"ratio/num}'`
    GCH=`zcat $meth | awk -v size=$GCH_site_num '{if($4 ~ /GC[ATC]/){ratio+=$5/$6;num+=1}}END{printf num"\t"num/size"\t"ratio/num}'`
    echo -e -n "$WCG\t$GCH\t" >>$stat_report

    ## conversion ratio, in vitro methylation effiency
    conv=`sed '1d' $meth_dir/$sp/$sp.methylpy_meth_summary.v3.txt | awk -v OFS="\t" '{print $2,$4,$5,$7}'` >>$stat_report
    echo -e -n "$conv\t" >>$stat_report

    ## Median_length, Median_quality,Longest reads from NanoStat_bam
    median_length=`cat $NanoStat_mCHlt80_bam  | grep "Median read length:" | awk -F ":" '{print $2}' | sed 's/[ \t,]//g'`
    median_quality=`cat $NanoStat_mCHlt80_bam  | grep "Median read quality:" | awk -F ":" '{print $2}' | sed 's/[ \t,]//g'`
    longest_len=`cat $NanoStat_mCHlt80_bam | sed -n '/Top 5 longest reads and their mean basecall quality score/{n;p}' | awk '{print $2}'`
    echo -e -n "$median_length\t$median_quality\t$longest_len\t" >>$stat_report

    ## CGI/Promoter(-500bp, 250bp) coverage (=100% length)
    bed=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.bed
    $bedtools_exe bamtobed -i $bam_mCHlt80 | grep -v "Lambda" > $bed 

    coverage_CGI=`$bedtools_exe intersect -a $CGI_ref -b $bed -wa -wb -f 1 | $bedtools_exe groupby -i - -g 1,2,3,4 -c 9 -o count | wc -l `
    coverage_promoter=`$bedtools_exe intersect -a $promoter_bed -b $bed -wa -wb -f 1 | $bedtools_exe groupby -i - -g 1,2,3,4 -c 9 -o count | wc -l `
    echo -e -n "$coverage_CGI\t$coverage_promoter\t" >> $stat_report

    ## CGI/Promoter(-500bp, 250bp) coverage (>90% length)
    coverage_CGI_90=`$bedtools_exe intersect -a $CGI_ref -b $bed -wa -wb -f 0.90 |$bedtools_exe groupby -i - -g 1,2,3,4 -c 9 -o count | wc -l `
    coverage_promoter_90=`$bedtools_exe intersect -a $promoter_bed -b $bed -wa -wb -f 0.90 | $bedtools_exe groupby -i - -g 1,2,3,4 -c 9 -o count | wc -l `
    echo -e -n "$coverage_CGI_90\t$coverage_promoter_90\t" >> $stat_report

    ## CGI/Promoter(-500bp,250bp) coverage (>=3CpG)
    CGI_CpG_bed=$bam_dir/$sp/$sp.promoter_CpG_coverage.tmp.bed
    $bedtools_exe intersect -a $CGI_CpG -b $bed -wa -wb | cut -f 4,9 | sort -k 1 -k 2 | $bedtools_exe groupby -i - -g 1,2 -c 2 -o count >$CGI_CpG_bed 
    CGI_3CpG=`awk '{if($3>=3)print $0}' $CGI_CpG_bed | $bedtools_exe groupby -i - -g 1 -c 2 -o count | awk 'END{print NR}'`
    echo -e -n "$CGI_3CpG\t" >> $stat_report && rm -f $CGI_CpG_bed

    promoter_CpG_bed=$bam_dir/$sp/$sp.promoter_CpG_coverage.tmp.bed
    $bedtools_exe intersect -a $promoter_CpG -b $bed -wa -wb | cut -f 4,9 | sort -k 1 -k 2 | $bedtools_exe groupby -i - -g 1,2 -c 2 -o count >$promoter_CpG_bed 
    promoter_3CpG=`awk '{if($3>=3)print $0}' $promoter_CpG_bed | $bedtools_exe groupby -i - -g 1 -c 2 -o count | awk 'END{print NR}'`
    echo -e -n "$promoter_3CpG\n" >> $stat_report && rm -f $promoter_CpG_bed

    echo "Basic stat end at: "`date +%Y-%m-%d,%H:%M:%S`  
}

## Step5-4(optional): plot single sample read length distribution
function do_read_length(){
	echo "=================================="
	echo "ComputeMatrix begin at: "`date +%Y-%m-%d,%H:%M:%S`
    bed=$bam_dir/$sp/$sp.clean_bismark_mm2.sort.rmdup.bed
    mkdir -p $info_dir/read_length

    awk '{print $3-$2}' $bed > $info_dir/read_length/$sp.read_len.txt
    echo "ComputeMatrix end at: "`date +%Y-%m-%d,%H:%M:%S`
}

