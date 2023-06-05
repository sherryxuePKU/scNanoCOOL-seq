import pysam
import pandas as pd
from collections import defaultdict
import logging
import re
import fnmatch
from functools import reduce

class HaploTag:
    """
    Usage: tag reads to haplotypes by phased SNVs

    :Example:
        >>> from clubcpg.ParseBam_haplotag import HaploTag
        >>> parser = HaploTag("cool_P12.CHlt80per.clean_bismark_mm2.sort.rmdup.bam")
        >>> reads = parser.parse_reads("chr7", 21420802, 'T', 'C'),
    """

    def __init__(self, bamfile):
        """
        Class used to read scNanoCOOL-seq reads from a BAM file, extract methylation, and convert into data frame

        :param bamfile: Path to bam file location
        """

        self.bamfile = bamfile

        self.OpenBamFile = pysam.AlignmentFile(bamfile, 'rb')
        # Check for presence of index file
        index_present = self.OpenBamFile.check_index()
        if not index_present:
            raise FileNotFoundError("BAM file index is not found. Please create it using samtools index")

    ## modified from SNPsplit subroutine score_bisulfite_SNPs()
    # @staticmethod
    def assign_tag(self, read_base: str, ref: str, snp: str, strand: str):
        if ref=="C":
            if snp=="T":
                if strand=="OT": hp_tag="unassigned"
            if strand=="OT":
                if read_base=="C" or read_base=="T": hp_tag="genome1"
                elif read_base==snp: hp_tag="genome2"
                else: hp_tag="unassigned"
            elif strand=="OB":
                if read_base==ref: hp_tag="genome1"
                elif snp=="G":
                    if read_base=="G" or read_base=="A": hp_tag="genome2"
                    else: hp_tag="unassigned"
                elif read_base==snp: hp_tag="genome2"
                else: hp_tag="unassigned"
        elif snp=="C":
            if ref=="T":
                if strand=="OT": hp_tag="unassigned"
            if strand=="OT":
                if read_base=="C" or read_base=="T": hp_tag="genome2"
                elif read_base==ref: hp_tag="genome1"
                else: hp_tag="unassigned"
            elif strand=="OB":
                if read_base==snp: hp_tag="genome2"
                elif ref=="G":
                    if read_base=="G" or read_base=="A": hp_tag="genome1"
                    else: hp_tag="unassigned"
                elif read_base==ref: hp_tag="genome1"
                else: hp_tag="unassigned"
        elif ref=="G":
            if snp=="A":
                if strand=="OB": hp_tag="unassigned"
            if strand=="OT":
                if read_base==ref: hp_tag="genome1"
                elif snp=="C":
                    if read_base=="C" or read_base=="T": hp_tag="genome2"
                    else: hp_tag="unassigned"
                elif read_base==snp: hp_tag="genome2"
                else: hp_tag="unassigned"
            elif strand=="OB":
                if read_base=="G" or read_base=="A": hp_tag="genome1"
                elif read_base==snp: hp_tag="genome2"
                else: hp_tag="unassigned"
        elif snp=="G":
            if ref=="A":
                if strand=="OB": hp_tag="unassigned"
            if strand=="OT":
                if base==snp: hp_tag="genome2"
                elif ref=="C":
                    if read_base=="C" or read_base=="T": hp_tag="genome1"
                    else: hp_tag="unassigned"
                elif read_base==ref: hp_tag="genome1"
                else: hp_tag="unassigned"
            elif strand=="OB":
                if read_base=="G" or read_base=="A": hp_tag="genome2"
                elif read_base==ref: hp_tag="genome1"
                else: hp_tag="unassigned"
        else:
            if read_base==ref: hp_tag="genome1"
            elif read_base==snp: hp_tag="genom2"
            else: hp_tag="unassigned"
        
        return hp_tag
    

    # Get reads from the bam file, extract methylation state
    def parse_reads(self, chromosome: str, pos:int, ref: str, snp: str):
        """
        :param chromosome: chromosome as "chr6"
        :param pos: coordinate of SNVs
        :return: List of reads and their haplotyping tags
        """
        reads = []
        for read in self.OpenBamFile.fetch(chromosome, pos, pos+1):
            reads.append(read)
        
        read_snps = []
        for read in reads:

            ## complementary bases
            comp_dict={"A":"T","T":"A","C":"G","G":"C"}

            pair=read.get_aligned_pairs(with_seq=True)

            ## (87, 21420853, 'C') 2nd element is 0-based(https://github.com/pysam-developers/pysam/issues/1180)
            pair_snp=reduce(lambda acc, x: acc + [x] if x[1]==(pos-1) else acc, pair, [])
            pair = pair_snp[0]

            ## SNPs located near indel
            if pair[0]==None: continue

            xg_tag=read.get_tag("XG")
            read_pos=pair[0]
            ref_pos=pair[1]
            read_base=pair[2]

            if xg_tag=="CT": strand="OT"
            if xg_tag=="GA": strand="OB"

            ## it seems that pysam(v0.19.0) have converted the bases in reverse strand by complement
            if read_base.islower(): true_base=read.query_sequence[read_pos]
            if read_base.isupper(): true_base=read_base

            # if read.flag==0:
            #     if read_base.islower(): true_base=read.query_sequence[read_pos]
            #     if read_base.isupper(): true_base=read_base
            # if read.flag==16:
            #     if read_base.islower(): true_base=comp_dict[read.query_sequence[read_pos]]
            #     if read_base.isupper(): true_base=comp_dict[read_base]

            hp_tag=self.assign_tag(true_base, ref, snp, strand)
            qname=read.query_name

            read_snps.append((qname, hp_tag))

        self.full_reads = reads
        self.read_snps = read_snps
        
        return read_snps
