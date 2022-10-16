#! /usr/bin/env python

"""
Usage: used to and draw tanghulu plot
Author: Xue Xiaohui
Last modification: 20220905
"""

import os
import re
import argparse
import logging
from clubcpg.ParseBam_NOMe_indel_snv import BamFileReadParser_NOMe_indel_snv
# from plotnine import *

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-b", "--input_bam",
                        help="Input bam file, coordinate sorted with index present")
arg_parser.add_argument("-r", "--input_fa",
                        help="Input fa file, indexed with fai",
                        default="/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/xuexiaohui/database/hg38/BWA/hg38.genome.fa")
arg_parser.add_argument("-c", "--input_coordinate",
                        help="Input region to be plot")
arg_parser.add_argument("-t", "--ctype",
                        help="Select type C site to be calc, e.g. WCG, GCH,CpG")
arg_parser.add_argument("-s", "--snv",
                        help="Input coordinate of snv to be plot", default=0)
arg_parser.add_argument("-o", "--outfile",
                        help="Output file to be saved")

args = arg_parser.parse_args()
input_bam = args.input_bam
input_fa = args.input_fa
coordinate = args.input_coordinate
ctype = args.ctype
snv = args.snv
outfile = args.outfile

coord = re.split(':|-',coordinate)

parser = BamFileReadParser_NOMe_indel_snv(input_bam,input_fa, quality_score=20)

'''
Debug:
    from clubcpg.ParseBam_NOMe_indel_snv import BamFileReadParser_NOMe_indel_snv
    parser = BamFileReadParser_NOMe_indel_snv("cool_P12.IGF2R_3.bam", "/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/xuexiaohui/database/hg38/BWA/hg38.genome.fa")
    reads = parser.parse_reads("chr6",160005926,160006426,"WCG",snv=160006174)
    reads = parser.parse_reads(chromosome=coord[0],start=int(coord[1]),stop=int(coord[2]),Ctype=ctype,snv=int(snv))
    print(coord[0],int(coord[1]),int(coord[2]),ctype,snv, type(snv))
    mtx = parser.create_matrix(reads)
    mtx.to_csv("IGF2R_3.SNV.test.txt", sep="\t", na_rep='NULL')
'''

try:
    reads = parser.parse_reads(coord[0],int(coord[1]),int(coord[2]),ctype,int(snv))
    mtx = parser.create_matrix(reads)
    mtx.to_csv(outfile, sep="\t", na_rep='NULL')
except ValueError as e:
    logging.debug(str(e))
