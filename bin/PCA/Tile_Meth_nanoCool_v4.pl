#!/usr/bin/perl -w
## This script is used to calculate the Tile Methylation(CpG, CHH CHG)
## Usage : perl $0 name_CpG.mfreq.txt TileSize Depth Type
## Xue Xiaohui 
## Latest Edition : 03182020

use warnings;
use strict;
#open INDEX,$ARGV[0] or die "plese input *.fa.fai";
if (@ARGV != 5 ){
	print "Usage : perl $0 <mfreq_methLevel.gz> <bin_length> <depth> <type> <outfile>\n";
	exit;
}

#open IN,$ARGV[0] or die "Please input methylKit file";## This is used for methylation site info of nanoCool WCG_meth_report
open IN,"cat $ARGV[0] | " or die $!; ## This is used for methylation site info in mfreq format
my $window = $ARGV[1]; ## bp
my $Depth = $ARGV[2] ; ## minimum reads for each site
my $Type = $ARGV[3] ; ## CpG | CHH | CHG
open OUT,">$ARGV[4]" or die "Please input the outfile name";

my %hash;
my ($start,$end);
while(<IN>){
	chomp;
	next if $_=~/chrBase/;
  next if $_=~/Lambda/;
	my @lines = split(/\s+/,$_);
	my($chr,$pos,$total,$meth,$ratio)= ($lines[0],$lines[1],$lines[5]+$lines[6],$lines[5], $lines[5]/($lines[5]+$lines[6]));
	next if ($total < $Depth);
	my $window_number = int($pos/$window);
	if(!exists $hash{$chr}{$window_number}){
		$hash{$chr}{$window_number} = "0\t0\t0\t0"; #[total_reads;meth_reads;meth_ratio_sum;CpG_number]
	}
	my @C_values = split(/\t/,$hash{$chr}{$window_number});
	my ($C_total,$C_meth,$C_ratio,$C_num) = ($C_values[0],$C_values[1],$C_values[2],$C_values[3]);
	$C_total += $total;
	$C_meth += $meth;
	$C_ratio += $ratio;
	$C_num += 1 ;
	$hash{$chr}{$window_number} = "$C_total\t$C_meth\t$C_ratio\t$C_num";

}

#print OUT "Chr\tStart\tEnd\tWindow_Number\t${Type}_total\t${Type}_meth\t${Type}_num\t${Type}_ratio\t${Type}_aveRatio\n";
my $average;
my $ave_ratio;
my $big;
foreach my $key1(sort keys %hash){
	my $key2;
	foreach $key2(sort {$a<=>$b} keys %{$hash{$key1}}){
		$big = $key2;
	}
	foreach $key2(0,1,2...$big){
		$start = 1+$window*$key2;
		$end = $window*($key2+1);
		if(exists $hash{$key1}{$key2}){
		my @words = split(/\t/,$hash{$key1}{$key2});
		my ($C_total,$C_meth,$C_ratio,$C_num) = ($words[0],$words[1],$words[2],$words[3]);
		next if ($C_total eq 0);
		$average = sprintf("%.3f" , $C_total==0 ? "-1" : $C_meth/$C_total);
		$ave_ratio = sprintf("%.3f" , $C_num==0 ? "-1" : $C_ratio/$C_num);
		#print OUT "$key1\t$start\t$end\t$key2\t$C_total\t$C_meth\t$C_num\t$average\t$ave_ratio\n";
		print OUT "$key1\t$start\t$end\t$average\n";
		}
	}
}
