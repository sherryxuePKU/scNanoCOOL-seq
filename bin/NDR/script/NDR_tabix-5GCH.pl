#!/usr/bin/perl -w
use strict;

open IN,"$ARGV[0]" or die $!;

while(<IN>){
	chomp;
	my @a = split /\s+/,$_;
		my $start=$a[1];
		my $end  =$a[2];
		my $tmp= `/home/tangfuchou_pkuhpc/lustre1/zhengyuxuan/Software/tabix-0.2.6/tabix $ARGV[1] $a[0]:$start-$end|wc -l`;
		if ( $tmp >= 5 ){
			print "$a[0]\t$start\t$end\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$tmp";
		}
	}
close IN;
 
