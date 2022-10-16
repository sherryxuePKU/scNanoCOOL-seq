#!/usr/bin/perl -w
#use strict;

# NOMe.genome_WCG/GCH format
# chr	start	end	strand	triple_base	meth_reads	unmeth_reads WCG/GCH	triple_base
# chr1	22461   22462   +       ACG     1       0       WCG     ACG

my %promoter;
open IN,"$ARGV[0]" or die $!; # ../../../../annotation/hg19.Intragenic.xls
while(<IN>){
	chomp;
	my @a = split(/\t/,$_);
	my $s5 = $a[1];
	my $s3 = $a[2];
	# if ($a[4] eq "-"){
	# 	$s3 = $a[1];
	# 	$s5 = $a[2];
	# }
	$promoter{$a[0]}{$s5} = $s3;
}
close IN;

#my $DEPTH = 3;
my $DEPTH = 3; ## added by xxh
my %CpG_meth;
my %CpG_num;
open IN,"$ARGV[1]" or die $!; # ../8cell.SingleCmet
while(<IN>){
	chomp;
	next if (/#/);
	my @a = split /\s+/,$_;
	#next if ($a[4] < $DEPTH);
  next if ($a[0] eq "chrM" or $a[0] =~ /random/);
	$CpG_meth{$a[0]}{$a[1]} = $a[5]; # added by xxh for NOMe-seq report
  $CpG_num{$a[0]}{$a[1]} = $a[5] + $a[6];
  #$CpG{$a[1]}{$a[2]} = $a[5]; #methylation level, modified in 20210201 
}
close IN;

my $bin_size = 50;
# my $bin_num = 200;

foreach my $chr (sort keys %promoter){
 foreach my $s5 (sort keys %{$promoter{$chr}}){
    my $s3 = $promoter{$chr}{$s5};
    my $mid = ($s5 + $s3)/2;
    my $length = abs($s5-$s3);
    # my $bin = int($length/$bin_num);
    print "$chr\t$s5";
        for(my $pos=$mid-1000;$pos+($bin_size - 1)<=$mid;$pos+=$bin_size){
            my ($sum,$n) = (0,0);
            foreach my $i ($pos..($pos+($bin_size - 1))){
                    if (exists $CpG_meth{$chr}{$i}){
                        $sum += $CpG_meth{$chr}{$i};
                        #$n ++;
                        $n += $CpG_num{$chr}{$i};
                    }
            }
            my $mean = ($n != 0) ? $sum/$n : "NA";
            print "\t$mean";
        }
        # my $pos = $s5;
        # my ($sum,$n) = (0,0);
        #     foreach my $j ($pos..$s3){
        #             if (exists $CpG_meth{$chr}{$j}){
        #                 $sum += $CpG_meth{$chr}{$j};
        #                 #$n ++;
        #                 $n += $CpG_num{$chr}{$j};
        #             }
        #     }
        #     my $mean = ($n != 0) ? $sum/$n : "NA";
        #     print "\t$mean";
        for(my $pos=$mid;$pos+($bin_size - 1)<=$mid+1000;$pos+=$bin_size){
            my ($sum,$n) = (0,0);
            foreach my $i ($pos..($pos+($bin_size - 1))){
                    if (exists $CpG_meth{$chr}{$i}){
                            $sum += $CpG_meth{$chr}{$i};
                            #$n ++;
                            $n += $CpG_num{$chr}{$i};
                    }
            }
            my $mean = ($n != 0) ? $sum/$n : "NA";
            print "\t$mean";
        }
        print "\n";
	}
}
