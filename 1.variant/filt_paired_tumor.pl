#!/usr/bin/perl
use strict;
use warnings;

my ($invcf) = @ARGV;

open my $FF, '-|', "zcat $invcf" or die "cannot open input file : $!";
while(my $line =<$FF>){
    if($line=~/^#/){
        print $line;
    }else{
		chomp $line;
		my @snp=split(/\t/,$line);
		next if $snp[0] eq "chrY" or $snp[0] eq "chrM";
		next if $snp[4] =~ /,/;
		my @ff2 = split /:/, $snp[9];
		my @DP2 = split /,/, $ff2[1];
		my $finalDP2 = $DP2[0] + $DP2[1];
		my $finalaltd2 = ($DP2[1] <= $DP2[0]) ? $DP2[1] : $DP2[0];
		my $vaf2 = $finalaltd2 / ($DP2[1] + $DP2[0]) if $DP2[1] + $DP2[0] > 0;
		my @ff1 = split /:/, $snp[10];
		my @DP1 = split /,/, $ff1[1];
		my $finalDP1 = $DP1[0] + $DP1[1];
		my $finalaltd1 = ($DP2[1] <= $DP2[0]) ? $DP1[1] : $DP1[0];
		my $vaf1 = $finalaltd1 / ($DP1[1] + $DP1[0]) if $DP1[1] + $DP1[0] > 0;
		if ($finalDP1 >= 30 && $finalDP2 >=30 && $vaf2 < 0.01 && $finalaltd1 >= 3) {
			print $line, "\n";
		}
	}
}

close $FF;
