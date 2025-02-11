#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV < 2) {
    die "Usage: $0 <input_vcf> <output_vcf>\n";
}

my ($invcf, $outvcf) = @ARGV;

open my $FF, '-|', "zcat $invcf" or die "cannot open input file : $!";
my @DP_values;

while(my $line =<$FF>){
    next if $line=~/^#/;
    chomp $line;
    my @snp=split(/\t/,$line);
    my $DP = $snp[7] =~ /DP4=([^;\s]+)/ ? $1 : "NA";
    next if $DP eq "NA";
    my @DPnew = split(/,/,$DP);
    my $DPhigh = $DPnew[0]+$DPnew[1]+$DPnew[2]+$DPnew[3];
    push @DP_values, $DPhigh;
}

close $FF;

my @sorted_DP_values = sort {$b <=> $a} @DP_values;

my $one_percent_index = int(scalar(@sorted_DP_values) * 0.01) - 1;

my $one_percent_threshold = $sorted_DP_values[$one_percent_index];

open my $FF2, '-|', "zcat $invcf" or die "cannot open input file : $!";
open my $WF, '|-', "gzip > $outvcf" or die "cannot open outputfile : $!";
while(my $line2 =<$FF2>){
    if($line2=~/^#/){
        print $WF $line2;
    }else{
    chomp $line2;
    my @snp2=split(/\t/,$line2);
    my $DP2 = $snp2[7] =~ /DP4=([^;\s]+)/ ? $1 : "NA";
    next if $DP2 eq "NA";
    my @DPnew2 = split(/,/,$DP2);
    my $DPhigh2 = $DPnew2[0]+$DPnew2[1]+$DPnew2[2]+$DPnew2[3];
    if ( $DPhigh2 < $one_percent_threshold ) {
        print $WF $line2,"\n";
    }
    } 
}

close $FF;
close $WF;
