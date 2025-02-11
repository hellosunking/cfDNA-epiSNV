#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV < 4) {
    die "Usage: $0 <invcf> <cosmic> <known_snp> <outvcf>\n";
}

my ($invcf, $ffilecos, $ffiledb, $outvcf) = @ARGV;

my %h;
my %d;

open my $FILEC, '-|', "zcat $ffilecos" or die "connot open filterfilecos : $!";
while(my $tmpc=<$FILEC>){
	chomp $tmpc;
	my @fc=split(/\t/,$tmpc);
        my $keyc=join("\t", @fc[0..3]);
	$h{$keyc} = 1;
	
}

close $FILEC;

open my $FILED, '-|', "zcat $ffiledb" or die "connot open filterfiledb : $!";
while(my $tmpd=<$FILED>){
        chomp $tmpd;
        my @fd=split(/\t/,$tmpd);
        my $keyd=join("\t", @fd[0..3]);
        if(!exists $h{$keyd}){
            $d{$keyd} = 2;
        }

}

close $FILED;

open my $FF, '-|', "zcat $invcf" or die "cannot open input file : $!";
open my $WF, '|-', "gzip > $outvcf" or die "cannot open outputfile : $!";
while(my $line =<$FF>){
    if($line=~/^#/){
        print $WF $line;
    }else{
        chomp $line;
        my @snp=split(/\t/,$line);
        my $keys=join("\t", $snp[0],$snp[1],$snp[3],$snp[4]);
        if(!exists $d{$keys} ){
            print $WF $line,"\n";
        }
    }
}

close $FF;
close $WF;
