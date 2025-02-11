#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV < 3) {
    die "Usage: $0 <input_vcf> <filter_file> <output_vcf>\n";
}

my ($invcf, $ffile, $outvcf) = @ARGV;

my %h;

open my $FILE, '<', $ffile or die "connot open filterfile : $!";
while(my $tmp=<$FILE>){
	chomp $tmp;
	my @f=split(/\t/,$tmp);
        my $key=join("\t", $f[0], $f[1]);
	$h{$key} = 1;
	
}

close $FILE;

open my $FF, '-|', "zcat $invcf" or die "cannot open input file : $!";
open my $WF, '|-', "gzip > $outvcf" or die "cannot open outputfile : $!";
while(my $line =<$FF>){
    if($line=~/^#/){
        print $WF $line;
    }else{
        chomp $line;
        my @snp=split(/\t/,$line);
        my $keys=join("\t", $snp[0],$snp[1]);
        if(!exists $h{$keys}){
            print $WF $line,"\n";
        }
    }
}

close $FF;
close $WF;
