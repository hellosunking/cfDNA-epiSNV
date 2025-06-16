#!/usr/bin/perl
use strict;
use warnings;

my ($invcf, $ffile) = @ARGV;

my %h;

open my $FILE, '-|', "zcat $ffile" or die "connot open filterfile : $!";
while(my $tmp=<$FILE>){
	chomp $tmp;
	next if $tmp=~/^#/;
	my @f=split(/\t/,$tmp);
	my $key=join("\t", @f[0..1]);
	$h{$key} = 1;
	
}
close $FILE;

open my $FF, '-|', "zcat $invcf" or die "cannot open input file : $!";
while(my $line =<$FF>){
    if($line=~/^#/){
        print $line;
    }else{
        chomp $line;
        my @snp=split(/\t/,$line);
        my $keys=join("\t", $snp[0],$snp[1]);
        if(!exists $h{$keys}){
            print $line,"\n";
        }
    }
}

close $FF;
