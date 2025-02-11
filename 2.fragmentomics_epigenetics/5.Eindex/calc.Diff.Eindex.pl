#!/usr/bin/perl
use strict;
use warnings;

my $sample = $ARGV[0];
my $filename = "$sample.Eindex.raw"; 
open my $fh, '<', $filename or die "can not open $filename: $!";

while (my $line1 = <$fh>) {
    chomp $line1;
    my $line2 = <$fh>;
    chomp $line2;
    
    my @fields1 = split /\t/, $line1;
    my $wt = $fields1[-1];
    
    my @fields2 = split /\t/, $line2;
    my $mut = $fields2[-1];
    
    my $diff = $mut - $wt;
    print "$sample\t$wt\t$mut\t$diff\n";
}

close $fh;

