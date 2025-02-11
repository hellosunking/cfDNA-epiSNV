#!/usr/bin/perl
use strict;
use warnings;

if( $#ARGV < 1 ) {
    print STDERR "\nUsage: $0 <in.bed.gz> <output.size>\n\n";
    exit 2;
}

my %size;

my $infile = $ARGV[0];
my $output = $ARGV[1];

open my $in, "zcat $infile |" or die "Failed to open file with zcat: $!";

while( <$in> ) {
    chomp;
    my @l = split /\t/;    ## chr10  77205  77390
    my $s = $l[2]-$l[1];
        $size{$s} ++;
}
close $in;

open OUT, ">$output" or die( "$!" );
foreach my $i ( 30..600 ) {
    print OUT join("\t", $i, $size{$i} || 0), "\n";
}
close OUT;
