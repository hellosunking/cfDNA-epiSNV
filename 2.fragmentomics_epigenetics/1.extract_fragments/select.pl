#!/usr/bin/perl
use strict;
use warnings;

my $minQual = 30;
$minQual += 33;

open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
        chomp;
        my @l = split /\t/;
        next if ord($l[4]) < $minQual;
        print join("\t", @l), "\n";
}
close IN;
