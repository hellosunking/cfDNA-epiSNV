#!/usr/bin/perl

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.fa> [motif.size=4] [discard.N=y|n]\n";
	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $size = $ARGV[1] || 4;
my $discardN = 1;
if( $#ARGV >= 2 ) {
	$discardN = 0 if $ARGV[2] =~ /^no?$/i;
}

my $all = 0;
my (%left, %right);
open IN, "$ARGV[0]" or die( "$!: $ARGV[0]" );
while( <IN> ) {
	my $seq = <IN>;
	next if $seq =~ /N/;
	chomp( $seq );
	my $a = uc substr( $seq, 0, $size );
	$left{$a} ++;
	my $b = uc substr( $seq, length($seq)-$size );
	$right{$b} ++;

	++ $all;
}
close IN;

print "#Motif\t#Left\t%Left\t#Right\t%Right\n";
foreach my $m (sort keys %left) {
	my $l = $left{$m}  || 0;
	my $r = $right{$m} || 0;
	next if $discardN && $m=~/N/i;
	print join("\t", $m, $l, $l/$all*100, $r, $r/$all*100), "\n";
}

