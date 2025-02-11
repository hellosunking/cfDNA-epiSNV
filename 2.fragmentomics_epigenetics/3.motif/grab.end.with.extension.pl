#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL
#

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.bed> <output.prefix> [mode=pe|se] [size.range=all] [autosome.only=n|y] [outer.ext=2] [inner.ext=4]\n\n";
	exit 2;
}

my $layout   = $ARGV[2] || 'PE';
my $sizeRange= $ARGV[3] || 'all';
my $autoOnly = $ARGV[4] || 'n';
my $out      = $ARGV[5] || 2;
my $in       = $ARGV[6] || 4;

my ($minSize, $maxSize) = (0, 1e9);

my $paired_end = 1;

if( $layout =~ /se/i ) {
	$paired_end = 0;
	print STDERR "INFO: Single-end mode is ON.\n";
} elsif( $sizeRange ne 'all' ) {
	if( $sizeRange =~ /^(\d+),(\d+)$/ ) {
		$minSize = $1;
		$maxSize = $2;
		print STDERR "INFO: Size range: $minSize - $maxSize.\n";
	} else {
		print STDERR "ERROR: incorrect size.range! Correct parameter should be all, or 100,500 thing.";
		exit 1;
	}
}

my $autoFlag = 0;
if( $autoOnly =~ /y/i ) {
	$autoFlag = 1;
	print STDERR "INFO: Autosome-only mode is ON.\n";
}

open L, ">$ARGV[1].left.outer$out.inner$in.bed" or die( "$!" );
#open R, ">$ARGV[1].right.outer$out.inner$in.bed" or die( "$!" );

if( $ARGV[0] =~/.bz2$/ ) {
	open IN, "pbzip2 -cd $ARGV[0] |" or die( "$!: $ARGV[0]" );
} elsif( $ARGV[0] =~/.gz$/ ) {
	open IN, "gzip -cd $ARGV[0] |" or die( "$!: $ARGV[0]" );
} else {
	open IN, "$ARGV[0]" or die( "$!: $ARGV[0]" );
}

#my %count;
#my $total = 0;
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##chr start end extra

	next if $autoFlag && $l[0]!~/^chr\d+$/;

	if( $paired_end ) {
		my $size = $l[2] - $l[1];
		next if $size < $minSize || $size > $maxSize;
	} else {	## single-end data, only use those mapped to Watson strand
		my $strand;
		if( $#l >= 5 ) {
			$strand = $l[5];
		} else {
			$strand = $l[-1];
		}

		if( $strand eq '+' ) {
			## watson strand data, keep it
		} elsif( $strand eq '.' || $strand eq '-' ) {
			next;
		} else {
			print STDERR "ERROR: incorrect strand for '$_'!\n";
			next;
		}
	}

	print L join("\t", $l[0], $l[1]-$out, $l[1]+$in, "$l[1]-$l[2]"), "\n";
#	print R join("\t", $l[0], $l[2]-$in, $l[2]+$out, "$l[1]-$l[2]"), "\n";
#	++ $total;
#	$count{$l[0]} ++;
}
close IN;

close L;
#close R;

#open COUNT, ">$ARGV[1].chr.count" or die( "$!" );
#print COUNT "#chr\tcount\t%fraction\n";
#foreach my $chr ( sort keys %count ) {
#	print COUNT join("\t", $chr,$count{$chr},$count{$chr}/$total*100), "\n";
#}
#close COUNT;

