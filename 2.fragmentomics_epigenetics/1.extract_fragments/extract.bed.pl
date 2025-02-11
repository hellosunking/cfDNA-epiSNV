#!/usr/bin/perl
use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.alleles> <output.prefix> [min.Qual=30]\n\n";
	exit 2;
}

my $minQual = $ARGV[2] || 30;
$minQual += 33;

my (%tsize, %nsize);
my ($tbed,  $nbed ) = ('', '');
my ($tcount, $ncount, $ecount, $lowQual) = ( 0, 0, 0, 0 );

open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##chr10	77205	77390	T	F	V300094330L3C005R0540237572	77363,T,G
	if( ord($l[4]) < $minQual ) {
		++ $lowQual;
		next;
	}
	my ($locus, $ref, $mut) = split /,/, $l[-1];
	if( $l[3] eq $ref ) {
		++ $ncount;
		if( $l[2] ) {
			my $s = $l[2]-$l[1];
			$nsize{$s} ++;
			$nbed .= "$l[0]\t$l[1]\t$l[2]\n";
		}
	} elsif( $l[3] eq $mut ) {
		++ $tcount;
		if( $l[2] ) {
			my $s = $l[2]-$l[1];
			$tsize{$s} ++;
			$tbed .= "$l[0]\t$l[1]\t$l[2]\n";
		}
	} else {
		++ $ecount;
	}
}

## bed
open OUT, '|-', "gzip > Mut_DNA.$ARGV[1].bed.gz" or die( "$!" );
print OUT $tbed;
close OUT;

open OUT, '|-', "gzip > Wt_DNA.$ARGV[1].bed.gz" or die( "$!" );
print OUT $nbed;
close OUT;
