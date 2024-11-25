#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
#use KSLIB::loadGenome qw/loadGenome/;
#use KSLIB::cigarUtil qw/fix_seq_from_CIGAR extract_size_from_CIGAR/;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.alleles> <output.prefix> [min.Qual=20]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $minQual = $ARGV[2] || 30;
$minQual += 33;	## phred 33

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
			$nbed .= "$l[0]\t$l[1]\t$l[2]\t$l[5]\n";
		}
	} elsif( $l[3] eq $mut ) {
		++ $tcount;
		if( $l[2] ) {
			my $s = $l[2]-$l[1];
			$tsize{$s} ++;
			$tbed .= "$l[0]\t$l[1]\t$l[2]\t$l[5]\n";
		}
	} else {
		++ $ecount;
	}
}
close IN;

## statistics
open OUT, ">$ARGV[1].allele.stat" or die( "$!" );
print OUT "File\t$ARGV[0]\n",
		  "Total\t", $tcount+$ncount+$ecount+$lowQual, "\n",
		  "Tumor\t", $tcount, "\n",
		  "Normal\t",$ncount, "\n",
		  "Error\t", $ecount, "\n",
		  "LowQual\t", $lowQual, "\n";
close OUT;

## size
open OUT, ">$ARGV[1].tumor.size" or die( "$!" );
foreach my $i ( 30..600 ) {
	print OUT join("\t", $i, $tsize{$i}||0), "\n";
}
close OUT;

open OUT, ">$ARGV[1].normal.size" or die( "$!" );
foreach my $i ( 30..600 ) {
	print OUT join("\t", $i, $nsize{$i}||0), "\n";
}
close OUT;

## bed
open OUT, ">$ARGV[1].tumor.bed" or die( "$!" );
print OUT $tbed;
close OUT;

open OUT, ">$ARGV[1].normal.bed" or die( "$!" );
print OUT $nbed;
close OUT;

