#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# First version:
# Modified date:

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <in.R1.fq> <in.R2.fq> <out.prefix>\n\n";
	exit 2;
}

open O1, ">$ARGV[2].R1.fq" or die( "$!" );
open O2, ">$ARGV[2].R2.fq" or die( "$!" );


open I1, "$ARGV[0]" or die( "$!" );
open I2, "$ARGV[1]" or die( "$!" );
my ($id1, $id2, $s1, $s2, $q1, $q2);

while( $id1=<I1> ) {
	$s1 = <I1>;
	<I1>;
	$q1 = <I1>;

	$id2 = <I2>;
	$s2 = <I2>;
	<I2>;
	$q2 = <I2>;

	my $umi = substr($s1, 0, 3) . substr($s2, 0, 3);
#	next if $umi=~/N/;	## should I discard these reads?

	chomp($id1);
	$id1 =~ s/\s.*//;
	chomp($id2);
	$id2 =~ s/\s.*//;

	print O1 "$id1#$umi\n", substr($s1, 5), "+\n", substr($q1, 5);
	print O2 "$id2#$umi\n", substr($s2, 5), "+\n", substr($q2, 5);
}
close I1;
close I2;

close O1;
close O2;

