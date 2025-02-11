#!/usr/bin/perl
#
# Author: Ahfyth
#

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <genome.fa> <in.bed> [in.bed ...]\n\n";
	exit 2;
}

my $fasta = shift;
my $g = load_genome( $fasta );

my $index = 0;
foreach my $ifile ( @ARGV ) {
	++ $index;

	open IN, "$ifile" or die( "$!" );
	my $ofile = $ifile;
	$ofile =~ s/bed$/fa/;
	open OUT, ">$ofile" or die( "$!" );
	print STDERR "\rLoading $index: $ifile => $ofile ... ";
	my $cnt = 0;
	while( <IN> ) {
		chomp;
		my @l = split /\t/;	##chr start(0-base) end(1-base) id score strand
		next unless exists $g->{$l[0]};
		my $id = "$l[0]:$l[1]:$l[2]";
		$id .= ":$l[5]" if $#l>=5;
		print OUT ">$id\n", substr($g->{$l[0]}, $l[1], $l[2]-$l[1]), "\n";
		++ $cnt;
	}
	close IN;
	close OUT;
	print STDERR "Done: $cnt lines written.";
}
print STDERR "\rDone: Totally $index BED files loaded successfully.\n";

sub load_genome {
	my $fasta = shift;
	my %g;

	my $chr = 'EMPTY';
	open FA, "$fasta" or die ($!);
	while( <FA> ) {
		if( /^>(\S+)/ ) {
			$chr = $1;
			$g{$chr} = '';
		} else {
			chomp;
			$g{$chr} .= uc $_;
		}
	}
	close FA;

	return \%g;
}

