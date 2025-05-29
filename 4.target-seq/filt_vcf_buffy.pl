#!/usr/bin/perl
use strict;
use warnings;

my ($invcf, $outvcf) = @ARGV;

open my $FF, '-|', "zcat $invcf" or die "cannot open input file : $!";
open my $WF, '|-', "bgzip > $outvcf" or die "cannot open outputfile : $!";
while( my $line=<$FF> ) {
	if( $line=~/^#/ ) {
		print $WF $line;
	} else {
		chomp $line;
		my @snp=split(/\t/,$line);
		my @ff = split /:/, $snp[9];
		my $DP = $ff[2];
		my @ddp = split /,/, $ff[1];
		my $altd = $ddp[1];
		my $refd = $ddp[0];
		my $finalaltd = ($altd <= $refd) ? $altd : $refd;
		my $fff = $finalaltd/$DP;
		if ($DP >= 300 && $fff >= 0.02 && $fff <= 0.3) {
			print $WF $line, "\n";
		}
	}
}

close $FF;
close $WF;
