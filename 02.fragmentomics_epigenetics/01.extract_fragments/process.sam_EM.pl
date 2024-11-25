#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
#

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <in.sam> <target.loci=chr:pos> [extra.info=pos] [mode=0|1]\n";
	print STDERR "\nBy default, I will output all fragments that covers the site.";
	print STDERR "\nWhen set mode=1, only PE reads with non-N will be kept.\n\n";
	exit 1;
}

my ($t_chr, $t_pos) = split /:/, $ARGV[1];
my $extra = $ARGV[2] || $t_pos;
my $mode  = $ARGV[3] || 0;	## default is relax-mode

my $MIN_ALIGN_SCORE = 60;

## PE/SE will be automatically detected
## load data
my %raw;
my %PEinfo;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##555626_R1	99	chr1	54339	23	36M	= 54579 276 TCAATTAAGTCC HHHHHHH XG:Z:GA
	next unless $l[2] eq $t_chr;
	next unless $l[4] >= $MIN_ALIGN_SCORE;
        #next unless $l[1] & 0x02;
        #next if     $l[1] & 0x100;
        #next unless $l[6] eq "=";
	$l[0] =~ /[#\s].*$/;##fix ID

	my $pos = $l[3];
        my $cover = $l[-1];
	my ($mseq, $mqual);

	## fix sequence using CIGAR
	my $cigar = $l[5];

	## fix clips
	$cigar =~ s/^\d+H//;
	$cigar =~ s/\d+H$//;
	next if $cigar =~ /S/;
	if( $cigar !~ /^\d+M$/ ) {	## there are indels here
		if( $cigar !~ /^[\dMIDN]+$/ ) {
			print STDERR "ERROR: unsupported CIGAR $cigar!\n";
			next;
		}

		$cigar =~ s/([MIDN])/$1:/g;
		my @info = split /:/, $cigar;
		($mseq, $mqual) = ('', '');	## modified seq and qual
		my $curr = 0;
		foreach my $m ( @info ) {
			if( $m =~ /^(\d+)I$/ ) {	## insertion: CAUTION!!! I WILL DISCARD IT!!!
				$curr += $1;
			} elsif ( $m =~ /^(\d+)[DN]$/ )	{ ## intron/deletion: CAUTION!!! I WILL ADD 'D' to indicate deletion !!!
				$mseq .= 'D' x $1;
				$mqual.= '!' x $1;
			} elsif( $m =~ /^(\d+)M$/ ) {
				$mseq .= substr($l[9],  $curr, $1);
				$mqual.= substr($l[10], $curr, $1);
				$curr += $1;
			}
		}
	} else {
		$mseq .=  $l[9];
		$mqual.= $l[10];
	}

#	print STDERR "$pos\t$mseq\t$mqual\n";
	push @{$raw{$l[0]}}, "$pos\t$mseq\t$mqual\t$cover";
	$PEinfo{$l[0]}=1 if $l[6] eq '=';
}
close IN;

## output format
## chr fragment-start fragment-end sequenced-allele sequenced-quality read.id target.locus
## for SE data, fragment-end will be recorded as 0
## process data, deal with PE/SE
my ($p1, $s1, $q1, $co1);
my ($p2, $s2, $q2, $co2);
foreach my $sid ( keys %raw ) {
	my $sam = $raw{$sid};

	if( $#$sam == 0 ) {	## SE or only 1 read is in the given region
#		print STDERR "Load $sid as unpaired\n";

		my ($pos, $mseq, $mqual, $cover) = split /\t/, $sam->[0];
		if( $pos <= $t_pos && $pos+length($mseq) > $t_pos ) {	## there is overalp
			if( $PEinfo{$sid} ) {	## PE, its pair is lost which is not good
				print STDERR "WARNING: Fragment $sid overlaps the target while NOT fully loaded! Extend your region!\n";
			} else {	## SE
				next if $mode;
				print join("\t", $t_chr, $pos-1, 0,
						substr($mseq, $t_pos-$pos, 1), substr($mqual, $t_pos-$pos, 1), $sid, $extra ), "\n"; ## 0 means SE data
			}
		}
	} elsif ( $#$sam == 1 ) {	## PE
#		print STDERR "Load $sid as paired\n";

		my @info1 = split /\t/, $sam->[0];
		my @info2 = split /\t/, $sam->[1];
                next if $info1[3] ne $info2[3];
		if( $info1[0] <= $info2[0] ) {
			($p1, $s1, $q1, $co1) = @info1;
			($p2, $s2, $q2, $co2) = @info2;
		} else {
			($p1, $s1, $q1, $co1) = @info2;
			($p2, $s2, $q2, $co2) = @info1;
		}
#		print STDERR "$sid\t$p1\t", $p2+length($s2), "\n";

		next unless $p1<=$t_pos && $p2+length($s2)>$t_pos;	## check overlap
#		print STDERR "Find an overlap\n";

		my ($mseq, $mqual);
		if( $p1 + length($s1) <= $p2 ) {   # no overlap between read1 and read2
			$mseq  = $s1;
			$mseq .= 'N' x ($p2-$p1-length($s1));
			$mseq .= $s2;

			$mqual  = $q1;
			$mqual .= '!' x ($p2-$p1-length($s1));
			$mqual .= $q2;

		} else {	# merge read1 and read2
			$mseq = '';
			$mqual = '';
			if( $p1 + length($s1) < $p2+length($s2) ) {
				my $k;
				$mseq  = substr($s1, 0, $p2-$p1);
				$mqual = substr($q1, 0, $p2-$p1);
				for( $k=$p2-$p1; $k!=length($s1); ++$k ) {   # overlapped region
					if( substr($q1,$k,1) ge substr($q2,$k+$p1-$p2,1) ) {
						$mseq  .= substr($s1, $k, 1);
						$mqual .= substr($q1, $k, 1);
					} else {
						$mseq  .= substr($s2, $k+$p1-$p2, 1);
						$mqual .= substr($q2, $k+$p1-$p2, 1);
					}
				}
				$mseq  .= substr($s2, length($s1)+$p1-$p2);
				$mqual .= substr($q2, length($s1)+$p1-$p2);
			} else {	## R1 completely contains R2, very rare
				$mseq  = $s1;
				$mqual = $q1;
			}
		}
		next if $mode && substr($mseq, $t_pos-$p1, 1) eq 'N';
                my $mm = substr($mseq, $t_pos-$p1, 1);
                next if $mm eq "T" && $co1 =~ /XG:Z:CT/;
                next if $mm eq "A" && $co1 =~ /XG:Z:GA/;
                my ($pp1, $pp2);
                if($co1 =~ /XG:Z:CT/) {
                    $pp1 = $p1;
                    $pp2 = $p2 + 25;
                } elsif($co1 =~ /XG:Z:GA/) {
                    $pp1 = $p1 - 25;
                    $pp2 = $p2;
                }
		print join("\t", $t_chr, $pp1-1, $pp2+length($s2)-1,
				substr($mseq, $t_pos-$p1, 1), substr($mqual, $t_pos-$p1, 1), $sid, $extra ), "\n";
	}
}

