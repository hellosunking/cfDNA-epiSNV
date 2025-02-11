#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV < 2) {
    die "Usage: $0 <Sid> <region>\n";
}

my $sid = $ARGV[0];
my $region = $ARGV[1];

my $n_file = "$sid.n.$region.CpG.meth.call.gz";
my $t_file = "$sid.t.$region.CpG.meth.call.gz";

my ($a, $b, $c, $d);
open my $n_fh, '-|', "zcat $n_file" or die "Cannot open $n_file: $!\n";
while (<$n_fh>) {
    next if $. == 1;
    my @fields = split;
    $a += $fields[3];
    $b += $fields[4];
    $c += $fields[7];
    $d += $fields[8];
}
close $n_fh;

my ($a_t, $b_t, $c_t, $d_t);
open my $t_fh, '-|', "zcat $t_file" or die "Cannot open $t_file: $!\n";
while (<$t_fh>) {
    next if $. == 1;
    my @fields = split;
    $a_t += $fields[3];
    $b_t += $fields[4];
    $c_t += $fields[7];
    $d_t += $fields[8];
}
close $t_fh;

my $wt = ($a + $c)/($a + $b + $c +$d);
my $mut = ($a_t + $c_t)/($a_t + $b_t + $c_t +$d_t);
my $diff = $mut - $wt;

print "$sid\t$wt\t$mut\t$diff\n";
