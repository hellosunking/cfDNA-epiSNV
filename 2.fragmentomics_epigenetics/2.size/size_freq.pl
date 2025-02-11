#!/usr/bin/perl

use strict;
use warnings;

my $file = $ARGV[0] or die "Usage: $0 <file_path>\n";

open(my $fh, '<', $file) or die "Cannot open $file: $!";

my $sum_filtered = 0;
my $sum_total = 0;

while (my $line = <$fh>) {
    chomp $line;
    my @columns = split("\t", $line);
    if ($columns[0] <= 150) {
        $sum_filtered += $columns[1];
    }
    $sum_total += $columns[1];
}

close($fh);

my $result = $sum_filtered / $sum_total;

print "$result\n";
