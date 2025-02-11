#!/usr/bin/perl
use strict;
use warnings;

my $file = $ARGV[0];

open my $fh, '<', $file;
my $header_line = <$fh>;
chomp $header_line;
my @header = split /\t/, $header_line;

my @data;
while (<$fh>) {
    chomp;
    my @row = split /\t/;
    push @data, \@row;
}
close $fh;

print join("\t", @header), "\n";
for my $row (0..$#data) {
    print $data[$row][0];
    for my $col (1..$#header) {
        my $value = $data[$row][$col] / sum_col($col, \@data);
        print "\t$value";
    }
print "\n";
}

sub sum_col {
    my ($col, $data_ref) = @_;
    my $sum = 0;
    foreach my $row (0..$#$data_ref) {
        $sum += $data_ref->[$row][$col];
    }
    return $sum;
}
