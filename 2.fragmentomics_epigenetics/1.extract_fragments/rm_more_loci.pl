#!/usr/bin/perl
use strict;
use warnings;

my @lines;
while (my $line = <STDIN>) {
    chomp $line;
    push @lines, $line;
}

my %id_count;
foreach my $line (@lines) {
    my @fields = split(/\t/, $line);
    $id_count{$fields[5]}++;
}

my %unique_ids;
foreach my $id (keys %id_count) {
    if ($id_count{$id} == 1) {
        $unique_ids{$id} = 1;
    }
}

foreach my $line (@lines) {
    my @fields = split(/\t/, $line);
    if (exists $unique_ids{$fields[5]}) {
        print "$line\n"; 
    }
}
