#!/usr/bin/perl
use strict;
use warnings;

my $input_file = $ARGV[0];
my %id_count;

open(my $fh, '<', $input_file) or die "Cannot open file: $!";
while (my $line = <$fh>) {
    chomp $line;
    my @fields = split(/\t/, $line);
    $id_count{$fields[5]}++;
}
close($fh);

my %unique_ids;
foreach my $id (keys %id_count) {
    if ($id_count{$id} == 1) {
        $unique_ids{$id} = 1;
    }
}

open($fh, '<', $input_file) or die "Cannot open file: $!";
while (my $line = <$fh>) {
    chomp $line;
    my @fields = split(/\t/, $line);
    if (exists $unique_ids{$fields[5]}) {
        print "$line\n"; 
    }
}
close($fh);
