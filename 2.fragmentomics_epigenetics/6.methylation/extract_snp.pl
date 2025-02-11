#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV < 2) {
    die "Usage: $0 <input> <output>\n";
}

my $input_file = $ARGV[0];
my $output_file = $ARGV[1];

my %n_data;
my %t_data;

open(my $fh, '<', $input_file) or die "Cannot open file $input_file: $!";

while (my $line = <$fh>) {
    chomp $line;
    my @fields = split(/\t/, $line);
    
    my @info = split(/,/, $fields[6]);

    my $joined_str = join("_", $fields[0], @info);

    if ($fields[3] eq $info[1]) {
        $n_data{$joined_str} = 1;
    }

    if ($fields[3] eq $info[2]) {
        $t_data{$joined_str} = 1;
    }
}

close($fh);

my %final_data;
foreach my $key (keys %n_data) {
    if (exists $t_data{$key}) {
        $final_data{$key} = 1;
    }
}

open(my $out_fh, '>', $output_file) or die "Cannot open file $output_file: $!";

foreach my $key (keys %final_data) {
    print $out_fh "$key\n";
}

close($out_fh);
