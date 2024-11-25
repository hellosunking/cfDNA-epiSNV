use strict;
use warnings;

my ($input_file, $output_file) = @ARGV;

open my $fh, '<', $input_file;

my %data;

while (<$fh>) {
    chomp;
    my ($sample, $value) = split;
    push @{ $data{$sample} }, $value;
}

close $fh;

open my $out, '>', $output_file;

foreach my $sample (sort keys %data) {
    my $sum   = 0;
    my $count = scalar @{ $data{$sample} };
    $sum += $_ for @{ $data{$sample} };
    my $mean = $sum / $count;
    print $out "$sample\t$mean\n";
}

close $out;

