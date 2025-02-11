#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV < 3) {
    die "Usage: $0 <input> <path> <output_prefix>\n";
}

my $input = $ARGV[0];
my $path = $ARGV[1];
my $output_prefix = $ARGV[2];

my %regions_promoter;
my %regions_genebody;

open my $promoter_fh, '<', "$path/promoter.region" or die "Cannot open $path/promoter.region: $!";
while (<$promoter_fh>) {
    chomp;
    my ($chr, $start, $end) = split;
    push @{ $regions_promoter{$chr} }, { start => $start, end => $end };
}
close $promoter_fh;

open my $genebody_fh, '<', "$path/genebody.region" or die "Cannot open $path/genebody.region: $!";
while (<$genebody_fh>) {
    chomp;
    my ($chr, $start, $end) = split;
    push @{ $regions_genebody{$chr} }, { start => $start, end => $end };
}
close $genebody_fh;

open my $input_fh, '-|', "zcat $input" or die "Cannot open $input: $!";

open my $promoter_out, '|-', "gzip > $output_prefix.promoter.CpG.meth.call.gz" or die "Cannot open $output_prefix.promoter.CpG.meth.call.gz for writing: $!";
open my $genebody_out, '|-', "gzip > $output_prefix.genebody.CpG.meth.call.gz" or die "Cannot open $output_prefix.genebody.CpG.meth.call.gz for writing: $!";
open my $inter_out, '|-', "gzip > $output_prefix.intergenic.CpG.meth.call.gz" or die "Cannot open $output_prefix.intergenic.CpG.meth.call.gz for writing: $!";

my $header = <$input_fh>;
print $promoter_out $header;
print $genebody_out $header;
print $inter_out $header;

while (<$input_fh>) {
    chomp;
    my ($chr, $pos, @other_info) = split;

    my $printed = 0;

    if (exists $regions_promoter{$chr}) {
        foreach my $region (@{ $regions_promoter{$chr} }) {
            if ($pos >= $region->{start} && $pos <= $region->{end}) {
                print $promoter_out "$_\n";
                $printed = 1;
                last;
            }
        }
    }

    if (exists $regions_genebody{$chr} && !$printed) {
        foreach my $region (@{ $regions_genebody{$chr} }) {
            if ($pos >= $region->{start} && $pos <= $region->{end}) {
                print $genebody_out "$_\n";
                $printed = 1;
                last;
            }
        }
    }

    if (!$printed) {
        print $inter_out "$_\n";
    }
}

close $input_fh;
close $promoter_out;
close $genebody_out;
close $inter_out;
