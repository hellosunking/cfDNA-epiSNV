#!/usr/bin/perl
use strict;
use warnings;

my ($fq1, $fq2, $prefix) = @ARGV;
die "Usage: $0 <fq1.gz> <fq2.gz> <prefix>\n" unless @ARGV == 3;

my $id_file = "$prefix.id";
open my $id_fh, '<', $id_file or die "Cannot open $id_file: $!\n";
my %id_hash;
while (my $line = <$id_fh>) {
    chomp $line;
    my @col = split(/\s+/, $line);
    my $name = '@' . $col[0];
    $id_hash{$name} = 1;
}
close $id_fh;

open my $fq1_fh, '-|', "zcat $fq1" or die "Cannot open $fq1: $!\n";
open my $fq2_fh, '-|', "zcat $fq2" or die "Cannot open $fq2: $!\n";

open my $out_fq1, '|-', "gzip > $prefix.R1.fq.gz" or die "Cannot open output file for R1: $!\n";
open my $out_fq2, '|-', "gzip > $prefix.R2.fq.gz" or die "Cannot open output file for R2: $!\n";

while (1){
        my $id_a = <$fq1_fh>; 
        my $id_b = <$fq2_fh>;
        last unless ($id_a and $id_b);
        my $seq_a = <$fq1_fh>; my $str_a = <$fq1_fh>; my $q_a = <$fq1_fh>;
        my $seq_b = <$fq2_fh>; my $str_b = <$fq2_fh>; my $q_b = <$fq2_fh>;

        chomp $id_a; chomp $id_b;
	my @a = split(/\s+/, $id_a);
	my @b = split(/\s+/, $id_b);
	next if $a[0] ne $b[0];

        if(exists $id_hash{$a[0]}){
		print $out_fq1 "$id_a\n", $seq_a, $str_a, $q_a;
		print $out_fq2 "$id_b\n", $seq_b, $str_b, $q_b;
	}
}

close $fq1_fh;
close $fq2_fh;
close $out_fq1;
close $out_fq2;
