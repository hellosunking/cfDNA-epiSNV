#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV < 3) {
    die "Usage: $0 <input> <id> <output_prefix>\n";
}

my $input_file = $ARGV[0];
my $id = $ARGV[1];
my $output = $ARGV[2];

my %name;
open my $IN, '<', $id or die("Cannot open $id: $!");
while( <$IN> ) {
    chomp;
    my @l = split /_/; 
    my $snp = join("\t", @l); 
    $name{$snp} = 1; 
}
close $IN;

open(my $fh, '<', $input_file) or die "Cannot open file $input_file: $!";
open(my $out_t, '>', "$output.mut.id") or die "Cannot open file $output.mut.id: $!";
open(my $out_n, '>', "$output.wt.id") or die "Cannot open file $output.wt.id: $!";

while (my $line = <$fh>) {
    chomp $line;
    my @fields = split(/\t/, $line);

    my @info = split(/,/, $fields[6]);

    my $snpnew = join("\t", $fields[0], @info);

    if (exists $name{$snpnew}) {
        if ($fields[3] eq $info[1]) {
            print $out_n "$fields[5]\n";
        }
        elsif ($fields[3] eq $info[2]) {
            print $out_t "$fields[5]\n";
        }
    }
}

close($fh);
close($out_t);
close($out_n);
