#!/usr/bin/perl
use strict;
use warnings;

my $line1 = <STDIN>;
chomp $line1;
my @fields1 = split /\t/, $line1;
my $chr_raw = $fields1[0];
my $pos_raw = $fields1[1];

my $pre_line = '';

while (my $line2 = <STDIN>) {
    chomp $line2;
    my @fields2 = split /\t/, $line2;
    my $chr_new = $fields2[0];
    my $pos_new = $fields2[1];

    if ($.==2) {
        if ($chr_new eq $chr_raw && $pos_new == $pos_raw +1) {
            print "$line1\n$line2\n";
        }
            $chr_raw = $chr_new;
            $pos_raw = $pos_new;
            $pre_line = $line2;            
    }else{            
        if ($chr_new eq $chr_raw && $pos_new == $pos_raw +1) {
                print "$pre_line\n$line2\n";
        }
        $chr_raw = $chr_new;
        $pos_raw = $pos_new;
        $pre_line = $line2;
    }
} 
