#!/usr/bin/perl
use strict;
use warnings;

my $sample = $ARGV[0];
die "Usage: $0 <sample_name>\n" unless defined $sample;

# Define file paths
my $alt_d = "$sample.Mut.all.D.dist";
my $alt_u = "$sample.Mut.all.U.dist";
my $ref_d = "$sample.Wt.all.D.dist";
my $ref_u = "$sample.Wt.all.U.dist";

# Function to calculate the sum of values for a specific range in a dist file
sub calculate_sum {
    my ($file) = @_;
    my $sum = 0;
    
    open my $fh, '<', $file or die "Could not open '$file': $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        my @fields = split /\t/, $line;
        if ($fields[0] >= -50 && $fields[0] <= 50) {
            $sum += $fields[1];
        }
    }
    close $fh;
    
    return $sum;
}

# Function to calculate the sum of the second column in a dist file
sub calculate_total {
    my ($file) = @_;
    my $sum = 0;
    
    open my $fh, '<', $file or die "Could not open '$file': $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        my @fields = split /\t/, $line;
        $sum += $fields[1];
    }
    close $fh;
    
    return $sum;
}

# Calculate the values for control group (CTR)
my $a = calculate_sum($alt_d);
my $b = calculate_total($alt_d);
my $c = calculate_sum($alt_u);
my $d = calculate_total($alt_u);
my $a1 = calculate_sum($ref_d);
my $b1 = calculate_total($ref_d);
my $c1 = calculate_sum($ref_u);
my $d1 = calculate_total($ref_u);

my $wt = ($a1+$c1)/($b1+$d1);
my $mut = ($a+$c)/($b+$d);
my $diff = $mut - $wt;

# Calculate and print final result (Alt/Ref ratios for both groups)
print join("\t", $sample, $wt, $mut, $diff) . "\n";
