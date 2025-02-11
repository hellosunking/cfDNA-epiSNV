use strict;
use warnings;

my $prefix = $ARGV[0];

my @files = glob("$prefix.test.rep*.txt");

my %data;

foreach my $file (@files) {
    open my $fh, '<', $file or die "can not open $file: $!";
    while (<$fh>) {
        chomp;
        my @fields = split /\t/;
        my $sample = $fields[0];
        my $value  = $fields[2];
        push @{ $data{$sample} }, $value;
    }

    close $fh;
}

foreach my $sample (sort keys %data) {
    my $sum   = 0;
    my $count = scalar @{ $data{$sample} };
    $sum += $_ for @{ $data{$sample} };
    my $mean = $sum / $count;
    print "$sample\t$mean\n";
}

