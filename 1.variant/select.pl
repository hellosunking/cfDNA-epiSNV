use strict;
use warnings;

if (@ARGV < 3) {
    die "Usage: $0 <input> <background> <output>\n";
}

my $file_a = $ARGV[0]; 
my $file_b = $ARGV[1];
my $output_file = $ARGV[2];

open(my $fh_b, '<', $file_b) or die "Cannot open file $file_b: $!";
my $header_line = <$fh_b>;
chomp $header_line;
my @headers = split(/\t/, $header_line);

my %samples_to_remove;
for my $i (1..$#headers) {
    $samples_to_remove{$headers[$i]} = 1;
}
close($fh_b);

open(my $fh_a, '<', $file_a) or die "Cannot open file $file_a: $!";
open(my $fh_out, '>', $output_file) or die "Cannot open file $output_file: $!";

my $header_line_a = <$fh_a>;
chomp $header_line_a;
print $fh_out "$header_line_a\n";

while (my $line = <$fh_a>) {
    chomp $line;
    my @fields = split(/\t/, $line);
    
    if (exists $samples_to_remove{$fields[0]}) {
        next;
    }
    
    print $fh_out "$line\n";
}

close($fh_a);
close($fh_out);
