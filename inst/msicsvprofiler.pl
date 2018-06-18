#!/usr/bin/env perl

=head1 NAME

msicsvprofiler.pl - Profile intensity and filled pixels per peak for large MSI CSV files

=head1 SYNOPSIS

msicsvprofiler.pl -csv msi.csv -skip 10 -delim="," -limit 100 -out intensitystats.tsv

msicsvprofiler.pl -help

=head1 DESCRIPTION

Parse a CSV file containing mass spectrometry imaging (MSI) data, where peaks
are columns and rows are spots. Column headers contain peak m/z values, and row
names are present.

The output is a tab-separated table with three columns: I<peaks> (m/z values per
mass peak), I<counts> (number of nonzero pixels per peak), and I<intensities>
(total mass intensities across all pixels per peak).

The result can be analyzed to decide on a peak filtering cutoff to be applied,
e.g. with the I<analyzeIntensityCutoffsCumul> function in the R package
I<mass2adduct>.

=cut

use strict;
use warnings;
use 5.020;
use Getopt::Long;
use Storable;
use Data::Dumper;
use Pod::Usage;

my $csv;
my $out = "test";

my $crlf;
my $delim=",";
my $skip = 0;
my $limit = -1;

my %hash;
my @peaks;

pod2usage(-verbose=>0) if !@ARGV;

GetOptions("csv=s" => \$csv,
           "out=s" => \$out,
           "crlf" => \$crlf,
           "delim=s" => \$delim,
           "skip=i" => \$skip,
           "limit=i" => \$limit,
           "help" => sub { pod2usage(-verbose=>1); },
           "man" => sub { pod2usage(-verbose=>2); },
           ) or pod2usage(-verbose=>0);

=head1 OPTIONS

=over 8

=item --csv I<FILE>

Name of CSV file to process

=item --out I<STRING>

Prefix for output file name

Default: "test"

=item --crlf

Input CSV file uses CRLF line endings (Windows/DOS-style)

Default: no

=item --delim I<STRING>

Record delimiter on each line of the input CSV file

Default: "," (comma-separated variables)

=item --skip I<INT>

Skip the first N lines of the file (e.g. headers, comments)

Default: 0

=item --limit I<INT>

Stop processing after N rows

Default: -1 (no limit)

=item --help

Short help message

=item --man

Full manual page

=back

=cut

# Open file
my $mode = $crlf ? "<:crlf" : "<"; # Check line ending type
open(my $fh, $mode, $csv) or die ("Cannot open file $csv: $!");
if ($skip > 0) { # Skip lines
    for (my $i=0; $i < $skip; $i++) {
        my $discard = <$fh>;
    }
}
# First line is header with peak names
my $header = <$fh>;
chomp $header;
@peaks = split /$delim/, $header;
shift @peaks; # First value is row label
say STDERR "Total peaks: ".scalar(@peaks);
# Start processing entries
my $counter = 0;
while (my $line = <$fh>) {
    last if ($limit > 0 && $counter > $limit);
    $counter ++;
    say STDERR "... $counter lines processed" if $counter%1000 == 0;
    chomp $line;
    my @split = split /$delim/, $line;
    shift @split;
    if (scalar @split != scalar @peaks) {
        say STDERR "Error: Number of entries doesn't match header in entry row $counter";
        say STDERR "Stop processing further entries";
        last;
    } else {
        for (my $i=0; $i <= $#peaks; $i++) {
            if ($split[$i] > 0) {
                $hash{$peaks[$i]}{'count'} ++;
                $hash{$peaks[$i]}{'totalint'} += $split[$i];
            }
        }
    }
}
close($fh);
say STDERR "Total $counter entries processed";

# Report intensity per peak profile
say STDERR "Reporting output to file $out.tsv";
open(my $fhout, ">", "$out.tsv") or die ("$!");
# Print header
say $fhout "# Input file $csv";
say $fhout "# Total peaks ".scalar (@peaks);
say $fhout "# Total lines $counter";
say $fhout join "\t", qw(peaks counts intensities);
foreach my $peak (@peaks) {
    my @out = ($peak);
    # Account for possible all-zeroes peaks (e.g. from subsetted data)
    push @out, $hash{$peak}{'count'} ? $hash{$peak}{'count'} : 0;
    push @out, $hash{$peak}{'totalint'} ? $hash{$peak}{'totalint'} : 0;
    say $fhout join "\t", @out;
}
close($fhout);
