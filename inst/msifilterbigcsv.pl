#!/usr/bin/env perl

=head1 NAME

msifilterbigcsv.pl - Subset large CSV file into hash image

=head1 SYNOPSIS

msifilterbigcsv.pl --csv input.csv --peaks peaklist --out output_prefix

msifilterbigcsv.pl --help

=head1 DESCRIPTION

Convert large CSV file into hash image, when file is too large to read into
memory even as a sparse matrix representation with I<msimunging.pl>.

Supply a shortlist of peaks to extract from the original CSV file, e.g. after
looking at the intensity distribution with I<msicsvprofiler.pl>. The original
data will then be filtered using this shortlist, and reported as a hashimage
which can be converted to CSV or triplet format with I<msimunging.pl>

=cut

use strict;
use warnings;
use 5.020;
use Storable;
use Pod::Usage;
use Getopt::Long;

my $csv;
my $peakshortlist;
my $delim = ",";
my $crlf;
my $skip = 0;
my $out = "test";

pod2usage(-verbose=>0) if !@ARGV;

GetOptions("csv=s"=>\$csv,
           "peaks=s"=>\$peakshortlist,
           "delim=s"=>\$delim,
           "crlf"=>\$crlf,
           "skip=i"=>\$skip,
           "out=s"=>\$out,
           "help"=> sub{ pod2usage(-verbose=>1); },
           "man"=> sub{ pod2usage(-verbose=>2); },
           ) or pod2usage(-verbose=>1);

=head1 OPTIONS

=over 8

=item --csv I<FILE>

CSV file

=item --peaks I<FILE>

List of peaks to extract from the CSV file

=item --delim I<STRING>

Delimiter used in CSV file.

Default: "," (comma-separated)

=item --crlf

CSV file uses CRLF line endings

Default: no

=item --skip I<INTEGER>

Number of lines (e.g. header, comments) in CSV file to skiop

Default: 0

=item --out I<STRING>

Prefix for output file name

Default: "test"

=item --help

Short help message

=item --man

Full manual page

=back

=cut

## MAIN ########################################################################

my $peakshref = read_peaklist($peakshortlist);
say STDERR "Processing CSV file";
my $filtered = process_csv($csv, $peakshref, $delim, $crlf, $skip);

say STDERR "Writing output to file $out.hashimage";
store $filtered, "$out.hashimage";

## SUBS ########################################################################

sub process_csv {
    my ($file, $peakshref, $delim, $crlf, $skip) = @_;
    my %hash; # Hash of 
    my @peakidx;
    my @peaknames;
    my %idx2peak;
    my @spots; 
    my $mode = $crlf? "<:crlf" : "<";
    open(my $fh, $mode, $file) or die ("$!");
    # Skip lines
    if ($skip > 0) {
        for (my $i=1; $i <= $skip; $i++) {
            my $discard = <$fh>;
        }
    }
    # Header line
    my $head = <$fh>;
    chomp $head;
    my @headsplit = split /$delim/, $head;
    shift @headsplit; # Discard rowlabel
    for (my $i=0; $i <= $#headsplit; $i++) {
        if (defined $peakshref->{$headsplit[$i]}) {
            push @peakidx, $i;
            push @peaknames, $headsplit[$i];
            $idx2peak{$i} = $headsplit[$i];
        }
    }
    say STDERR scalar(@peaknames)." peak names matched in CSV file header";
    # Process entries
    my $countfilledpixels=0;
    my $linecount=0;
    while (my $line = <$fh>) {
        $linecount ++;
        say STDERR "... $linecount lines processed" if $linecount%1000 == 0;
        chomp $line;
        my @split = split /$delim/, $line;
        my $spotname = shift @split;
        push @spots, $spotname;
        foreach my $idx (@peakidx) {
            if ($split[$idx] > 0) {
                $hash{$idx2peak{$idx}}{$spotname} = $split[$idx];
                $countfilledpixels ++;
            }
        }
    }
    close($fh);
    my $dim = scalar(@peaknames)*scalar(@spots);
    my %info = ("Input file" => $file,
               "Total peaks" => scalar @peaknames,
               "Total spots" => scalar @spots,
               "Matrix size" => $dim,
               "Filled pixels" => $countfilledpixels,
               "Filled percentage" => 100*($countfilledpixels/$dim),
               );
    my %out = ("info" => \%info,
               "values" => \%hash,
               "peaks" => \@peaknames,
               "spots" => \@spots);
    return (\%out);
}

sub read_peaklist {
    my ($peaklist) = @_;
    my %hash;
    open(my $fh, "<", $peaklist) or die ("$!");
    while (my $line = <$fh>) {
        chomp $line;
        $hash{$line} ++;
    }
    close($fh);
    say STDERR "Read shortlist of ".scalar(keys %hash)." peaks";
    return \%hash;
}