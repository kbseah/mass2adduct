#!/usr/bin/env perl

=head1 NAME

msimunging.pl - Work with large CSV files in mass spectrometry imaging

=head1 SYNOPSIS

msimunging.pl -csv msi.csv -skip 10 -outfmt hashimage -out msi_reformat

msimunging.pl -hash msi.hashimage -filterby freq -freqmin 0.01 -outfmt hashimage -out msi_filtered

msimunging.pl -hash msi.hashimage -filterby top - top 5000 -outfmt hashimage -out msi_filtered

msimunging.pl -help

=head1 DESCRIPTION

Convert between different representations of mass spectrometry imaging data.

Mass spectrometry imaging produces large data sets. The basic data required are
mass peaks (m/z values) and the intensity of each peak at each pixel (spot).
This can be represented in different ways, including:

1. Matrix of intensities arranged as peaks (columns) vs. spots (rows)

2. Table of triplets: peak, spot, intensity

3. Hash (dictionary) of intensity values keyed by peaks (primary key) and spots
(secondary key).

Method 1 is simple to understand, but results in unwieldy files when the numbers
of peaks is large, which is often the case in high-resolution data sets.

Methods 2 and 3 are more indirect, but are more efficient representations when
many mass peaks have zero values for most pixels, i.e. sparse matrices, because
the zeroes are implicit and do not take up disk space.

This script converts between plain-text CSV files (e.g. exported from MSI
software like SCILS or MSIreader) and the triplet and hash formats, to reduce
disk space for efficient processing.

The hash format is a Perl Storable file that can be read directly into a new
Perl hash object with the Storable module.

The triplet format reports five lists in separate files: column numbers, row
numbers, and intensity values for nonzero pixels, and peak masses and spot names
corresponding to the columns and rows respectively. The first three lists can
be used to create a dgTMatrix object in R with the I<sparseMatrix> command from
the CRAN package I<Matrix>.

This script also reports basic parameters such as the dimension of the matrix
and the number of spots.

In addition to format conversion, it also can do simple peak filtering, either
by the proportion of nonzero pixels, or by total intensity.

=cut

use strict;
use warnings;
#use Data::Dumper; # For diagnostics
use Getopt::Long;
use Pod::Usage;
use Storable;
use 5.020;
use POSIX qw (floor);

my $matrix;
my $out = "test";
my $outfmt = 'hashimage';
my $stored;
my $filterby;
my $freqmin = 0.01;
my $toppeaks = 5000;
my $crlf;
my $delim = ",";
my $skip = 0;
my $quiet;

pod2usage(-verbose=>0) if !@ARGV;

GetOptions("csv=s" => \$matrix,
           "hash=s" => \$stored,
           "delim=s" => \$delim,
           "crlf" => \$crlf,
           "skip=i" => \$skip,
           "filterby=s" => \$filterby,
           "freqmin=f" => \$freqmin,
           "top=i" => \$toppeaks,
           "out=s" => \$out,
           "outfmt=s" => \$outfmt,
           "quiet" => \$quiet,
           "help" => sub { pod2usage(-verbose=>1); },
           "man|m" => sub { pod2usage(-verbose=>2); },
           ) or pod2usage(-verbose=>1);

=head1 OPTIONS

=head2 INPUT

Input can be either in CSV or Perl Storable hashimage formats.

=over 8

=item --csv I<FILE>

CSV file of MSI data

=item --delim I<STRING>

Delimiter for CSV file.

Default: "," (comma)

=item --crlf

Input CSV file has CRLF line endings (Windows/DOS-style)

Default: No, CSV file has Unix-style line endings

=item --skip I<INTEGER>

Skip the first N lines of the input CSV file. E.g. if there is header information
or comment fields that are not part of the CSV-delimited input.

Default: 0

=item --hash I<FILE>

Perl Storable hash image produced by this script using the I<--outfmt hashimage>
option.

If this is provided, any input to I<--csv>, I<--delim>, I<--crlf>, I<--skip>
will be ignored.


=back

=head2 FILTERING

=over 8

=item --filterby I<STRING>

What method for filtering? Options: "freq" (frequency of filled pixels), or
"top" (take the top N peaks by total intensity)

Default: None (use all peaks)

=item --freqmin I<NUMERIC>

If option I<--filterby freq> is used. Filter peaks to include only those which
appear in this minimum fraction of all pixels. Should be a number between 0 and 1.

Default: None 

=item --top I<INTEGER>

If option I<--filterby top> is used. Take the top N peaks by total intensity.

Default: 4000

=item 

=back

=head2 OUTPUT

=over 8

=item --outfmt I<STRING>

Output format. Options: "csv", "triplet", "hashimage"

Default: "hashimage"

=item --out I<STRING>

Prefix for output file names

Default: test

=back

=head2 HELP

=over 8

=item --help

Short help message

=item --man

Detailed manual page

=back

=cut

## MAIN ########################################################################

## READ INPUT

my $inhref;
if ($stored) {
    say STDERR "Retrieving previously-hashed matrix from Storable hash image $stored";
    $inhref = retrieve ($stored);
    my @infofields = ('Input file',
                      'Total peaks',
                      'Total spots',
                      'Matrix size',
                      'Filled pixels',
                      'Filled percentage');
    foreach my $key (@infofields) {
        say STDERR join "\t", ($key, $inhref->{'info'}{$key});
    }
} else {
    say STDERR "Reading data from CSV file $matrix";
    $inhref = csv2hash($matrix,$skip,$delim,$crlf);
}

## FILTER PEAK LIST
#print Dumper \%hash;
my $outputref;
if (defined $filterby) {
    if ($filterby eq 'freq' && defined $freqmin) {
        say STDERR "Filtering peaks at $freqmin minimum pixel frequency";
        my $newhref = peakFilterFreq($inhref,$freqmin); # 1% filter
        $outputref = $newhref;
    } elsif ($filterby eq 'top' && defined $toppeaks) {
        say STDERR "Filtering peaks to retain top $toppeaks peaks by total intensity";
        my $newhref = peakFilterTop($inhref,$toppeaks);
        $outputref = $newhref;
    } else {
        say STDERR "Unrecognized option for --filterby, input matrix will not be filtered";
        $outputref = $inhref;
    }
} else {
    say STDERR "Input matrix will not be filtered";
    $outputref = $inhref;
}

## WRITE OUTPUT

if ($outfmt eq 'triplet') {
    say STDERR "Writing files for triplet representation of intensity matrix";
    hash2triplets($outputref,$out);
} elsif ($outfmt eq 'csv') {
    say STDERR "Writing CSV representation of intensity matrix";
    hash2csv($outputref,$out);
} else {
    say STDERR "Writing hash image representation of intensity matrix";
    store $outputref, "$out.hashimage";
}

## SUBS ########################################################################

sub hash2csv {
    my ($href,
        $outprefix) = @_;
    my @peaklist = @{$href->{'peaks'}};
    my @spotlist = @{$href->{'spots'}};
    my $outfile = "$outprefix.csv";
    open(my $fh, ">", $outfile) or die ("$!");
    # Write header line
    say $fh join ",", ('', @peaklist);
    # Report record lines
    my $counter = 0;
    foreach my $spot (@spotlist) {
        my @outline;
        push @outline, $spot;
        foreach my $peak (@peaklist) {
            if (defined $href->{'values'}{$peak}{$spot}) {
                push @outline, $href->{'values'}{$peak}{$spot};
            } else {
                push @outline, 0;
            }
        }
        say $fh join ",", @outline;
        $counter ++;
        say STDERR "... $counter lines written" if $counter%1000==0;
    }
    close($fh);
}

sub hash2triplets {
    my ($href,
        $outprefix) = @_;
    my @rows;
    my @cols;
    my @vals;
    my @peaklist = @{$href->{'peaks'}};
    my @spotlist = @{$href->{'spots'}};
    my %spothash;
    for (my $i=0; $i <= $#spotlist; $i++) {
        $spothash{$spotlist[$i]} = $i;
    }
    for (my $col=0; $col <= $#peaklist; $col++) {
        my $peak = $peaklist[$col];
        foreach my $spot (keys %{$href->{'values'}{$peak}}) {
            my $row = $spothash{$spot};
            push @rows, $row;
            push @cols, $col;
            push @vals, $href->{'values'}{$peak}{$spot};
        }
    }
    array2file(\@rows,"$outprefix\_rows.list");
    array2file(\@cols,"$outprefix\_cols.list");
    array2file(\@vals,"$outprefix\_vals.list");
    array2file(\@peaklist,"$outprefix\_peaknames.list");
    array2file(\@spotlist,"$outprefix\_spotnames.list");
}

sub peakFilterTop {
    my ($href,
        $top) = @_;
    my $peaknum = scalar @{$href->{'peaks'}};
    if ($top > $peaknum) {
        say STDERR "Number of peaks fewer than number to be filtered, not filtering";
        return $href;
    } else {
        say STDERR "Performing peak filtering, retaining $top top peaks by intensity";
        my %peakintensitysum;
        # Sum up intensities per peak
        my $peaksprocessed = 0;
        foreach my $peak (@{$href->{'peaks'}}) {
            $peaksprocessed ++;
            say STDERR "... $peaksprocessed peaks processed" if $peaksprocessed%1000 == 0;
            foreach my $spot (keys %{$href->{'values'}{$peak}}) {
                $peakintensitysum{$peak} += $href->{'values'}{$peak}{$spot};
            }
        }
        # Sort by total intensity per peak, and keep top N peaks
        my @peakstop = sort {$peakintensitysum{$b} <=> $peakintensitysum{$a}} (keys %peakintensitysum);
        @peakstop = @peakstop[0..$top-1];
        # Re-sort peaks by their masses
        @peakstop = sort {$a <=> $b} @peakstop;
        my %newvals = %{$href->{'values'}}{@peakstop};
        my @spots = @{$href->{'spots'}};
        my $dim = $top * scalar (@spots);
        my %newinfo = ("Input file" => $href->{'info'}{'Input file'}." filtered",
                       "Total peaks" => $top,
                       "Total spots" => scalar @spots,
                       "Matrix size" => $dim,
                       #"Filled pixels" => $countfilledpixels,
                       #"Filled percentage" => 100*($countfilledpixels/$dim),
                       );
        my %newhash = ("info" => \%newinfo,
                       "values" => \%newvals,
                       "peaks" => \@peakstop,
                       "spots" => \@spots);
        return \%newhash;
    }
    
}

sub peakFilterFreq {
    # Filter peaks by minimum fraction of spots they must appear in
    # Return new peaklist
    my ($href,
        $minfreq) = @_;
    if ($minfreq > 1 || $minfreq < 0) {
        say STDERR "Error: Peak filtering min fraction must be between 0 and 1";
    } else {
        say STDERR "Performing peak filtering";
        # Find minimum number of spots
        my $minspots = floor($minfreq * scalar @{$href->{'spots'}});
        my @peakfiltered;
        my $countallpeaks = 0;
        my $countfiltpeaks = 0;
        my $countfilledpixels = 0;
        foreach my $key (@{$href->{'peaks'}}) {
            #say STDERR scalar keys %{$href->{$key}};
            $countallpeaks ++;
            say STDERR "... $countallpeaks peaks processed" if $countallpeaks%10000 == 0;
            if (scalar keys %{$href->{'values'}{$key}} > $minspots) {
                push @peakfiltered, $key;
                $countfiltpeaks ++;
                $countfilledpixels += scalar keys %{$href->{'values'}{$key}};
            }
        }
        say STDERR "Total of $countfiltpeaks meet $minfreq criterion";
        my %newvals = %{$href->{'values'}}{@peakfiltered};
        my @spots = @{$href->{'spots'}};
        my $dim = $countfiltpeaks * scalar @spots;
        my %newinfo = ("Input file" => $href->{'info'}{'Input file'}." filtered",
                       "Total peaks" => $countfiltpeaks,
                       "Total spots" => scalar @spots,
                       "Matrix size" => $dim,
                       "Filled pixels" => $countfilledpixels,
                       "Filled percentage" => 100*($countfilledpixels/$dim),
                       );
        my %newhash = ("info" => \%newinfo,
                       "values" => \%newvals,
                       "peaks" => \@peakfiltered,
                       "spots" => \@spots);
        return \%newhash; # Return filtered peaklist
    }
}

sub csv2hash {
    my ($file,  # CSV File
        $skip,  # Initial lines to skip
        $delim, # Delimiter for CSV file
        $crlf,
        ) = @_;
    
    my @peaks;
    my @spots;
    my %hash;
    
    my $mode = $crlf ? "<:crlf" : "<";
    open(my $fh, $mode, $file) or die ("$!");
    # Skip lines if requested
    if ($skip > 0) {
        for (my $linecounter = 1; $linecounter <= $skip; $linecounter ++) {
            my $discard = <$fh>;
        }
    }
    # First line encountered afterwards is header line
    my $header = <$fh>;
    chomp $header;
    @peaks= split /$delim/, $header;
    shift @peaks; # Get rowname off the header line
    
    # Process remaining lines
    my $counter = 0;
    my $row = 0; # 0-based numbering
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ m/^#/; # Skip additional comment lines
        next if $line eq '';
        my @split = split /$delim/, $line;
        my $spot = shift @split;
        push @spots, $spot;
        for (my $col=0; $col <= $#split; $col++) {
            if ($split[$col] > 0) {
                $hash{$peaks[$col]}{$spot} = $split[$col];
                $counter ++;
            }
        }
        $row ++;
        say STDERR "... $row rows processed" if ($row%500 == 0);
    }
    close($fh);
    
    my $dim = scalar(@peaks) * scalar(@spots);
    my $filledpc = 100 * $counter / $dim;
    my %info = ("Input file" => $file,
                "Total peaks" => scalar @peaks,
                "Total spots" => scalar @spots,
                "Matrix size" => $dim,
                "Filled pixels" => $counter,
                "Filled percentage" => $filledpc);
    my %outhash = ("info" => \%info,
                   "peaks" => \@peaks,
                   "spots" => \@spots,
                   "values" => \%hash);
    return \%outhash;
}

sub array2file {
    # Write an array to file, delimited by newlines
    my ($aref, $filename) = @_;
    open(my $fh, ">", $filename) or die ("$!");
    foreach my $val (@$aref) {
        print $fh "$val\n";
    }
    close($fh);
}

