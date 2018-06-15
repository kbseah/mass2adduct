# mass2adduct - Finding molecular adducts in mass spectrometry data

In mass spectrometry imaging, adducts can form between target molecules (e.g. metabolites) and other substances such as matrix or salt ions. This package presents tools for counting and identifying possible adducts in MS data, and accompanies Janda et al. (in prep.).

## Install package

You will need the [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package, if you don't have it already:

```R
install.packages("devtools") # Install from CRAN
library(devtools) # Load the package
```

Install mass2adduct from Github:

```R
install_github("kbseah/mass2adduct")
```

## Import data

### Import MS imaging intensity data

MSI data exported from the MSiReader software with the "intensity export" function. Import the data into R as a data.frame:

```R
dat.msi <- readMSI("msi.csv", sep=";", remove.zeroes=TRUE)
```

### Import list of masses directly

If you do not intend to use pixel intensity values for downstream analyses, you can import a list of mass values as a plain text file (one per line):

```R
dat <- scan("mass_example",what=numeric())
```

## Tabulate all pairwise mass differences

Take all possible pairs of masses and calculate the mass difference for each pair. These mass differences represent potential molecular adducts. 

```R
dat.msi.diff <- diffTabulateMSI(d=dat.msi,override.limit=TRUE) # From MSI data
dat.diff <- diffTabulate(d=dat) # From simple list of mass values
```

Output is a data.frame with three columns: the two parent masses and the difference between them.

## Bin mass differences into a histogram

The calculated mass differences are misleadingly precise, because measurement error and uncertainty are not taken into account. They should be binned into a histogram with a user-specified bin width, that depends on the known mass precision of your instrument. In this example, we use the value of 0.01:

```R
dat.diffhist <- diffHist (d = dat.diff, width = 0.01, plot = TRUE) # Produces a histogram plot
```

The resulting object `dat.diffhist` is a standard R histogram object. You can "zoom" into specific regions of the plot with the `xlim` parameter in the R `plot` command.

```R
plot (dat.diffhist, xlim=c(100,150))
```

## Count occurrences of known adducts

The package comes with a small data set of common molecular adducts, called `adducts` (use the command `help(adducts)` to see a description). Users can supply their own custom sets of adducts as long as they are in the same format.

The following function looks for known adducts by finding the closest-matching bin in the mass difference histogram produced above. It reports the number of counts (i.e. how many pairs of MS peaks have that mass difference) and the quantile.

Note that quantile values will usually be quite high because the majority of mass differences have zero to few counts.

```R
adductMatch (dat.diffhist, n=10) # Show the top ten matches to known adducts
```

## List the highest-counted mass differences and any matches to known adducts

The procedure above (`adductMatch`) only looks for known adducts. What about peaks in the histogram which may not have a good match with the reference list? The function `topAdducts` performs the complementary function: Rank mass differences by the number of times they are observed, and report any matches to known adducts.

```R
topAdducts(dat.diffhist,n=10) # Report top ten hits
```

## Find mass peaks associated with specific mass difference

Find parent and derivative mass peaks associated with a specific molecular adduct, e.g. DHB-H2O (mass 136.01600)

```R
dat.DHBH2O.peaks <- diffGetPeaks (dat.diff, mass=136.01600, width=0.01)
```

Output is a data.frame with three columns, each row representing a putative pair of parent and derivative ion masses.

## Test for significant correlations between mass peaks

Test for spatial correlations between mass peaks in MS imaging data (imported with the `readMSI` function).

### Parent ions vs putative derivative ions

For example, we wish to see if what we believe to be pairs of parent- and derivative-ion masses tend to occur together. If they are truly related by molecular adduct formation, then their abundances should be correlated. These putative pairs are found with the `diffGetPeaks` function (above).

```R
dat.DHBH2O.corr <- corrPairsMSI (d=dat.msi, pairs = dat.DHBH2O.peaks)
```

Output is a data.frame with p-values for each ion pair. By default the cutoff for significance is p=0.05 with Bonferroni correction.

### Target mass of interest vs all masses

Compare only a single target mass value against all masses in the MSI data

```R
dat.350.988.corr <- corrSinglePeakMSI(d=dat.msi, peak = 350.988)
```

## Reformatting large CSV files

Plain-text CSV files of MSI data exported by software such as SCILS or MSIreader can be large, on the order of several Gb. For many MSI data sets, a lot of this is "wasted" because the majority of peaks are only detected in a minority of pixels. Low-abundance peaks can be thrown out through filtering, but the majority of entries in the data matrix may still be zeroes.

Reading these large files directly into R is inefficient and may not be possible because of insufficient memory. To be more economical, the data can be represented in a "triplet" format rather than matrix.

A Perl script `msimunging.pl` (in the `inst/` subdirectory of the package source) is provided to do this format conversion before importing the files into R. It can convert between plain CSV format, the so-called "triplet" format, and a binary format readable by Perl (mostly used for testing). Note that it is nonetheless still limited by the amount of memory you have on your system (e.g. if your CSV file is 40 Gb and your RAM is 32 Gb, then it will probably not work).

Conversion to triplet format output gives five files, which correspond to the inputs required by the `readMSI()` function in the R package.

The script also provides the option to filter peaks by their frequency of nonzero pixels.

Usage information is shown by typing `msimunging.pl --help`.


## Help and documentation

Further documentation for all the above functions can be found with the `help` or `?` functions in R. Overview of the entire package can be viewed with `help(mass2adduct)`.

Please report any problems to the package maintainer, either by email or with the issue tracker on Github.

## Contact

Package maintainer: Brandon Seah (kbseah@mpi-bremen.de)

Authors: Moritz Janda, Manuel Liebeke

[Department of Symbiosis, Max Planck Institute for Marine Microbiology](https://www.mpi-bremen.de/en/Department-of-Symbiosis.html)


