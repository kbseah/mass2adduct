#!/usr/bin/env Rscript
# For an MSI dataset, and a parent ion mass within that dataset
# Find putative adduct peaks that have significant spatial correlation with the
# parent ion, and plot distribution maps for each of the parent and putative
# adduct peaks

# Load optparse package and parse options first before loading other libraries
suppressPackageStartupMessages(library(optparse));

options <- list(
  make_option(
    c("-i", "--imzml"), type="character", action="store",
    help="File name prefix of ImzML format file, without the .imzML prefix. Both the imzML and ibd files should be present."
    ),
  make_option(
    c("-b","--binwidth"), type="numeric", action="store", default="0.001",
    help="Mass resolution for peak binning, in m/z units."
  ),
  make_option(
    c("-t","--triplet"), type="character", action="store",
    help="File name prefix for output files from msimunging.pl script, without the file name suffixes"
  ),
  make_option(
    c("-p","--parentmass"), type="numeric", action="store",
    help="Parent mass for which to find adducts and generate plots"
  ),
  make_option(
    c("-o","--out"), type="character", action="store",
    help="Output plot file name"
  )
)

parser <- OptionParser(
  option_list=options,
  usage="usage: %prog [options]",
  description="Plot MS intensity image for a specified parent mass and putative adducts"
)

# Check that necessary arguments are supplied, otherwise quit

conf <- parse_args(parser, positional_arguments = TRUE);
if (length(conf$options$imzml) != 1) {
  cat ("ERROR: File not specified to option --imzml \n")
  print_help(parser)
  quit(status=2)
} else if (length(conf$options$triplet) != 1) {
  cat ("ERROR: File prefix not specified to option --triplet \n")
  print_help(parser)
  quit(status=2)
} else if (length(conf$options$parentmass) != 1) {
  cat ("ERROR: Parent mass not specified to option --parentmass \n")
  print_help(parser)
  quit(status=2)
}

# Rename input options
ImzMLfile <- conf$options$imzml
binwidth <- conf$options$binwidth # usually 0.001
triplet <- conf$options$triplet # five files that result from msimunging.pl script
parentmass <- conf$options$parentmass
outfile <- conf$options$out

if (length(conf$options$out) != 1) {
  outfile <- paste(as.character(parentmass),
                   "_correlation_image.png",
                  sep = "")
  cat (paste("No output file name provided. Defaulting to ",outfile,"\n",sep=""))
}

# Load libraries
suppressPackageStartupMessages(library(Cardinal))
library(mass2adduct)

# import data
ImzML <- readImzML(ImzMLfile,
                   attach.only = FALSE,
                   as="MSImagingExperiment",
                   resolution = binwidth, # dependent on instrument setup
                   units = "mz") # or ppm
d <- msimat(rows = paste(c(triplet, "_rows.list"), collapse=""),
            cols = paste(c(triplet, "_cols.list"), collapse=""),
            vals = paste(c(triplet, "_vals.list"), collapse=""),
            peaks = paste(c(triplet, "_peaknames.list"), collapse=""),
            spots = paste(c(triplet, "_spotnames.list"), collapse=""))

#-------------------------------------------------------------------------

# mass2adduct analysis
## filter for top 10000 peaks
d.filter <- filterPeaks(d,"topX",10000)
## calculate mass differences
d.diff <- massdiff(d.filter)
## match mass differences to adducts
d.diff.annot <- adductMatch(d.diff,add=adducts2)
## calculate correlation
d.diff.annot.cor <- corrPairsMSI(d,d.diff.annot)
## delete not significant correlations
d.diff.annot.cor.sig <- subset(d.diff.annot.cor,Significance==1)
## order the list based on the mass of the parent ion
d.diff.annot.cor.sig <- d.diff.annot.cor.sig[with(d.diff.annot.cor.sig,order(A)),]
summary(d.diff.annot.cor)

#-------------------------------------------------------------------------

# imaging of correlations
## subset significant correlations for parental mass
d.parentmass <- d.diff.annot.cor.sig[which(d.diff.annot.cor.sig$A == parentmass),]
## plot images with significant correlation to parental mass
options(Cardinal.dark = FALSE) # For dark mode set Cardinal.dark=TRUE and text_col <- "white"
text_col <- "black"
png(filename = outfile,
    width = 3000,
    height = 3000,
    pointsize = 12,
    res = 200)
if (nrow(d.parentmass)>9) {
  par(mfrow = c(4, 4),
      mar = c(5, 5, 5, 5),
      oma = c(0, 0, 5, 0),
      cex = 0.8)
} else {
  par(mfrow = c(3, 3),
      mar = c(5, 5, 5, 5),
      oma = c(0, 0, 5, 0),
      cex = 0.8)
}
print(image(ImzML,
            mz = parentmass,
            plusminus = 0.001,
            colorscale = viridis,
            colorkey = TRUE,
            layout = FALSE))
mtext("M+H",
      col = text_col,
      side = 3,
      line = 1,
      adj = 0,
      cex = 1.5)
for (i in 1:nrow(d.parentmass)) {
  print(image(ImzML,
              mz = d.parentmass$B[i],
              plusminus = 0.001,
              colorscale = viridis,
              colorkey = TRUE,
              layout = FALSE))
  mtext(as.character(d.parentmass$matches[i]),
        col = text_col,
        side = 3,
        line = 1,
        adj = 0,
        cex = 1.5)
}
mtext(paste("Significant correlations of ",parentmass, sep=""),
      col = text_col,
      side = 3,
      line = 0,
      outer = TRUE,
      cex = 2)
dev.off()
