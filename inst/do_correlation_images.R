#!/usr/bin/env Rscript
# For an MSI dataset, and a parent ion mass within that dataset
# Find putative adduct peaks that have significant spatial correlation with the
# parent ion, and plot distribution maps for each of the parent and putative
# adduct peaks

args = commandArgs(trailingOnly=TRUE)

# Get input filename prefix
ImzMLfile <- args[1]
binwidth <- as.numeric(args[2]) # usually 0.001
triplet <- args[3] # five files that result from msimunging.pl script
parentmass <- as.numeric(args[4])
outfile <- args[5]

# load packages
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
png(filename = paste(as.character(parentmass), # Experimental design for file name
                     "_correlation_image.png",
                     sep = ""),
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
