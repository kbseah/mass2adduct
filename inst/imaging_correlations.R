# imaging correlations based on mass2adduct results

## load packages
library(Cardinal) # Analysis tool
library(mass2adduct)
library(plot3D)   # for external colorkey of MS images
library(viridis)  # color palette

## read data
SDHBp_file_spots <- "MetaMix_SDHBp_60_900_pos_A28_100um_89x74" # file name of *.imzML and *.ibd, both need to be present
folder           <- "/msi_datasets" # path of msi data
SDHBp_spots <- readImzML(SDHBp_file_spots, 
                         folder,
                         attach.only = FALSE,
                         #as="MSImagingExperiment",
                         resolution = 0.005,
                         units = "mz")
d <- msimat(rows = "intensity_matrices/MetaMix_SDHBp_munged_rows.list",
            cols = "intensity_matrices/MetaMix_SDHBp_munged_cols.list",
            vals = "intensity_matrices/MetaMix_SDHBp_munged_vals.list",
            peaks = "intensity_matrices/MetaMix_SDHBp_munged_peaknames.list",
            spots = "intensity_matrices/MetaMix_SDHBp_munged_spotnames.list")

# filter for top 10000 peaks
d.filter <- filterPeaks(d,
                        "topX",
                        10000)

# calculate mass differences
d.diff <- massdiff(d.filter)

# match mass differences to adducts
d.diff.annot <- adductMatch(d.diff,
                            add = adducts2) # for extended adduct list use "add = adducts"

# calculate correlation
d.diff.annot.cor <- corrPairsMSI(d,
                                 d.diff.annot)

# delete not significant correlations
d.diff.annot.cor <- subset(d.diff.annot.cor,
                           Significance == 1)

# delete all matches for different matrices
d.diff.annot.cor <- subset(d.diff.annot.cor,
                           matches != "CHCA") # DHB was used, CHCA correlations not needed

# order the list based on the mass of the parent ion
d.diff.annot.cor <- d.diff.annot.cor[with(d.diff.annot.cor,
                                          order(A)),]
summary(d.diff.annot.cor)

# imaging of correlations
png(filename = "MetaMix_SDHBp_corr_Carnitine.png",
    width = 2000,
    height = 2000,
    pointsize = 12,
    bg = "white",
    res = 200)
par(mfrow = c(2,2))
image(SDHBp_spots,
      mz = 162.112,
      plusminus = 0.001,
      col.regions = viridis(n = 100),
      colorkey = TRUE)
mtext("M+H",
      col = "black",
      side = 3,
      line = 1,
      adj = 0,
      cex = 1.5)
image(SDHBp_spots,
      mz = 184.094,
      plusminus = 0.001,
      col.regions = viridis(n = 100),
      colorkey = TRUE)
mtext("M+Na",
      col = "black",
      side = 3,
      line = 1,
      adj = 0,
      cex = 1.5)
image(SDHBp_spots,
      mz = 200.068,
      plusminus = 0.001,
      col.regions = viridis(n = 100),
      colorkey = TRUE)
mtext("M+K",
      col = "black",
      side = 3,
      line = 1,
      adj = 0,
      cex = 1.5)
image(SDHBp_spots,
      mz = 298.128,
      plusminus = 0.001,
      col.regions = viridis(n = 100),
      colorkey = TRUE)
mtext("M+DHB",
      col = "black",
      side = 3,
      line = 1,
      adj = 0,
      cex = 1.5)
dev.off()