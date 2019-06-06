# Peak picking with Cardinal

## Load package
library(Cardinal)

## import imzML data
file <- "msi_datasets/MetaMix_SDHBp_60_900_pos_A28_100um_89x74" # file name of *.imzML and *.ibd, both need to be present
msiset <- readImzML(file,
                   attach.only = FALSE,
                   as="MSImageSet",
                   resolution = 0.0005, # 0.0005 results in a bin width of 0.001 Da
                   units = "mz")

## apply peakpick function
peaks <- peakPick(msiset,
                  method = "adaptive",
                  SNR = 10)

## creat list of all peaks in all pixels
peaks_list <- pData(mzData(imageData(peaks)))

## unlist peaks
peaks_unlist <- unlist(peaks_list)
head(peaks_unlist)

## keep only unique peaks
peaks_unique <- unique(peaks_unlist)

## sort by increasing m/z
peaks_unique_sorted <- sort(peaks_unique,
                            decreasing = FALSE)

## save as text file
write.table(peaks_unique_sorted,
            file = "peaklist.txt",
            row.names = FALSE,
            col.names = FALSE,
            sep = ",")
