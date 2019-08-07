#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load package
library("Cardinal")

# file name
ImzMLfile <- args[1]
binwidth <- args[2] # usually 0.001
chunksize <- args[3] # 50000
i_threshold <- args[4] # 100 for not normalized datasets
outfile <- args[5]

#-------------------------------------------------------------

# Import MSi dataset (ImzML and ibd!)
msi <- readImzML(ImzMLfile,
                 attach.only = FALSE,
                 as="MSImagingExperiment",
                 resolution = binwidth, # dependent on instrument setup
                 units = "mz") # or ppm

# define fragmentation
pixel       <- as.character(seq(1,ncol(msi),1))
mz          <- as.vector(mz(msi))

# create chunksize list
chunksize_list <- list()
chunksize_list[[1]] <- 1:chunksize
for (i in 2:floor(length(mz)/chunksize)) {
  chunksize_list[[i]] <- chunksize_list[[i-1]]+chunksize
}
chunksize_list[[length(chunksize_list)+1]] <- (floor(length(mz)/chunksize)*chunksize+1):length(mz)

# create mz list
mz_list <- list()
for (i in 1:length(chunksize_list)) {
  mz_list[[i]] <- mz[chunksize_list[[i]]]
}

# create minor matrices
matrix_list <- list()
for (i in 1:length(mz_list)) {
  matrix_list[[i]]            <- spectra(msi)[chunksize_list[[i]],]
  colnames(matrix_list[[i]])  <- pixel
  rownames(matrix_list[[i]])  <- mz_list[[i]]
  matrix_list[[i]]            <- matrix_list[[i]][rowSums(matrix_list[[i]] > i_threshold) >= 1, ] # strips matrix by deleting peaks (rows) under threshold
  matrix_list[[i]]            <- t(matrix_list[[i]])
  mz_list[[i]]                <- as.numeric(colnames(matrix_list[[i]]))
  matrix_list[[i]]            <- rbind(mz_list[[i]], matrix_list[[i]])
}

# fuse matrices
rm(msi, mz, chunksize_list, mz_list)
intensity_matrix <- do.call(cbind, matrix_list)
rm(matrix_list)

# Save data as *.csv
write.table(intensity_matrix,
            file = paste(outfile, ".csv"),
            row.names = TRUE,
            na = "",
            col.names = FALSE,
            sep = ",")



