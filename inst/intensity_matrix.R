# Load package
library("Cardinal")

# Import MSi dataset (ImzML and ibd!)
file <- "MetaMix_SDHBp_60_900_pos_A28_100um_89x74" # file name of *.imzML and *.ibd, both need to be present
msiset <- readImzML(file,
                   attach.only = FALSE,
                   as="MSImageSet",
                   resolution = 0.0005, # 0.0005 results in a bin width of 0.001 Da
                   units = "mz")
dim(msiset) # get dimensions of the dataset

# Create pixel vector
pixel <- as.character(seq(1,
                      ncol(msiset),
                      1))

# Create m/z vector
mz <- as.vector(mz(msiset))

# Build intensity matrix
df <- spectra(msiset)[,1]
for (i in 2:ncol(spectra(msiset))) {
  df <- cbind(df, spectra(msiset)[,i])
}
threshold <- 100 # set threshold
colnames(df)  <- pixel
rownames(df)  <- mz
df2           <- df[rowSums(df > threshold) >= 1, ] # strip matrix by deleting peaks under threshold
df3           <- t(df2)
mz2           <- as.numeric(colnames(df3))
df4           <- rbind(mz2, df3)

# Save data as *.csv
write.table(df4,
            file = file, # as previous filename
            row.names = TRUE,
            na = "",
            col.names = FALSE,
            sep = ",")