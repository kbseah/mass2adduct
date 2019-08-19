#' Convert Cardinal MSProcessedImagingExperiment object to msimat object
#'
#' Convert imaging data in MSProcessedImagingExperiment, which is stored
#' internally as a sparse matrix (class sparse_matc object) to an msimat object.
#'
#' @param d MSProcessedImagingExperiment; MSI data processed by Cardinal, must
#'        be processed with peakBin.
#' @export

cardinal2msimat.MSProcessedImagingExperiment <- function(d) {

  # Check if peaks have been binned, otherwise the keys will not correspond to
  # peaks
  if ("peakBin" %in% names(d@processing)) {
    # Check that data are indeed in sparse matrix format
    if (class(spectra(d)) == "sparse_matc") {
      # Initialize vectors for containing data
      peaks <- spectra(d)@keys # Binned mass values
      spots <- 1:dim(d)['Pixels'] # Pixel name values
      cols <- vector()
      rows <- vector()
      vals <- vector()
      # For each spot, parse the relevant vector
      for (i in 1:length(spots)) {
        vals <- c(vals, spectra(d)@data$values[[i]])
        rows <- c(rows, rep(i-1, length(spectra(d)@data$values[[i]])))
        peaksvec <- spectra(d)@data$keys[[i]]
        peaksidx <- sapply(peaksvec, function(p) which (peaks==p) - 1)
        cols <- c(cols, peaksidx) # -1 to convert to 0-based indexing
      }
    } else {
      stop("Matrix not stored internally as sparse_matc, which should be the case for MSProcessedImagingExperiment")
    }
  }
  else {
    stop("Only applicable to binned data. Please apply peakBin() processing first.")
  }

  # Convert to sparseMatrix object
  tsm <- Matrix::sparseMatrix(i=rows,
                              j=cols,
                              x=vals,
                              index1=FALSE,
                              giveCsparse=FALSE,
                              check=TRUE)
  peakintensities <- Matrix::colSums(tsm)
  # Store as list and declare class as msimat
  data <- list(mat=tsm,
               peaks=peaks,
               spots=spots,
               peakintensities=peakintensities)
  class(data) <- "msimat"
  # Return msimat object
  return(data)
}
