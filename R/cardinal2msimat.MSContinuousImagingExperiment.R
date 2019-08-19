#' Convert Cardinal MSContinuousImagingExperiment object to msimat object
#'
#' Convert imaging data in MSContinuousImagingExperiment, which is stored
#' internally as a sparse matrix (class sparse_matc object) to an msimat object.
#'
#' @param d MSContinuousImagingExperiment; MSI data processed by Cardinal, must
#'        be processed with peakBin.
#' @export

cardinal2msimat.MSContinuousImagingExperiment <- function(d) {
  # Check if peaks have been binned, otherwise the keys will not correspond to
  # peaks
  if ("peakBin" %in% names(d@processing)) {
    # Check that data are indeed in sparse matrix format
    if (class(spectra(d)) == "matrix") {
      # Initialize vectors for containing data
      peaks <- mz(d) # Binned mass values
      spots <- 1:dim(d)['Pixels'] # Pixel name values
      data <- t(spectra(d))
      peakintensities <- colSums(data)
      # Store as list and declare class as msimat
      out <- list(mat=as.matrix(data),
                  peaks=as.numeric(peaks),
                  spots=spots,
                  peakintensities=peakintensities)
      class(out) <- "msimat"
      # Return msimat object
      return(out)
    } else {
      stop("Matrix not stored internally as sparse_matc, which should be the case for MSProcessedImagingExperiment")
    }
  } else {
    stop("Only applicable to binned data. Please apply peakBin() processing first.")
  }
}
