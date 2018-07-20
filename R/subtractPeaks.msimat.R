#' Remove peaks from an msimat object within specified mass accuracy
#'
#' Remove peaks, e.g. known contaminants, from an msimat object. The peaks are
#' specified by their m/z values, and peak intensities are not taken into
#' account. Mass accuracies (expressed as ppm values) for both the msimat
#' dataset and the removal peaklist must be specified, and will be used to
#' determine the margin of error for excluding a given peak.
#'
#' @param input numeric; msimat object containing peaks to be removed
#' @param v numeric; List of peaks to be removed from d
#' @param ppm1 numeric; Mass accuracy for peaklist v1, in ppm
#' @param ppm2 numeric; Mass accuracy for peaklist v2, in ppm
#'
#' @return msimat object
#'
#' @seealso \code{\link{subtractPeaks.numeric}}
#' @seealso \code{\link{filterPeaks}}
#' @export

subtractPeaks.msimat <- function(input,
                                 vec,
                                 ppm1,
                                 ppm2) {
  idx <- subtractPeaks(input$peaks,vec,ppm1,ppm2,index=TRUE)
  out <- list(mat=input$mat[,idx], # Subset of original matrix
              peaks=input$peaks[idx],
              spots=input$spots, # Unchanged
              peakintensities=input$peakintensities[idx])
  class(out) <- "msimat"
  return(out)
}