#' Convert Cardinal objects to msimat object
#'
#' Convert imaging data imported with Cardinal v2+ to msimat format for analysis
#' with mass2adduct. The data must be already pre-processed and peak-binned.
#'
#' @param d MSProcessedImagingExperiment; MSI data in Cardinal format, must
#'        be already pre-processed with peakBin.
#'
#' @return Object of class msimat
#'
#' @seealso \code{\link{msimat}}
#' @export

cardinal2msimat <- function(d) UseMethod("cardinal2msimat")
