#' Sum total pixel intensities per peak and report as data frame
#'
#' @param d msidf; MSI data imported with \code{\link{readMSI}}
#'
#' @return data.frame of peaks and intensities
#' @export

sumPeakIntensities.msidf <- function(d) {
    sums <- colSums(d)
    peaks <- names(d)
    out <- data.frame(peaks=peaks, intensities=sums)
    return(out)
}