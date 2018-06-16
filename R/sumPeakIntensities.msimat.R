#' Sum total pixel intensities per peak and report as data frame
#'
#' @param d msitsm; MSI data imported with \code{\link{readMSI}}
#'
#' @return data.frame of peaks and intensities
#' @export

sumPeakIntensities.msimat <- function(d) {
    sums <- Matrix::colSums(d[["mat"]])
    peaks <- d[["peaks"]]
    out <- data.frame(peaks=peaks, intensities=sums)
    return(out)
}
