#' Sum total pixel intensities per peak and report as data frame
#'
#' @param d msimat; MSI data imported with \code{\link{readMSI}}
#'
#' @return data.frame of peaks and intensities
#' @export

sumPeakIntensities.msimat <- function(d) {
    sums <- Matrix::colSums(d[["mat"]])
    counts <- Matrix::colSums(d[["mat"]] != 0) # Count nonzero values per column
    peaks <- d[["peaks"]]
    out <- data.frame(peaks=peaks, counts=counts, intensities=sums)
    return(out)
}
