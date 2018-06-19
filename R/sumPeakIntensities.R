#' Sum total pixel intensities per peak and report as data frame
#'
#' @param d msimat; MSI data imported with \code{\link{msimat}}
#'
#' @seealso \code{\link{analyzeIntensityCutoffsCumul}} and
#'          \code{\link{analyzeIntensityCutoffsDistr}} to visualize and explore
#'          possible peak filtering cutoffs, given the output from this function
#' @return data.frame of peaks and intensities
#' @export

sumPeakIntensities <- function(d) {
    if (class(d) != "msimat") {
        stop("Input should be an object of class msimat")
    }
    sums <- Matrix::colSums(d[["mat"]])
    counts <- Matrix::colSums(d[["mat"]] != 0) # Count nonzero values per column
    peaks <- d[["peaks"]]
    out <- data.frame(peaks=peaks, counts=counts, intensities=sums)
    return(out)
}
