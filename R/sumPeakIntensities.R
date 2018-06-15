#' Sum total pixel intensities per peak and report as data frame
#'
#' @param d msitsm or msidf; MSI data imported with \code{\link{readMSI}}
#'
#' @return data.frame of peaks and intensities
#' @export

sumPeakIntensities <- function(d) UseMethod("sumPeakIntensities")