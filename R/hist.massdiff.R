#' Histogram method for massdiff object
#'
#' Bin mass differences from a massdiff object produced by
#' \code{\link{massdiff}} into a histogram, and plot if requested.
#'
#' @param diff massdiff; Output from function \code{\link{massdiff}},
#'        containing three numeric vectors (A, B, diff) representing mass
#'        peaks and their differences
#' @param widthFunc character; function to use to bin mass differences into
#'        histogram. (default "equal", other options to be added)
#' @param width numeric; bin width in m/z units (default 0.01)
#' @param ... Other options to be passed to hist()
#'
#' @return Object of classes hist and massdiffhist
#' @seealso \code{\link{massdiff}} to generate the mass difference list
#' @seealso \code{\link{plot.massdiffhist}} for plotting and annotating the
#'          resulting object
#' @export

hist.massdiff <- function(diff, widthFunc="equal", width=0.01, ...) {
  if (!is.numeric(diff$diff)) {
    stop("Input mass difference list must be numeric")
  }
  # Calculate number of breaks for histogram, integer value
  minval <- floor(min(diff$diff,na.rm=TRUE))
  maxval <- ceiling(max(diff$diff,na.rm=TRUE))
  if (widthFunc == "equal") { # equal bin widths
    breaks <- round((maxval - minval)/width, digits=0)
  } # other options TBD
  output <- hist(diff$diff, breaks=breaks, plot=FALSE, ...)
  class(output) <- c("histogram","massdiffhist")
  return(output)
}
