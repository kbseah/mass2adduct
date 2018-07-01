#' Identify peak on mass spectrum plot
#'
#' Interactively identify peak on mass spectrum plot produced by
#' \code{\link{plot.msimat}}. Click on the plot to identify the closest mass
#' peak. Right-click to stop on most systems.
#'
#' @param d msimat object used to plot the mass spectrum
#' @param plot logical; display label beside identified peak with the peak mass
#'        value?
#' @param give.names logical; return peak masses instead of their indices
#' @param ... Other arguments to \code{\link{identify}}
#'
#' @return Vector of peak masses if \code{give.names=TRUE}, otherwise indices
#'         for the identified peaks
#' @seealso \code{\link{plot.msimat}}
#' @export

identify.msimat <- function(d,
                            plot=TRUE,
                            give.names=TRUE,
                            ...) {
  idx <- identify(x=d$peaks,
                  y=d$peakintensities,
                  labels=as.character(d$peaks),
                  ...)
  if (give.names) {
    return(as.character(d$peaks[idx]))
  } else {
    return(idx)
  }

}
