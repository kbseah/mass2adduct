#' Plot method for object msimat
#'
#' Plots a graph of total intensity per mass peak
#'
#' @param d msimat; MSI imaging object created by function \code{\link{msimat}}
#' @param type What type of plot to be drawn, passed to \code{\link{plot}}
#' @param xlab Label for x axis
#' @param ylab Label for y-axis
#' @param ... Other parameters to pass to \code{plot}
#'
#' @export

plot.msimat <- function(d, type="l", xlab="m/z", ylab="Total intensity", ...) {
  plot(x=d[["peaks"]],
       y=d[["peakintensities"]],
       type=type,
       xlab=xlab,
       ylab=ylab,
       ... )
}
