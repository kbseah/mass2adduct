#' Points method for object msimat
#'
#' Points overlay of graph of total intensity per mass peak
#'
#' @param d msimat; MSI imaging object created by function \code{\link{msimat}}
#' @param type What type of plot to be drawn, passed to \code{\link{points}}
#' @param col color for plot
#' @param ... Other parameters to pass to \code{plot}
#'
#' @seealso \code{\link{plot.msimat}}
#' @export

points.msimat <- function(d, type="h", col="blue", ...) {
  points(x=d[["peaks"]],
         y=d[["peakintensities"]],
         type=type,
         col=col,
         ... )
}
