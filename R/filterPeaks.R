#' Filter peak list of MSI data
#'
#' Perform peak filtering of MSI dataset by one of three methods. "topX" retains
#' the top X peaks by intensity. "XofTop" retains all peaks above fraction X of
#' the highest peak. "XofTotal" retains the top peaks (by intensity) that
#' cumulatively account for 1-X of the total intensity. For "topX", the value of
#' parameter X should be an integer >= 1, whereas for the latter two, X should
#' be between 0 and 1.
#'
#' @param d either an \code{\link{msimat}} object, or a data.frame with three
#'        columns: \code{peaks}, \code{counts}, and \code{intensities}
#'        representing the mass spectrum of the MSI data.
#' @param how character; Method to use for filtering (see above)
#' @param x numeric; Parameter for the filtering method.
#' @param index logical; If input is a dataframe, and \code{index=TRUE}, then
#'        return will be a vector of indices rather than a data frame.
#'
#' @return Either a data.frame or msimat object, representing subset of the
#'         input data, or a vector of indices if \code{index=TRUE} and the input
#'         is a data.frame.
#' @export

filterPeaks <- function(d,
                        how=c("topX","XofTop","XofTotal"),
                        x=NULL,
                        index=NULL) UseMethod("filterPeaks")
