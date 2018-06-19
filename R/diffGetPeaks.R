#' Find pairs of mass peaks corresponding to a specific mass
#'
#' From a set of all possible mass pairs, find all pairs corresponding to a
#' specific mass difference, which might represent a molecular transformation
#' of interest.
#'
#' @param diff massdiff; output from \code{\link{massdiff}}
#' @param by string; Subset the mass by mass difference (putative adduct
#'                   transformation) ("diff") or by parent ion ("parent")?
#' @param mass numeric; mass of putative adduct or parent ion in m/z units
#' @param width numeric; range of values to use
#'
#' @return Subset of the original massdiff object where the mass difference or
#'         parent ion is close to the specified mass of interest (within the
#'         specified width)
#'
#' @export

diffGetPeaks <- function(diff, by=c("diff","parent"), mass=NULL, width=0.001) {
    if (!"massdiff" %in% class(diff)) {
        stop ("Input to parameter diff must be an object of class massdiff")
    }
    if (is.null (mass)) {
        stop("Mass difference not specified")
    }
    idx <- diffGetPeaksIndex(diff=diff,by=by,mass=mass,width=width)
    newA <- diff$A[idx]
    newB <- diff$B[idx]
    newdiff <- diff$diff[idx]
    output <- list(A=newA,
                   B=newB,
                   diff=newdiff)
    class(output) <- c("massdiff","data.frame")
    return (output)
}
