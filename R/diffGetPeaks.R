#' Find pairs of mass peaks corresponding to a specific mass
#'
#' @param diff massdiff; output from diffTabulate
#' @param by string; Subset the mass by mass difference (putative adduct
#'                   transformation) ("diff") or by parent ion ("parent")?
#' @param mass numeric; mass of putative adduct or parent ion in m/z units
#' @param width numeric; range of values to use
#'
#' @return subset of the original massdiff object where the mass difference or
#'         parent ion is close to the specified mass of interest (within the
#'         specified width)
#' @export

diffGetPeaks <- function(diff, by="diff", mass=NULL, width=0.001) {
    if (class(diff) != "massdiff") {
        cat ("Error: Input to parameter diff must be an object of class massdiff\n")
    } else {
        if (! is.null (mass)) {
            idx <- diffGetPeaksIndex(diff=diff,by=by,mass=mass,width=width)
            newA <- diff$A[idx]
            newB <- diff$B[idx]
            newdiff <- diff$diff[idx]
            output <- list(A=newA,
                           B=newB,
                           diff=newdiff)
            class(output) <- "massdiff"
            return (output)
        } else {
            cat ("Error: Mass difference not specified\n")
        }
    }
}
