#' Find pairs of mass peaks corresponding to a specific mass difference (putative adduct)
#'
#' @param d data.frame; output from diffTabulate
#' @param by string; Subset the mass by mass difference (putative adduct
#'                   transformation) ("diff") or by parent ion ("parent")?
#' @param mass numeric; mass of putative adduct or parent ion in m/z units
#' @param width numeric; range of values to use
#'
#' @return subset of the original data frame d where the mass difference or
#'         parent ion is close to the specified mass of interest (within the
#'         specified width)

diffGetPeaksIndex <- function(d, by="diff", mass=NULL, width=0.01) {
    if (! is.null (mass)) {
        diffLow <- mass - width/2
        diffUpp <- mass + width/2
        if (by=="diff") {
            output <- which(d$diff > diffLow & d$diff <= diffUpp)
        } else if (by == "parent") {
            output <- which(d$A > diffLow & d$A <= diffUpp)
        }
        return (output)
    } else {
        cat ("Error: Mass difference not specified\n")
    }
}
