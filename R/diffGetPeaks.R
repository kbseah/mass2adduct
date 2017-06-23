#' Find pairs of mass peaks corresponding to a specific mass difference
#'
#' @param d data.frame; output from diffTabulate
#' @param mass numeric; mass difference in m/z units
#' @param width numeric; range of values to use
#'
#' @return data.frame of peak pairs where the difference is within the
#' width of the specified mass of interest
#' @export

diffGetPeaks <- function(d, mass=NULL, width=0.01) {
    if (! is.null (mass)) {
        diffLow <- mass - width/2
        diffUpp <- mass + width/2
        output <- d[which(d$diff > diffLow & d$diff <= diffUpp),]
        return (output)
    } else {
        cat ("Error: Mass difference not specified\n")
    }
}
