#' Find pairs of mass peaks corresponding to a specific mass
#'
#' From a set of all possible mass pairs, find all pairs corresponding to a
#' specific mass difference, which might represent a molecular transformation
#' of interest.
#'
#' @param diff massdiff; output from \code{\link{massdiff}}
#' @param mass numeric; mass of putative adduct or parent ion in m/z units
#' @param ppm numeric; mass resolution of the original measurements, to
#'        determine the cutoff for the match (default: 2)
#' @param width numeric; Instead of a sliding ppm cutoff, use an absolute m/z 
#'        difference cutoff to call matches (default: NULL)
#'
#' @return Subset of the original massdiff object where the mass difference or
#'         parent ion is close to the specified mass of interest (within the
#'         specified width)
#'
#' @export

diffGetPeaks <- function(diff, mass=NULL, ppm=2, width=NULL) {
  if (!"massdiff" %in% class(diff)) {
    stop ("Input to parameter diff must be an object of class massdiff")
  }
  if (is.null (mass)) {
    stop("Mass difference not specified")
  }
  if (!is.null(width)) {
    idx <- diffGetPeaksIndex(diff=diff,by="diff",mass=mass,width=width)
  } else {
    Ad <- diff$A * ppm * 1e-6
    Bd <- diff$B * ppm * 1e-6
    diff$delta <- sqrt(Ad**2 + Bd**2) # Uncertainties add in quadrature
    matchdiff <- abs(diff$diff - mass)
    idx <- which(matchdiff < diff$delta)
  }
  
  newA <- diff$A[idx]
  newB <- diff$B[idx]
  newdiff <- diff$diff[idx]
  output <- data.frame(A=newA,
                       B=newB,
                       diff=newdiff)
  class(output) <- c("massdiff","data.frame")
  return (output)
}
