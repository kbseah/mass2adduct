#' Find pairs of mass peaks corresponding to a specific mass
#'
#' From a set of all possible mass pairs, find all pairs corresponding to a
#' specific mass difference, which might represent a molecular transformation
#' of interest.
#'
#' There are two options for setting a cutoff to match a given mass to a given
#' target value. The cutoff can be either proportional to the mass (and
#' expressed as parts per million \code{ppm}) or a flat absolute value (in
#' milliDaltons \code{mDa}). ppm is used by default. If a value for mDa is
#' specified, then any value given to ppm is ignored. The ppm or mDa values are
#' usually reported by the peak-picking software used.
#'
#' @param diff massdiff; output from \code{\link{massdiff}}
#' @param mass numeric; mass of putative adduct or parent ion in m/z units
#' @param ppm numeric; the peak width (uncertainty) of the calculated peaks in
#'        proportional units, as parts per million (default: 2)
#' @param mDa numeric; the peak width (uncertainty) of the calculated peaks in
#'        absolute mass units, as milliDaltons (default: NULL)
#'
#' @return Subset of the original massdiff object where the mass difference or
#'         parent ion is close to the specified mass of interest within the
#'         uncertainty specified
#'
#' @export

diffGetPeaks <- function(diff, mass=NULL, ppm=2, mDa=NULL) {
  if (!"massdiff" %in% class(diff)) {
    stop ("Input to parameter diff must be an object of class massdiff")
  }
  if (is.null (mass)) {
    stop("Mass difference not specified")
  }
  if (!is.null(mDa)) {
    Ad <- rep(mDa * 1e-3, times=length(diff$A))
    Bd <- rep(mDa * 1e-3, times=length(diff$B))
    #idx <- diffGetPeaksIndex(diff=diff,by="diff",mass=mass,width=width)
  } else {
    Ad <- diff$A * ppm * 1e-6
    Bd <- diff$B * ppm * 1e-6
  }

  diff$delta <- sqrt(Ad**2 + Bd**2) # Uncertainties add in quadrature
  matchdiff <- abs(diff$diff - mass)
  idx <- which(matchdiff < diff$delta)
    
  newA <- diff$A[idx]
  newB <- diff$B[idx]
  newdiff <- diff$diff[idx]
  output <- data.frame(A=newA,
                       B=newB,
                       diff=newdiff)
  class(output) <- c("massdiff","data.frame")
  return (output)
}
