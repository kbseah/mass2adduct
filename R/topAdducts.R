#' Find top adduct candidates from a mass difference histogram
#'
#' Identify mass differences with the most observations and report closest-
#' matching known adducts. This also reports mass differences which do not have
#' a close match to a known adduct, unlike \code{\link{adductMatch}}
#'
#' @param hist histogram of a massdiff object
#' @param add data.frame of adducts (default: "adducts" dataset of mass2adduct)
#' @param n numeric; number of top hits to report (default 20)
#' @param use.bw logical; adduct is reported as match only if it falls in the
#'        same histogram bin (default: TRUE)
#' @param threshold numeric; maximum distance between adduct mass and mass
#'        difference peak to report adduct as a match. use.bw must be FALSE
#'        (default: NULL)
#'
#' @return data.frame of mass bins, counts, closest-matching adduct if known
#' @seealso \code{\link{adductMatch}} to find the closest-matching mass
#'          differences for a list of known adducts (the converse of the
#'          current function)
#' @export

topAdducts <- function(hist, add=adducts, n=20, use.bw=TRUE, threshold=NULL) {
  if (!"massdiffhist" %in% class(hist)) {
    stop("Input must be an object of class massdiffhist")
  }
  if (use.bw) { # Check is using binwidth as the threshold for assigning adduct
    if (!is.null(threshold)) {
      warning("Parameter \"threshold\" ignored when use.bw is TRUE")
    }
    # Get width of histogram bins
    binwidth <- hist$mids[2] - hist$mids[1] # Width of histogram bins
    threshold <- binwidth/2
  } else {
    if (is.null(threshold)) {
      warning("Parameter \"threshold\" not defined although use.bw is FALSE ... setting threshold to 0.005 by default")
      threshold <- 0.005
    }
  }
  # Sorted indices of mass difference bins by counts
  indsort <- order(hist$counts,decreasing=TRUE)[1:n]
  mids <- hist$mids[indsort]
  addIndices <- sapply(mids, function(x) {
      index <- which.min(abs(add$mass - x))
      if ( abs(add$mass[index] - x) < threshold ) {
        return (index)
      } else {
        return (NA)
      }
    })
  output <- data.frame (mids,
                        counts=hist$counts[indsort],
                        density=hist$density[indsort],
                        adduct.name=as.character(add$name[addIndices]),
                        adduct.formula=as.character(add$formula[addIndices]),
                        adduct.mass=add$mass[addIndices])
  return(output)
}
