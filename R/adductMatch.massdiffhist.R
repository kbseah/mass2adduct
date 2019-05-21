#' Report closest-matching adduct type for a mass difference histogram
#'
#' Given a histogram of mass differences (from applying \code{hist} on a massdiff
#' object), report the bins most closely matching a list of potential molecular
#' transformations (adduct types). The built-in data sets \code{\link{adducts}}
#' and \code{\link{adducts2}} are examples of such lists of potential adducts.
#' This function is designed for a quick exploratory overview of an MSI data set,
#' especially when the number of peaks is large. To report each peak pair that
#' has a match to the list of potential adducts, apply \code{\link{adductMatch}}
#' directly to a massdiff object.
#'
#' @param x massdiffhist object produced by \code{\link{hist.massdiff}}
#' @param add data.frame of adduct masses (default: "adducts" dataset in package)
#' @param density logical; return density instead of counts? (default: FALSE)
#'
#' @return data.frame of adducts with counts/density from closest-matching bin
#'         and their corresponding quantiles in mass difference bins
#' @seealso \code{\link{adductMatch.massdiff}} reports closest-matching adducts
#'          for each individual peak pair.
#' @seealso \code{\link{topAdducts}} to find the closest-matching adducts for
#'          the mass differences with the highest counts (the converse of the
#'          current function)
#' @export

adductMatch.massdiffhist <- function(x,add=mass2adduct::adducts, density=FALSE) {
  binwidth <- x$mids[2] - x$mids[1] # Width of histogram bins
  indices <- sapply(add$mass, function(z) {
    index <- which.min(abs(x$mids - z))
    # Check for values that fall outside the closest bin!
    if ( abs(x$mids[index] - z) < binwidth/2 ) {
      return (index)
    } else {
      return (NA)
    }
  })
  if (density) { # Report densities
    # Get quantiles for bin of each adduct
    diffqnts <- ecdf (x$density)
    quantiles <- diffqnts(x$density[indices])
    output <- data.frame(add, density=x$density[indices], quantiles)
  } else { # Report counts
    diffqnts <- ecdf (x$counts)
    quantiles <- diffqnts (x$counts[indices])
    output <- data.frame (add, counts=x$counts[indices],quantiles)
  }
  return(output)
}
