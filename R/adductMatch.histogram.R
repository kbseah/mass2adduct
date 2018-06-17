#' Find closest-matching mass difference bin for known adducts
#'
#' For use with internal dataset "adducts" and other similar data
#'
#' @param x histogram object produced by \code{\link{hist.massdiff}}
#' @param add data.frame of adduct masses (default: "adducts" dataset in package)
#' @param density logical; return density instead of counts? (default: FALSE)
#'
#' @return data.frame of adducts with counts/density from closest-matching bin
#'         and their corresponding quantiles in mass difference bins
#' @seealso \code{\link{topAdducts}} to find the closest-matching adducts for
#'          the mass differences with the highest counts (the converse of the
#'          current function)
#' @export

adductMatch.histogram <- function(x,add=adducts,density=FALSE) {
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
