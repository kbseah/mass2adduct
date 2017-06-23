#' Find closest-matching mass difference bin for known adducts
#'
#' For use with internal dataset "adducts" and other similar data
#'
#' @param hist histogram object produced by diffHist()
#' @param add data.frame of adduct masses (default: "adducts" dataset in package)
#' @param density logical; return density instead of counts? (default: FALSE)
#'
#' @return data.frame of adducts with counts/density from closest-matching bin
#'         and their corresponding quantiles in mass difference bins
#' @seealso \code{\link{topAdducts}} to find the closest-matching adducts for
#'          the mass differences with the highest counts (the converse of the
#'          current function)
#' @export

adductMatch <- function(hist,add=adducts,density=FALSE) {
    if (class(hist) != "histogram") {
        cat("Error: Object hist has to be a histogram\n")
    } else {
        binwidth <- hist$mids[2] - hist$mids[1] # Width of histogram bins
        indices <- sapply(add$mass, function(x) {
            index <- which.min(abs(hist$mids - x))
            # Check for values that fall outside the closest bin!
            if ( abs(hist$mids[index] - x) < binwidth/2 ) {
                return (index)
                } else { return (NA) }})
        if (density) { # Report densities
            # Get quantiles for bin of each adduct
            diffqnts <- ecdf (hist$density)
            quantiles <- diffqnts(hist$density[indices])
            output <- data.frame(add, density=hist$density[indices], quantiles)
        } else { # Report counts
            diffqnts <- ecdf (hist$counts)
            quantiles <- diffqnts (hist$counts[indices])
            output <- data.frame (add, counts=hist$counts[indices],quantiles)
        }
        return(output)
    }
}