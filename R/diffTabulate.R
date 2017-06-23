#' Tabulate pairwise combinations of mass peaks
#'
#' Takes all pairwise combinations of masses from a simple numeric vector
#' that represent the mass values from a mass spectrometry data set.
#'
#' The number of possible pairs grows combinatorially and so the calculation
#' can take a very long time for more than 5000 mass values. Use the
#' \code{override.limit} parameter to override the default limit of 5000.
#' 
#' @param d numeric vector representing m/z values of mass peaks
#' @param override.limit logical; continue with computation even if vector
#'        length is >5000 (default: FALSE)
#'
#' @return data.frame of all pairs of masses and their respective differences
#' @seealso \code{\link{diffTabulateMSI}} to tabulate mass differences from
#'          MS imaging data object imported with \code{\link{readMSI}} function,
#'          \code{\link{diffHist}} to plot histogram of mass differences
#' @export

diffTabulate <- function(d, override.limit=FALSE) {
    # Check that d is numeric 
    if (!is.numeric(d)) {
        cat ("Error: Input should be numeric\n")
    } else {
        len <- length(d)
        cat ("Input is numeric vector of length", len, "with", choose(len,2), "possible pairs\n")
        if (len > 5000 & override.limit==FALSE) { # Warning if vector is too long to compute in reasonable time
            cat ("Warning: Input is longer than 5000 elements, computation may be long and memory-intensive\n")
            cat ("Use option override.limit=TRUE to compute anyway\n")
        } else {
            y <- combn(d, 2, function(x) return (c(x[1], x[2], round(abs(x[2] - x[1]), digits = 4))))
            # Format results as data frame with column labels
            y <- as.data.frame (t(y))
            names(y) <- c("A","B","diff")
            return(y)
        }
    }
}