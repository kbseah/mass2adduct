#' Tabulate pairwise combinations of mass peaks
#'
#' Takes all pairwise combinations of masses from a simple numeric vector
#' that represent the mass values from a mass spectrometry data set.
#'
#' The number of possible pairs grows combinatorially and the calculation will
#' fail for >65536 values (the function uses a 16-bit integer internally).
#' 
#' @param d numeric vector representing m/z values of mass peaks
#'
#' @return data.frame of all pairs of masses and their respective differences
#' @seealso \code{\link{diffTabulateMSI}} to tabulate mass differences from
#'          MS imaging data object imported with \code{\link{readMSI}} function,
#'          \code{\link{diffHist}} to plot histogram of mass differences
#' @export

diffTabulate <- function(d) {
    # Check that d is numeric 
    if (!is.numeric(d)) {
        cat ("Error: Input should be numeric\n")
    } else {
        len <- length(d)
        cat ("Input is numeric vector of length", len, "with", choose(len,2), "possible pairs\n")
        if (len > 65536) { # Warning if vector is too long 
            cat ("Warning: Input is longer than 65536 elements\n")
            cat ("... perhaps you should filter your peak list\n")
        } else {
            y <- comb2M(d) # Alternative to R built-in combn function
            # Format results as data frame with column labels
            y <- as.data.frame (t(y))
            names(y) <- c("A","B","diff")
            return(y)
        }
    }
}