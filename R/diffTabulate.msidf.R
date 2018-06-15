#' Tabulate pairwise combinations of mass peaks from imported MSI data
#'
#' Takes all pairwise combinations of mass peaks from MS imaging data imported
#' with the function \code{\link{readMSI}}. 
#'
#' The number of possible pairs grows combinatorially and so the calculation
#' can take a very long time for more than 5000 mass values. Use the
#' \code{override.limit} parameter to override the default limit of 5000.
#' 
#' @param d data.frame; MS imaging data imported with \code{readMSI}
#' @param override.limit logical; continue with computation even if vector
#'        length is >5000 (default: FALSE)
#'
#' @return data.frame of all pairs of masses and their respective differences
#' @seealso \code{\link{diffTabulate}} for tabulating from a simple numeric
#'          list of masses, \code{\link{diffHist}} to plot histogram of mass
#'          differences, \code{\link{readMSI}} to read MS imaging data
#' @export

diffTabulate.msidf <- function(d) {
    # Check that d is a data.frame
    if (class(d) == "msidf") {
        # Get list of column names as numeric vector (should be peak mass values)
        peaklist <- as.numeric(names(d))
        # Tabulate pairs with diffTabulate
        output <- diffTabulate.numeric(d=peaklist)
        return(output)
    } else {
        cat ("Error: Input should be a data.frame\n")
    }
}