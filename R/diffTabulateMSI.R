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

diffTabulateMSI <- function(d, override.limit=FALSE) {
    # Check that d is a data.frame
    if (!is.data.frame(d)) {
        cat ("Error: Input should be a data.frame\n")
    } else {
        # Get list of column names as numeric vector (should be peak mass values)
        d.names <- as.numeric(names(d))
        # Tabulate pairs with diffTabulate
        output <- diffTabulate(d=d.names, override.limit=override.limit)
        return(output)
    }
}