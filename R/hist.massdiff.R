#' Histogram method for massdiff object
#'
#' @param diff massdiff; Output from function \code{\link{diffTabulate}},
#'        containing three numeric vectors (A, B, diff) representing mass
#'        peaks and their differences
#' @param widthFunc character; function to use to bin mass differences into
#'        histogram. (default "equal", other options to be added)
#' @param width numeric; bin width in m/z units (default 0.01)
#' @param plot logical; whether histogram graphic should be plotted (default FALSE)
#' @param ... Parameters to be passed to hist()
#'
#' @return Object of class hist
#' @seealso \code{\link{diffTabulate}} to generate the mass difference list
#' @export

hist.massdiff <- function(diff, widthFunc="equal", width=0.01, plot=FALSE, ...) {
    if (!is.numeric(diff$diff)) {
        cat("Error: Input mass difference list must be numeric\n")
    } else {
        # Calculate number of breaks for histogram, integer value
        minval <- floor(min(diff$diff,na.rm=TRUE))
        maxval <- ceiling(max(diff$diff,na.rm=TRUE))
        if (widthFunc == "equal") { # equal bin widths
            breaks <- round((maxval - minval)/width, digits=0)
        } # other options TBD
        output <- hist(diff$diff, breaks=breaks, plot=plot, ...)
        return(output)
    }
}