#' Report closest-matching adduct type for a list of mass differences
#'
#' Given a set of mass differences produced by function \code{\link{massdiff}},
#' report for each pair of masses the closest-matching molecular transformation
#' from a list of such potential adduct types. The built-in data sets
#' \code{\link{adducts}} and \code{\link{adducts2}} are examples of such lists
#' of potential adducts.
#'
#' @param x massdiff; Mass differences produced by \code{\link{massdiff}}
#' @param add data.frame of adduct masses (default: "adducts" dataset in package)
#' @param width numeric; the precision to which to match the mass value (default: 0.001)
#'
#' @return Object of class massdiff, with additional element "matches" reporting
#'         the closest matches
#' @seealso \code{\link{adductMatch.histogram}} gives an overview for large data
#'          sets by matching closest adducts to mass differences already binned
#'          into a histogram.
#' @seealso \code{\link{topAdducts}} to find the closest-matching adducts for
#'          the mass differences with the highest counts (the converse of the
#'          current function)
#'
#' @examples
#' d <- msimat(csv=system.file("extdata","msi.csv",package="mass2adduct"),sep=";")
#' d.diff <- massdiff(d) # Calculate mass differences from imported MSI data
#' d.diff.annotate <- adductMatch(d.diff,add=adducts2) # Report closest matches
#'
#' @export

adductMatch.massdiff <- function(x,add=adducts,width=0.001) {
    #matches <- rep(NA,length(x$diff))
    indices <- vector()
    matches <- vector()
    for (i in 1:length(add$mass)) {
        # Using a loop because number of adducts are few
        idx <- diffGetPeaksIndex(diff=x, by="diff", mass=add$mass[i], width=width)
        indices <- c(indices, idx)
        matches <- c(matches, rep(as.character(add$name[i]),length(idx)))
    }
    output <- data.frame(A=x$A[indices],
                         B=x$B[indices],
                         diff=x$diff[indices],
                         matches=matches)
    class(output) <- c("massdiff","data.frame")
    row.names(output) <- NULL
    return(output)
}
