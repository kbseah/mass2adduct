#' Find closest-matching known adducts for a set of mass differences
#'
#' For use with internal dataset "adducts" and other similar data
#'
#' @param x massdiff; Mass differences produced by \code{\link{diffTabulate}}
#' @param add data.frame of adduct masses (default: "adducts" dataset in package)
#' @param width numeric; the precision to which to match the mass value (default: 0.01)
#'
#' @return data.frame of mass differences given to argument diff, with an
#'         additional column reporting the closest matches
#' @seealso \code{\link{topAdducts}} to find the closest-matching adducts for
#'          the mass differences with the highest counts (the converse of the
#'          current function)
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
    return(output)
}
