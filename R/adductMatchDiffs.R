#' Find closest-matching known adducts for a set of mass differences
#'
#' For use with internal dataset "adducts" and other similar data
#'
#' @param diff data.frame of mass differences produced by \code{\link{diffTabulate}}
#' @param add data.frame of adduct masses (default: "adducts" dataset in package)
#' @param width numeric; the precision to which to match the mass value (default: 0.01)
#'
#' @return data.frame of mass differences given to argument diff, with an
#'         additional column reporting the closest matches
#' @seealso \code{\link{topAdducts}} to find the closest-matching adducts for
#'          the mass differences with the highest counts (the converse of the
#'          current function)
#' @export

adductMatchDiffs <- function(diff,add=adducts,width=0.01) {
    if (class(diff) != "data.frame") {
        cat("Error: Object hist has to be a data.frame produced by function diffTabulate\n")
    } else {
        diff$match <- rep(NA,length(diff$diff))
        for (i in 1:length(add$mass)) {
            # Using a loop because number of adducts are few
            idx <- diffGetPeaksIndex(diff=diff, by="diff", mass=add$mass[i], width=width)
            diff$match[idx] <- add$name[i]
        }
        return (diff)
    }
}
