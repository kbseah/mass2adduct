#' Count peaks corresponding to putative adduct types
#'
#' Adduct matching with \code{\link{adductMatch}} or visualizing the massdiff
#' histogram with \code{\link{hist.massdiff}} will give an overview of the
#' number of ion pairs corresponding to transformations of interest. However,
#' this does not directly indicate the number of actual ions that participate in
#' the transformations, because a given ion can participate in multiple
#' transformations either as parent or adduct.
#'
#' This function takes an annotated massdiff result, where adducts of interest
#' have already been matched to the mass difference values, and tabulates the
#' numbers of parent or adduct ions per adduct type. The total may sum to more
#' than 100%, because a given ion may represent more than one type of
#' chemical transformation.
#'
#' @param diff Object of class massdiff that has been annotated with adduct
#'        types by \code{\link{adductMatch}} function.
#' @param which Whether to enumerate putative adducts or parent ions
#'
#' @return data.frame with number of peaks corresponding to each putative
#'         chemical transformation type.
#' @export

countPeaksFromDiff <- function(diff,
                               which=c("adduct","parent")
                               ) {
    if (which == "adduct") {
        vec <- data.frame(mass=diff$B,matches=diff$matches)
    } else {
        vec <- data.frame(mass=diff$A,matches=diff$matches)
    }
    out <- plyr::ddply(vec,"matches",function(x) data.frame(peaks=length(unique(x$mass))))
    total <- length(unique(vec$mass))
    out <- rbind(out,data.frame(matches="all",peaks=total))
    return(out)
}