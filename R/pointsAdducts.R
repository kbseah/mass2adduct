#' Overlay adduct parental and adduct masses on plot of msimat object
#'
#' On a mass spectrum plot (from an msimat object), plot points representing
#' either parental or adduct masses as an overlay. 
#'
#' @param d msimat; Object of class msimat created by function
#'        \code{\link{msimat}}
#' @param diff data.frame; Mass difference correlation table produced by
#'        function \code{\link{corrPairsMSI}}
#' @param which string; Either plot "parental" or "adduct" peaks as overlay
#' @param signif logical; Only plot points which have statistically significant
#'        correlation
#'
#' @seealso \code{\link{corrPairsMSI}} calculate correlations for mass pairs
#' @seealso \code{\link{massdiff}} tabulate all possible mass pairs
#' @export

pointsAdducts <- function(d, diff, which=c("adduct","parental"), signif=TRUE, pch=20, cex=0.5, col="red", ...) {
  if (signif) {
    diff <- subset(diff, Significance == 1)
  }
  if (which[1] == "adduct") {
    shortlist <- diff$B
  } else if (which[1] == "parental") {
    shortlist <- diff$A
  }
  # It is possible that there are duplicates in shortlist because there may be
  # more than one parent/adduct pair with a mass difference within the
  # tolerated range
  peaksshortlist <- d[["peaks"]][which(d[["peaks"]] %in% shortlist)]
  intensities <- d[["peakintensities"]][which(d[["peaks"]] %in% shortlist)]
  points(x=peaksshortlist,
         y=intensities,
         pch=pch,
         cex=cex,
         col=col,
         ...)
}
