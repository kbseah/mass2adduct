#' Overlay adduct parental and adduct masses on plot of msimat object
#'
#' On a mass spectrum plot (from an msimat object), plot points representing
#' either parental or adduct masses as an overlay. By default only those flagged
#' as having statistically significant spatial correlation will be plotted. The
#' massdiff object should have been processed with \code{\link{adductMatch}} and
#' \code{\link{corrPairsMSI}} to annotate with adduct matches and correlation
#' test results. 
#'
#' @param d msimat; Object of class msimat created by function
#'        \code{\link{msimat}}
#' @param diff data.frame; Mass difference correlation table produced by
#'        function \code{\link{corrPairsMSI}}
#' @param which string; Either plot "parent" or "adduct" peaks as overlay
#' @param invert logical; Overlay the points which are NOT adducts
#' @param signif logical; Only plot points which have statistically significant
#'        correlation
#' @param pch Plot character to use, passed to \code{points}
#' @param cex Character expansion parameter, passed to \code{points}
#' @param col Color for overlay points, passed to \code{points}
#' @param ... Other parameters to pass to \code{points}
#'
#' @seealso \code{\link{corrPairsMSI}} calculate correlations for mass pairs
#' @seealso \code{\link{massdiff}} tabulate all possible mass pairs
#' @export

pointsAdducts <- function(d, diff, which=c("adduct","parent"), invert=FALSE, signif=TRUE, pch=20, cex=0.5, col="red", ...) {
  if (signif) {
    diff <- diff[which(diff$Significance == 1),]
  }
  if (which[1] == "adduct") {
    shortlist <- diff$B
  } else if (which[1] == "parent") {
    shortlist <- diff$A
  }
  # It is possible that there are duplicates in shortlist because there may be
  # more than one parent/adduct pair with a mass difference within the
  # tolerated range
  if (invert) {
    peaksshortlist <- d[["peaks"]][which(!d[["peaks"]] %in% shortlist)]
    intensities <- d[["peakintensities"]][which(!d[["peaks"]] %in% shortlist)]
  } else {
    peaksshortlist <- d[["peaks"]][which(d[["peaks"]] %in% shortlist)]
    intensities <- d[["peakintensities"]][which(d[["peaks"]] %in% shortlist)]
  }
  points(x=peaksshortlist,
         y=intensities,
         pch=pch,
         cex=cex,
         col=col,
         ...)
}
