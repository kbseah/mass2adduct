#' Calculate correlations of pairs of mass peaks from MSI data in chunks
#'
#' Calculate correlations of pairs of mass peaks from MSI data (imported to
#' R with the \code{\link{readMSI}} function). The list of pairs is supplied
#' as a data.frame with the parameter \code{pairs}, but is broken up into "chunks"
#' for processing, to avoid going over a specified memory limit.
#'
#' This function is otherwise equivalent to \code{\link{corrPairsMSI}} with
#' option how="parallel". The memory requirement is estimated conservatively but
#' has only been tested empirically (i.e. the process might go over the limit)
#' in practice.
#'
#' Example usage scenario for this function: Mass differences for all pairwise
#' combinations of masses are tabulated with \code{\link{diffTabulate}}, and the
#' peak pairs corresponding to a specific adduct of interest are extracted with
#' \code{\link{diffGetPeaks}}. Check which peaks are significantly correlated to
#' each other with this function.
#'
#' Correlation is calculated with \code{\link{cor.test}} function from the
#' \code{stats} package, and by default uses the Pearson methone (one-sided, i.e.
#' positive correlations). The Bonferroni correction is applied to the p-values
#' before assessing significance, but the original p-value is reported in column
#' P.value.
#'
#' @param d msimat; MSI data with peaks as columns and pixels as rows,
#'        output from \code{\link{readMSI}} function
#' @param diff massdiff; List of mass differences, parent and putative adduct
#'        ions, as produced by function \code{\link{diffTabulate}}
#' @param p.val numeric; p-value cutoff (before Bonferroni correction) (default: 0.05)
#' @param method string; Method to use for \code{\link{cor.test}} (default: "pearson")
#' @param alternative string; Which alternative for \code{\link{cor.test}} (default: "greater")
#' @param ncores integer; Number of cores if using parallel version. Default is
#'        total number of cores minus one.
#' @param mem.limit integer; Memory limit estimate (in Gb). (default: 5)
#' @param ... Other parameters to pass to \code{\link{cor.test}}
#'
#' @return Object of class data.frame with the following fields:
#'
#'         A - First peak mass in each pair
#'
#'         B - Second peak mass in the pair
#'
#'         Estimate - Estimated correlation
#'
#'         P.value - P-value for the correlation
#'
#'         Significance - Whether p-value meets cutoff (specified by "p.val", with
#'                        Bonferroni correction)
#'
#' @seealso \code{\link{corrSinglePeakMSI}} to calculate correlations of all
#'          peaks in a dataset vs. a single peak of interest
#'
#' @export

corrPairsMSIchunks <- function(d,
                               diff,
                               p.val=0.05,
                               method="pearson",
                               alternative="greater",
                               ncores=NULL,
                               mem.limit=5,
                               ...) UseMethod("corrPairsMSIchunks")
