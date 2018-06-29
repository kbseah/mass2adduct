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
#' @param ppm numeric; the mass resolution of the original measurement
#'        (default: 2)
#'
#' @return Object of class massdiff, with additional element "matches" reporting
#'         the closest matches
#' @seealso \code{\link{adductMatch.massdiffhist}} gives an overview for large data
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

adductMatch.massdiff <- function(x, add=adducts, ppm=2) {
  # For each mass pair calculate the mass difference tolerance
  Ad <- x$A * ppm * 1e-6
  Bd <- x$B * ppm * 1e-6
  x$delta <- sqrt(Ad**2 + Bd**2) # Uncertainties add in quadrature
  
  indices <- vector()
  matches <- vector()
  for (i in 1:length(add$mass)) {
    # Using a loop because number of adducts are few
    matchdiff <- abs(x$diff - add$mass[i])
    idx <- which(matchdiff < x$delta)
    indices <- c(indices, idx)
    matches <- c(matches, rep(as.character(add$name[i]),length(idx)))
  }
  output <- data.frame(A=x$A[indices],
                       B=x$B[indices],
                       diff=x$diff[indices],
                       delta=x$delta[indices],
                       matches=matches)
  class(output) <- c("massdiff","data.frame")
  row.names(output) <- NULL
  return(output)
}
