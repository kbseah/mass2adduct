#' Tabulate pairwise combinations of mass peaks
#'
#' Takes all pairwise combinations of masses from either an MSI data set that
#' was imported with the function \code{\link{msimat}}, or from a simple
#' numeric vector of masses.
#'
#' The number of possible pairs grows combinatorially and the calculation will
#' fail for >65536 values (the function uses a 16-bit integer internally).
#'
#' @param d Either an MSI data set of class msimat, imported with the function
#'          \code{\link{msimat}}, or a numeric vector representing m/z values
#'          of mass peaks
#'
#' @return Object of class \code{massdiff} with three elements:
#'
#'                    A - m/z values for peak pairs; for a given pair, A will
#'                        be lower than B. Requries that the peak list in the
#'                        input be sorted ascending numerically by mass.
#'
#'                    B - m/z values for peak pairs
#'
#'                    diff - their respective differences
#'
#' @seealso \code{\link{hist.massdiff}} to plot histogram of mass differences
#' @seealso \code{\link{adductMatch}} to match mass differences to possible
#'          molecular transformations
#' @seealso \code{\link{corrPairsMSI}} to test whether given peak pairs are
#'          correlated, using MSI intensity data from massdiff and msimat objects
#'
#' @examples
#' d <- msimat(csv=system.file("extdata","msi.csv",package="mass2adduct"),sep=";")
#' d.diff <- massdiff(d) # d is an msimat object
#' print(d.diff) # Show a summary of massdiff object
#'
#' @export

massdiff <- function(d) UseMethod("massdiff")
