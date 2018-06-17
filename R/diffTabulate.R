#' Tabulate pairwise combinations of mass peaks
#'
#' Takes all pairwise combinations of masses from either an MSI data set that
#' was imported with the function \code{\link{readMSI}}, or from a simple
#' numeric vector of masses.
#'
#' The number of possible pairs grows combinatorially and the calculation will
#' fail for >65536 values (the function uses a 16-bit integer internally).
#'
#' @param d Either an MSI data set of class msidf or msitsm, imported with the
#'          \code{\link{readMSI}} function from the mass2adduct package, or
#'          numeric vector representing m/z values of mass peaks
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
#' @seealso \code{\link{diffTabulateMSI}} to tabulate mass differences from
#'          MS imaging data object imported with \code{\link{readMSI}} function,
#'          \code{\link{diffHist}} to plot histogram of mass differences
#' @export

diffTabulate <- function(d) UseMethod("diffTabulate")
