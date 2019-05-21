#' Find closest-matching known adducts for a set of mass differences or histogram thereof
#'
#' For use with internal dataset "adducts" and other similar data
#'
#' @param x Object of class massdiff from function \code{\link{massdiff}} or
#'        histogram from running \code{\link{hist}} on a massdiff object
#' @param add data.frame; Table of known chemical transfomrations and their
#'        respective masses. (default: built-in data set "adducts")
#' @param ... Other options to pass to underlying methods
#'
#' @return data.frame of mass differences given to argument diff, with an
#'         additional column reporting the closest matches
#'
#' @seealso \code{\link{adductMatch.massdiffhist}} method for histogram objects
#' @seealso \code{\link{adductMatch.massdiff}} method for massdiff objects
#' @seealso \code{\link{topAdducts}} to find the closest-matching adducts
#'          for the mass differences with the highest counts (the converse of the
#'          current function)
#'
#' @export

adductMatch <- function(x, add=mass2adduct::adducts, ...) UseMethod("adductMatch")
