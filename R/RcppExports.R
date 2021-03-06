# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' C++ version of combn function
#'
#' @param x Numeric vector
#' @return Numeric matrix with three rows: all pairwise combinations of values
#'         in vector x (rows 1 and 2), and their arithmetic difference (row 3)
#' @export
comb2M <- function(x) {
    .Call('_mass2adduct_comb2M', PACKAGE = 'mass2adduct', x)
}

