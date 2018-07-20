#' Remove peaks from a peaklist within specified mass accuracy
#'
#' @param input numeric; Peaklist to be subtracted from
#' @param vec numeric; List of peaks to be removed from v1
#' @param ppm1 numeric; Mass accuracy for peaklist v1, in ppm
#' @param ppm2 numeric; Mass accuracy for peaklist v2, in ppm
#' @param index logical; Report indices of peaks to be retained if TRUE, else
#'        report values
#' @return indices of vector v1 to be retained if \code{index=TRUE}, otherwise
#'         the values of vector v1 corresponding to those indices.
#' @export

subtractPeaks.numeric <- function(input,
                                  vec,
                                  ppm1,
                                  ppm2,
                                  index=TRUE) {
  inbounds <- function(x, v, ppm1, ppm2) {
    # x is a single numeric value
    vecdiff <- abs(v - x)
    vecmin <- v[which.min(vecdiff)]
    vecerr <- sqrt((x*ppm1*1e-6)**2 + (vecmin*ppm2*1e-6)**2)
    if (min(vecdiff) <= vecerr) {
      return (FALSE) # Do not retain if there is a match
    } else {
      return (TRUE) # Retain
    }
  }
  out <- lapply(input, function(i) inbounds(i,vec,ppm1,ppm2))
  out <- unlist(out)
  if (index) {
    return(which(out))
  } else {
    return(input[which(out)])
  }
}
