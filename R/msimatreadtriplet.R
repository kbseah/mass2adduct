#' Read triplet files to msimat object
#'
#' Internal function used by \code{\link{msimat}} to read triplet files
#'
#' @param rows filename; File of row numbers, 0-based (see note on triplet
#'        representation above), one per line
#' @param cols filename; File of column numbers, 0-based, one per line
#' @param vals filename; File of intensity values, one per line
#' @param peaks filename; File of peak mass values, one per line. Order in list
#'        should correspond to the column numbers given to argument cols
#' @param spots filename; File of pixel spot names, one per line. Order in list
#'        should correspond to the row numbers given to argument rows
#' @return Object of class msimat
#' @keywords internal

msimatreadtriplet <- function(rows,
                              cols,
                              vals,
                              peaks,
                              spots) {
  inrow <- scan(rows,what=numeric(),sep="\n")
  incol <- scan(cols,what=numeric(),sep="\n")
  inval <- scan(vals,what=numeric(),sep="\n")
  inpeaks <- scan(peaks,what=numeric(),sep="\n")
  inspots <- scan(spots,what=character(),sep="\n")
  tsm <- Matrix::sparseMatrix(i=inrow,
                              j=incol,
                              x=inval,
                              index1=FALSE,       # 0-based numbering
                              giveCsparse=FALSE,  # TsparseMatrix
                              check=TRUE)
  peakintensities <- Matrix::colSums(tsm)
  data <- list(mat=tsm,
               peaks=inpeaks,
               spots=inspots,
               peakintensities=peakintensities)
  class(data) <- "msimat"
  return(data)
}