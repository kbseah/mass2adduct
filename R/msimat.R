#' Read mass spectrometry imaging data
#'
#' Import mass spectrometry imaging data in either CSV or triplet format.
#' CSV file is for example exported from MSiReader using the Intensity export
#' function or from other software. For CSV file, column headers should be
#' ion masses, and table values should be all numeric or there will be an
#' error.
#'
#' Alternatively the data can be in triplet form, with separate lists of row,
#' column, and intensity values. Two additional lists hold the peak mass values
#' and the pixel spot names. This is useful for importing very large but sparse
#' data matrices. If a filename is specified to option csv, then options row,
#' cols, vals, peaks, and spots are ignored.
#'
#' @param csv filename; CSV-formatted table of intensities arranged by peak masses
#'        (columns, with column headers) and pixel spots (rows, with row names)
#' @param sep character; column separator in CSV file (default: ",")
#' @param remove.zeroes logical; remove columns where all values are zero? (default: TRUE)
#' @param rows filename; File of row numbers, 0-based (see note on triplet
#'        representation above), one per line
#' @param cols filename; File of column numbers, 0-based, one per line
#' @param vals filename; File of intensity values, one per line
#' @param peaks filename; File of peak mass values, one per line. Order in list
#'        should correspond to the column numbers given to argument cols
#' @param spots filename; File of pixel spot names, one per line. Order in list
#'        should correspond to the row numbers given to argument rows
#' @param header logical; (default: TRUE)
#' @param row.names numeric; column in import table specifying row names (default: 1)
#' @param check.names logical; (default: FALSE)
#' @param reformat.large.csv logical; Convert ##############################################################################################################
#' @param ... Parameters to be passed to read.csv()
#'
#' @return Object of class msimat
#'
#' @examples
#' d <- msimat(csv=system.file("extdata","msi.csv",package="mass2adduct"),sep=";")
#' print(d) # Summary of contents
#' plot(d) # Mass spectrum of total intensity per peak
#'
#' @export

msimat <- function(csv=NULL,
                   sep=",",
                   rows=NULL,
                   cols=NULL,
                   vals=NULL,
                   peaks=NULL,
                   spots=NULL,
                   remove.zeroes=TRUE,
                   header=TRUE,
                   check.names=FALSE,
                   row.names=1,
                   reformat.large.csv=FALSE,
                   ...) {
  if (is.null(csv)) {
    if (is.null(rows) | is.null(cols) | is.null(vals) | is.null(peaks) | is.null(spots)) {
      stop("Required input not specified")
    }
    out <- msimatreadtriplet(rows=rows,
                             cols=cols,
                             vals=vals,
                             peaks=peaks,
                             spots=spots)
    return(out)
  }
  else {
    out <- msimatreadcsv(csv=csv,
                         sep=sep,
                         remove.zeroes=remove.zeroes,
                         header=header,
                         check.names=check.names,
                         row.names=row.names,
                         ...)
    return(out)
  }
}
