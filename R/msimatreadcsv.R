#' Read CSV file to msimat object
#'
#' Internal function used by \code{\link{msimat}} to read CSV file
#'
#' @param csv filename; CSV-formatted table of intensities arranged by peak masses
#'        (columns, with column headers) and pixel spots (rows, with row names)
#' @param sep character; column separator in CSV file (default: ",")
#' @param remove.zeroes logical; remove columns where all values are zero? (default: TRUE)
#' @param header logical; (default: TRUE)
#' @param row.names numeric; column in import table specifying row names (default: 1)
#' @param check.names logical; (default: FALSE)
#' @param ... Parameters to be passed to read.csv()
#' @return Object of class msimat
#' @keywords internal

msimatreadcsv <- function (csv,
                           sep,
                           remove.zeroes,
                           header,
                           check.names,
                           row.names,
                           ...) {
  data <- read.csv (file=csv,
                    sep=sep,
                    header=header,
                    check.names=check.names,
                    row.names=row.names,
                    ...)
  # Check that all values are numeric
  data.isnumeric <- sapply(data, function(x) is.numeric(x))
  if (all(data.isnumeric)) {
    if (remove.zeroes) { # Remove columns that are only zero
      data.colsums <- sapply(data, function(x) sum(x))
      index.data.colzero <- which(data.colsums==0)
      index.data.notzero <- which(data.colsums!=0)
      message("Removing ", length(index.data.colzero), " columns with only zeroes")
      data <- data[index.data.notzero]
    }
    peaks <- names(data)
    spots <- row.names(data)
    peakintensities <- colSums(data)
    out <- list(mat=as.matrix(data),
                peaks=as.numeric(peaks),
                spots=spots,
                peakintensities=peakintensities)
    class(out) <- "msimat"
    return(out)
  } else {
    # Error message if non-numeric values found in data input
    stop("Some values in the input data are non-numeric; please check the input file")
  }
}