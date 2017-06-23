#' Read mass spectrometry imaging data
#'
#' Import mass spectrometry imaging data in CSV format, exported from MSiReader
#' using the Intensity export function. Note that the column headers should be
#' ion masses, and that table values should be all numeric or there will be an
#' error.
#'
#' @param file filename; CSV export from MSiReader software using Intensity export
#' @param sep character; column separator in CSV file (default: ";")
#' @param remove.zeroes logical; remove columns where all values are zero? (default: TRUE)
#' @param header logical; (default: TRUE)
#' @param row.names numeric; column in import table specifying row names (default: 1)
#' @param check.names logical; (default: FALSE)
#' @param ... Parameters to be passed to read.csv()
#'
#' @return Object of class data.frame
#' @export

readMSI <- function(file, sep=";", remove.zeroes=TRUE, header=TRUE, check.names=FALSE, row.names=1, ...) {
    data <- read.csv (file=file, sep=sep, header=TRUE, check.names=FALSE, row.names=row.names, ...)
    # Check that all values are numeric
    data.isnumeric <- sapply(data, function(x) is.numeric(x))
    if (all(data.isnumeric)) {
        if (remove.zeroes) { # Remove columns that are only zero
            data.colsums <- sapply(data, function(x) sum(x))
            index.data.colzero <- which(data.colsums==0)
            index.data.notzero <- which(data.colsums!=0)
            data.nozero <- data[index.data.notzero]
            cat ("Removing", length(index.data.colzero), "columns with only zeroes\n")
            return(data.nozero)
        } else {
            return(data)
        }
    } else {
        # Error message if non-numeric values found in data input
        cat ("Error: Some values in the input data are non-numeric; please check the input file \n")
    }
    
}
