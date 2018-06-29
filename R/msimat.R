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
#' It is also possible to convert a CSV file into triplet form before import, by
#' calling the Perl script \code{msimunging.pl}, included in the package source,
#' either outside of R, or during the import (see description of parameter
#' \code{reformat.large.csv} below). The conversion may take a long time, and
#' hence it may be more practical to convert large CSV files first before
#' interactive data analysis.
#'
#' @param csv filename; CSV-formatted table of intensities arranged by peak
#'        masses (columns, with column headers) and pixel spots (rows, with row
#'        names)
#' @param sep character; column separator in CSV file (default: ",")
#' @param remove.zeroes logical; remove columns where all values are zero?
#'        (default: TRUE)
#' @param rows filename; File of row numbers, 0-based (see note on triplet
#'        representation above), one per line
#' @param cols filename; File of column numbers, 0-based, one per line
#' @param vals filename; File of intensity values, one per line
#' @param peaks filename; File of peak mass values, one per line. Order in list
#'        should correspond to the column numbers given to argument cols
#' @param spots filename; File of pixel spot names, one per line. Order in list
#'        should correspond to the row numbers given to argument rows
#' @param header logical; (default: TRUE)
#' @param row.names numeric; column in import table specifying row names
#'        (default: 1)
#' @param check.names logical; (default: FALSE)
#' @param reformat.large.csv filename prefix; If this is not null, then the CSV
#'        file specified to option \code{csv} will be converted to triplet
#'        format and imported as a sparse matrix. The character string passed to
#'        this option will be used as a filename prefix for the output files.
#'        This option can be used if the CSV file is too large to be imported
#'        as-is, and the majority of data in the matrix are zeroes. If not (a
#'        "dense" matrix), then the data may need to be reduced before import,
#'        e.g. by sampling a subset of the pixels, or using a more stringent
#'        filtering cutoff for the peaks.
#' @param crlf logical; Whether input CSV file uses CRLF line endings (Windows/
#'        DOS format).
#' @param script.path Path to \code{msimunging.pl} script, if 
#'        \code{reformat.large.csv} is not NULL. By default looks in the
#'        installation path for the package
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
                   reformat.large.csv=NULL,
                   crlf=FALSE,
                   script.path=system.file("msimunging.pl",package="mass2adduct"),
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
  }
  else {
    if (!is.null(reformat.large.csv)) {
      # Convert CSV to triplet format with msimunging.pl
      command <- "perl"
      command.params <- paste(script.path,
                              "--csv", csv,
                              paste(c("--delim=\"",sep,"\""),collapse=""),
                              "--outfmt triplet",
                              "--out", reformat.large.csv
                              )
      if (crlf) { # If file has CRLF line endings
        command.params <- paste(command.params, "--crlf")
      }
      if (remove.zeroes) { # Remove zeroes if necessary
        command.params <- paste(command.params, "--filterby nozero")
      }
      tripletfiles <- system2(command,
                              command.params,
                              stderr=NULL,
                              stdout=TRUE)
      if (length(tripletfiles < 5)) {
        # Throw error message if not all files are reported
        stop("Conversion from CSV to triplet format failed")
      } else {
        message("CSV file converted to triplet format, output files in current working folder")
      }
      # Import into R as msimat object
      out <- msimatreadtriplet(rows=paste(c(reformat.large.csv,"_rows.list"),collapse=""),
                               cols=paste(c(reformat.large.csv,"_cols.list"),collapse=""),
                               vals=paste(c(reformat.large.csv,"_vals.list"),collapse=""),
                               peaks=paste(c(reformat.large.csv,"_peaknames.list"),collapse=""),
                               spots=paste(c(reformat.large.csv,"_spotnames.list"),collapse=""))
    } else {
      out <- msimatreadcsv(csv=csv,
                           sep=sep,
                           remove.zeroes=remove.zeroes,
                           header=header,
                           check.names=check.names,
                           row.names=row.names,
                           ...)
    }
  }
  return(out)
}
