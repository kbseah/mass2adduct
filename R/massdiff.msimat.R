#' @export

massdiff.msimat <- function(d) {
  # Extract peaklist from msimat object and apply massdiff.numeric
  peaklist <- d[["peaks"]]
  if (!is.numeric(peaklist)) {
    stop("Input peak list is not numeric")
  }
  output <- massdiff.numeric(d=peaklist)
  return(output)
}
