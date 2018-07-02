#' Filter peak list of MSI data
#' @export

filterPeaks.msimat <- function(d,
                               how=c("topX","XofTop","XofTotal"),
                               x=NULL,
                               index=NULL) {
  df <- data.frame(peaks=d$peaks,
                   intensities=d$peakintensities)
  outidx <- filterPeaks(d=df,
                        how=how,
                        x=x,
                        index=TRUE)
    
  outidx <- outidx[order(d$peaks[outidx],decreasing=FALSE)]
  newpeaks <- d$peaks[outidx]
  newmat <- d$mat[,outidx]
  newpeakintensities <- d$peakintensities[outidx]
  
  out <- list(mat=newmat,
              peaks=newpeaks,
              spots=d$spots,
              peakintensities=newpeakintensities)
  class(out) <- "msimat"
  return(out)
}