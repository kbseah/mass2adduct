#' Filter peak list of MSI data
#' @export

filterPeaks.data.frame <- function(d,
                                   how=c("topX","XofTop","XofTotal"),
                                   x=NULL,
                                   index=FALSE) {
  # Check input
  if (is.null(x)) {
    stop("No value supplied for parameter x")
  }
  if (is.null(d$intensities)) {
    stop("Data frame should have column named \"intensities\"")
  }
  if (is.null(d$peaks)) {
    stop("Data frame should have column named \"peaks\"")
  }
  if (! how[1] %in% c("topX","XofTop","XofTotal")) {
    stop("Invalid filtering method")
  }
  
  # topX method
  if (how[1]=="topX") {
    if (x < 1) {
      stop("Value of parameter x must be an integer >= 1")
    }
    #out <- d[order(d$intensities,decreasing=TRUE)[1:x],]
    outidx <- order(d$intensities,decreasing=TRUE)[1:x]
  } else if (how[1]=="XofTop") {
    if (x >= 1 | x <= 0) {
      stop("Value of parameter x should be between 0 and 1")
    }
    #out <- d[which(d$intensities >= x*max(d$intensities)),]
    outidx <- which(d$intensities >= x*max(d$intensities))
  } else if (how[1]=="XofTotal") {
    cutoff <- x*sum(d$intensities)
    reorderidx <- order(d$intensities, decreasing=F)
    cumulint <- cumsum(d$intensities[reorderidx])
    #out <- d[reorderidx[which(cumulint>cutoff)],]
    outidx <- reorderidx[which(cumulint>cutoff)]
  }
  
  if (index) {
    return (outidx)
  } else {
    out <- d[outidx,]
    return(out)
  }
  
}