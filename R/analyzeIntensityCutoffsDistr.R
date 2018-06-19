#' Calculate effects of applying intensity-level cutoffs on mass spectrum
#'
#' Given a table of total intensities per mass peak, produced by function
#' \code{\link{sumPeakIntensities}}, calculate the number of peaks retained
#' if a minimum intensity cutoff is applied (either as absolute value or as
#' percentage of the maximum intensity value observed). To do the inverse, i.e.
#' find the corresponding intensity value for a desired number of peaks, use
#' the \code{\link{quantile}} function.
#'
#' @param df data.frame of three columns: peaks, counts, intensities, produced
#'        by function \code{\link{sumPeakIntensities}} or by external script
#'        msicsvprofiler.pl
#' @param by string; Apply cutoff by intensities or counts (default: intensities)
#' @param value numeric; Minimum intensity/counts cutoff to apply 
#' @param pc numeric; Percent intensity/counts cutoff to apply (ignored if value
#'        supplied) (Default: 1)
#' @param plot logical; Plot empirical CDF?
#' @param report.peaks logical; Return peaks that are above cutoff? (default: no)
#'
#' @return data.frame summarizing the number of peaks above the specified cutoff
#'         and a plot of the distribution with cutoff overlaid as a vertical
#'         line. If option \code{report.peaks} is used, then the a vector of
#'         peaks above the cutoff is returned instead.
#'
#' @seealso \code{\link{analyzeIntensityCutoffsCumul}} does a similar analysis,
#'          but orders the peaks by intensity and takes the cumulative sum.
#' @seealso \code{\link{sumPeakIntensities}} function to generate the data.frame
#'          of total intensities and pixel counts per peak from a msimat object.
#' @export

analyzeIntensityCutoffsDistr <- function(df,
                                         by=c("intensities","counts"),
                                         value=NULL,
                                         pc=1,
                                         plot=TRUE,
                                         report.peaks=FALSE
                                         ) {
    if (class(df) != "data.frame") {
        stop("Input parameter df must be a data.frame")
    } 
    numpeaks <- length(df[["peaks"]])
    if (!is.null(value)) {
        cutoff <- value
        pc <- value / max(df[[by[1]]])
    } else if (pc <= 100 & pc >= 0) {
        cutoff <- (pc/100) * max(df[[by[1]]])
    }
    peakquant <- ecdf(df[[by[1]]])(cutoff)
    peaksabove <- numpeaks * (1 - peakquant)
    out <- data.frame(numpeaks,
                      cutoff,
                      pc,
                      peakquant,
                      peaksabove
                      )
    names(out) <- c("Total peaks",
                    "Cutoff value",
                    "Cutoff percent",
                    "Quantile",
                    "Peaks above cutoff")
    if (plot) {
        plot(ecdf(x=df[[by[1]]]),
             main=paste(c("Empirical CDF by",by[1]),collapse=" "),
             xlab=by[1]
            )
        abline(v=cutoff,col="grey")
    }
    if (report.peaks) {
        peaksshortlist <- df[["peaks"]][which(df[[by[1]]] > cutoff)]
        print(out)
        return(peaksshortlist)
    } else {
        return(out)
    }
}
