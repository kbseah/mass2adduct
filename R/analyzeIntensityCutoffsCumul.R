#' Calculate effects of applying cumulative intensity-level cutoffs on mass spectrum
#'
#' Given a table of total intensities per mass peak, produced by function
#' \code{\link{sumPeakIntensities}}, calculate the number of peaks retained
#' if a minimum intensity cutoff is applied (either as absolute value or as
#' percentage of the TOTAL intensity value observed) to the cumulative sum of
#' the intensity values. This has similarities to the N50 metric applied to genome
#' assemblies.
#'
#' @param df data.frame of three columns: peaks, counts, intensities, produced
#'        by function \code{\link{sumPeakIntensities}}, or by external script
#'        msicsvprofiler.pl if matrix is too large to import into R
#' @param by string; Apply cutoff by intensities or counts (default: "intensities")
#' @param value numeric; Minimum intensity/counts cutoff to apply (default: NULL)
#' @param pc numeric; Percent intensity/counts cutoff to apply (ignored if value
#'        supplied). For example, 50% would find the top peaks arranged by
#'        intensity, whose combined total intensity accounts for 50% of the total
#'        intensity in the data (default: 1)
#' @param plot logical; Plot distribution? (default: TRUE)
#' @param log.plot logical; Plot cumulative values with log scale (default: TRUE)
#' @param report.peaks logical: Report list of peaks above the cutoff, as vector
#'        (default: FALSE)
#'
#' @return data.frame summarizing the number of peaks above the specified cutoff
#'         and a plot of the distribution with cutoff overlaid as a vertical
#'         line. If option \code{report.peaks} is used, then the a vector of
#'         peaks above the cutoff is returned instead.
#'
#' @seealso \code{\link{analyzeIntensityCutoffsDistr}} does a similar analysis
#'          but without taking cumulative sums, instead using the empirical CDF
#'          of the values.
#' @seealso \code{\link{sumPeakIntensities}} function to generate the data.frame
#'          of total intensities and pixel counts per peak from a msimat object.
#' @export

analyzeIntensityCutoffsCumul <- function(df,
                                         by=c("intensities","counts"),
                                         value=NULL,
                                         pc=1,
                                         plot=TRUE,
                                         log.plot=TRUE,
                                         report.peaks=FALSE
                                         ) {
    if (class(df) != "data.frame") {
        stop("Input parameter df must be a data.frame\n")
    }
    numpeaks <- length(df[["peaks"]])
    if (!is.null(value)) {
        cutoff <- value
        pc <- value / sum(df[[by[1]]])
    } else if (pc <= 100 & pc >= 0) {
        cutoff <- (pc/100) * sum(df[[by[1]]])
    }
    reorderidx <- order(df[[by[1]]], decreasing=F)
    cumulint <- cumsum(df[[by[1]]][reorderidx])
    peaksabove <- length(which(cumulint > cutoff))
    out <- data.frame(numpeaks,
                      max(cumulint),
                      cutoff,
                      pc,
                      peaksabove
                      )
    names(out) <- c("Total peaks",
                    "Total value",
                    "Cutoff value",
                    "Cutoff percent",
                    "Peaks above cutoff")
    if (plot) {
        if (log.plot) {
            yvals <- log10(cumulint)
            ycutoff <- log10(cutoff)
            ylabel <- "Cumulative values (log10)"
        } else {
            yvals <- cumulint
            ycutoff <- cutoff
            ylabel <- "Cumulative values"
        }
        plot(x=1:numpeaks,
             y=yvals,
             main=paste(c("Cumulative sum by ",by[1]),collapse=" "),
             xlab="Number of peaks",
             ylab=ylabel,
             type="l"
            )
        abline(h=ycutoff,col="grey")
    }
    if (report.peaks) {
        peaksshortlist <- df[["peaks"]][reorderidx[which(cumulint>cutoff)]]
        print(out)
        return(peaksshortlist)
    } else {
        return(out)
    }
}
