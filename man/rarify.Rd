\name{rarify}
\alias{rarify}

\title{
Subsample via rarification
}
\description{
Subsamples when a sample size is over a given amount.
}
\usage{
rarify(inbug, sample.ID, abund, subsiz)
}
\arguments{
  \item{inbug}{
The data frame.
}
  \item{sample.ID}{
The taxa list.
}
  \item{abund}{
The abundances.
}
  \item{subsiz}{
Subsample size.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (inbug, sample.ID, abund, subsiz) 
{
    start.time = proc.time()
    outbug <- inbug
    sampid <- unique(inbug[, sample.ID])
    nsamp <- length(sampid)
    outbug[, abund] <- 0
    for (i in 1:nsamp) {
        isamp <- sampid[i]
        flush.console()
        onesamp <- inbug[inbug[, sample.ID] == isamp, ]
        onesamp <- data.frame(onesamp, row.id = seq(1, dim(onesamp)[[1]]))
        samp.expand <- rep(x = onesamp$row.id, times = onesamp[, 
            abund])
        nbug <- length(samp.expand)
        ranvec <- runif(n = nbug)
        samp.ex2 <- samp.expand[order(ranvec)]
        if (nbug > subsiz) {
            subsamp <- samp.ex2[1:subsiz]
        }
        else {
            subsamp <- samp.ex2
        }
        subcnt <- table(subsamp)
        newsamp <- onesamp
        newsamp[, abund] <- 0
        newsamp[match(newsamp$row.id, names(subcnt), nomatch = 0) > 
            0, abund] <- as.vector(subcnt)
        outbug[outbug[, sample.ID] == isamp, abund] <- newsamp[, 
            abund]
    }
    outbug
  }
}

\keyword{ ibi }
