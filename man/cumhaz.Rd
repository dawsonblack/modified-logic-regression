\name{cumhaz}
\alias{cumhaz}
\title{Cumulative hazard transformation}


\description{Transforms survival times using the cumulative hazard function. }

\usage{cumhaz(y, d)}


\arguments{

\item{y}{vector of nonnegative survival times}

\item{d}{vector of censoring indicators, should be the same length
as \code{y}. If \code{d} is missing the data is assumed to be uncensored.}}


\value{
A vector of transformed survival times.}



\note{The primary use of doing a cumulative
hazard transformation is that after such a transformation, exponential
survival models yield results that are often very much comparable to
proportional hazards models. In our implementation of
Logic Regression, however, exponential
survival models run much faster than proportional hazards models when there
are no continuous separate covariates.
}



\references{

Ruczinski I, Kooperberg C, LeBlanc ML (2003).  Logic Regression,
\emph{Journal of Computational and Graphical Statistics}, \bold{12}, 475-511.
}


\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fhcrc.org}.
}


\seealso{
\code{\link{logreg}}
}


\examples{
data(logreg.testdat)
#
# this is not survival data, but it shows the functionality
yy <- cumhaz(exp(logreg.testdat[,1]), logreg.testdat[, 2])
# then we would use
# logreg(resp=yy, cens=logreg.testdat[,2], type=5, ...
# insted of
# logreg(resp=logreg.testdat[,1], cens=logreg.testdat[,2], type=4, ...
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}

