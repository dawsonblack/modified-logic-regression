\name{plot.logregmodel}
\alias{plot.logregmodel}
\title{Plots for Logic Regression}


\description{Makes plots for an object of class \code{logregmodel}
fitted by \code{logreg}.}


\usage{\method{plot}{logregmodel}(x, pscript=FALSE, title=TRUE, nms, ...)}


\arguments{

\item{x}{object of class \code{logregmodel}, typically a part of an
object of class \code{logreg}, which is the result of the function
\code{logreg}. }

\item{pscript}{if \code{TRUE} all plots will be stored in postscript
files with distinct names. }

\item{title}{if \code{TRUE} this generates a title for some plots,
typically listing the number of trees and the model size. }

\item{nms}{names of variables. If \code{nms} is provided variable names will
be plotted, otherwise indices will be used. }

\item{...}{graphical parameters can be given as arguments to plot.}
}


\value{The fitted trees are plotted. }


\references{

Ruczinski I, Kooperberg C, LeBlanc ML (2003).  Logic Regression,
\emph{Journal of Computational and Graphical Statistics}, \bold{12}, 475-511.

Ruczinski I, Kooperberg C, LeBlanc ML (2002).  Logic Regression -
methods and software.  \emph{Proceedings of the MSRI workshop on
Nonlinear Estimation and Classification} (Eds: D. Denison, M. Hansen,
C. Holmes, B. Mallick, B. Yu), Springer: New York, 333-344.

Selected chapters from the dissertation of Ingo Ruczinski, available from
\url{https://research.fredhutch.org/content/dam/stripe/kooperberg/ingophd-logic.pdf}}


\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fredhutch.org}.
}


\seealso{
\code{\link{logreg}},
\code{\link{logregmodel}},
\code{\link{plot.logreg}},
\code{\link{logreg.testdat}}
}


\examples{
data(logreg.savefit1)
# myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 1000)
# logreg.savefit1 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21],
#                type = 2, select = 1, ntrees = 2, anneal.control = myanneal)
# plot(logreg.savefit1)
plot(logreg.savefit1$model) # does the same
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}
