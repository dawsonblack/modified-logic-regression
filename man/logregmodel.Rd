\name{logregmodel}
\alias{logregmodel}
\title{Format of class logregmodel}


\description{This help file contains a description of the format of
class logregmodel. }


\usage{logregmodel()}


\value{
An object of class logregtree has the following components:
\item{ntrees}{the number of trees in the current model.}
\item{nleaves}{the number of leaves for the fitted model.}
\item{coef}{the coefficients for this model.}
\item{score}{the score of the fitted model.}
\item{trees}{a list of \code{ntrees} objects of class
\code{logregtree}.}
}


\references{

Ruczinski I, Kooperberg C, LeBlanc ML (2003).  Logic Regression,
\emph{Journal of Computational and Graphical Statistics}, \bold{12}, 475-511.

Ruczinski I, Kooperberg C, LeBlanc ML (2002).  Logic Regression -
methods and software.  \emph{Proceedings of the MSRI workshop on
Nonlinear Estimation and Classification} (Eds: D. Denison, M. Hansen,
C. Holmes, B. Mallick, B. Yu), Springer: New York, 333-344.

Kooperberg C, Ruczinski I, LeBlanc ML, Hsu L (2001). Sequence Analysis
using Logic Regression, \emph{Genetic Epidemiology}, \bold{21},
S626-S631.

Selected chapters from the dissertation of Ingo Ruczinski, available
from \url{http://bear.fhcrc.org/~ingor/logic/documents/myphd-logic.pdf}}


\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fhcrc.org}.
}


\seealso{
\code{\link{logreg}},
\code{\link{plot.logregmodel}},
\code{\link{print.logregmodel}},
\code{\link{logregtree}}
}


\examples{
logregmodel()       # displays this help file
help(logregmodel)   # equivalent
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}

