\name{print.logregtree}
\alias{print.logregtree}
\title{Prints Logic Regression Formula}


\description{Prints formulas for objects fitted by \code{logreg}. }

\usage{\method{print}{logregtree}(x, nms, notnms, pstyle, ...)}


\arguments{

\item{x}{object of class \code{logregtree}, typically a part of an
object of class \code{logreg}, which is the result of the function
\code{logreg}, or a matrix with five columns
(see \code{logregtree}). }

\item{nms}{names of variables. If \code{nms} is provided variable names will
be printed, otherwise 
indices will be used.}

\item{notnms}{names of complements of the variables. If 
\code{notnms} is not provided ``not'' will be added before the
variable names.}

\item{pstyle}{parenthesis style. If \code{pstyle = 1} (the default)
rules are more compact than if \code{pstyle = 2}.}

\item{...}{other options are ignored}
}


\value{
A text representation of the tree will be printed.
}


\references{

Ruczinski I, Kooperberg C, LeBlanc ML (2003).  Logic Regression,
\emph{Journal of Computational and Graphical Statistics}, \bold{12}, 475-511.

Ruczinski I, Kooperberg C, LeBlanc ML (2002).  Logic Regression -
methods and software.  \emph{Proceedings of the MSRI workshop on
Nonlinear Estimation and Classification} (Eds: D. Denison, M. Hansen,
C. Holmes, B. Mallick, B. Yu), Springer: New York, 333-344.}

\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fhcrc.org}.
}


\seealso{
\code{\link{logreg}}, 
\code{\link{logregtree}}, 
\code{\link{print.logreg}}, 
\code{\link{print.logregmodel}}, 
\code{\link{logreg.testdat}} 
}


\examples{
data(logreg.savefit1)
#
# myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 1000)
# logreg.savefit1 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21],
#                type = 2, select = 1, ntrees = 2, anneal.control = myanneal)
print(logreg.savefit1$model$trees[[1]]) 
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}

