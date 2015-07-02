\name{logreg.savefit1}
\alias{logreg.savefit1}
\alias{logreg.savefit2}
\alias{logreg.savefit3}
\alias{logreg.savefit4}
\alias{logreg.savefit5}
\alias{logreg.savefit6}
\alias{logreg.savefit7}
\title{Sample results for Logic Regression}


\description{
The \code{logreg.savefit} objects are the results of fitting
\code{\link{logreg}} with various options. The examples in
the functions of the LogicReg packages all use far fewer
iterations than is needed. (The number of iterations was reduced
to provide quick results for bug-checking.) The number of
iterations in the \code{logreg.savefit} objects are more
reasonable (though they would still be small for larger
problems). Otherwise the arguments used to fit the
\code{logreg.savefit} objects are identical as those used
in the examples of \code{\link{logreg}}. The
\code{logreg.savefit} objects are used for examples involving
things like plotting, printing, and predicting.
}

\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fredhutch.org}}

\examples{
data(logreg.savefit1)
print(logreg.savefit1$call)
data(logreg.savefit2)
print(logreg.savefit2$call)
data(logreg.savefit3)
print(logreg.savefit3$call)
data(logreg.savefit4)
print(logreg.savefit4$call)
data(logreg.savefit5)
print(logreg.savefit5$call)
data(logreg.savefit6)
print(logreg.savefit6$call)
data(logreg.savefit7)
print(logreg.savefit7$call)
}

\keyword{datasets}


\seealso{
\code{\link{logreg}},
\code{\link{logreg.testdat}}
}
