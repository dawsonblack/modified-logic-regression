\name{print.logreg}
\alias{print.logreg}
\title{Prints Logic Regression Output}


\description{Prints formulas for objects fitted by \code{logreg}. }

\usage{print.logreg(x, nms, notnms, pstyle, ...)}


\arguments{

\item{x}{object of class \code{logreg}, typically the result of the
function \code{logreg}. }

\item{nms}{names of variables. If \code{nms} is provided variable names will
be printted, otherwise \code{x$binnames} will be used. If that does
not exist indices will be used.}

\item{notnms}{names of complements of the variables. If 
\code{notnms} is not provided ``not'' will be added before the
variable names.}

\item{pstyle}{parenthesis style. If \code{pstyle = 1} (the default)
rules are more compact than if \code{pstyle = 2}.}

\item{...}{other options are ignored}
}


\value{
If \code{x$select} equals 1 or 2 the fitted logic rule(s)
are generated as a text string. Scores, and if
\code{x$select} equals 2 or 6 modelsizes, are also provided. 
If
\code{x$select} equals 4 or 5 a summary of the permutation test(s) is printed.
If
\code{x$select} equals 3 a summary of the cross validation is printed.
If \code{x$select} is equal to 7 an error message is generated.
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
\code{\link{print.logregmodel}}, 
\code{\link{print.logregtree}}, 
\code{\link{logreg.testdat}} 
}


\examples{
data(logreg.savefit1,logreg.savefit2,logreg.savefit3,logreg.savefit4,
     logreg.savefit5,logreg.savefit6)
#
# fit a single model
# myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 1000)
# logreg.savefit1 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21],
#                type = 2, select = 1, ntrees = 2, anneal.control = myanneal)
# the best score should be in the 0.96-0.98 range
print(logreg.savefit1)
#
# fit multiple models
# myanneal2 <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 0)
# logreg.savefit2 <- logreg(select = 2, ntrees = c(1,2), nleaves =c(1,7),
#                oldfit = logreg.savefit1, anneal.control = myanneal2)
print(logreg.savefit2)
# After an initial steep decline, the scores only get slightly better
# for models with more than four leaves and two trees.
#
# cross validation
# logreg.savefit3 <- logreg(select = 3, oldfit = logreg.savefit2)
print(logreg.savefit3)
# 4 leaves, 2 trees should give the best test set score
#
# null model test
# logreg.savefit4 <- logreg(select = 4, anneal.control = myanneal2, oldfit = logreg.savefit1)
print(logreg.savefit4)
# A summary of the permutation test
#
# Permutation tests
# logreg.savefit5 <- logreg(select = 5, oldfit = logreg.savefit2)
print(logreg.savefit5)
# A table summarizing the permutation tests
#
# a greedy sequence
# logreg.savefit6 <- logreg(select = 6, ntrees = 2, nleaves =c(1,12), oldfit = logreg.savefit1)
print(logreg.savefit6)
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}

