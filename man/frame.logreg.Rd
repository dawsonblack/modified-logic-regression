\name{frame.logreg}
\alias{frame.logreg}
\title{Constructs a data frame for one or more Logic Regression
models}


\description{Evaluates all components of one or more Logic Regression
models fitted by a single call to \code{logreg}.}


\usage{frame.logreg(fit, msz, ntr, newbin, newresp, newsep, newcens, newweight)}


\arguments{

\item{fit}{object of class \code{logreg}, that resulted from applying
the function \code{logreg} with 
\code{select = 1} (single model fit),
\code{select = 2} (multiple model fit), or
\code{select = 6} (greedy stepwise fit).}

\item{msz}{if \code{frame.logreg} is executed on an object of class
\code{logreg}, that resulted from applying the function \code{logreg}
with \code{select = 2} (multiple model fit)
or \code{select = 6} (greedy stepwise fit) all logic trees for all
fitted models are returned. To restrict the model size and the number
of trees to some models, specify \code{msz} and \code{ntr}
(for \code{select = 2}) or just \code{msz} (for \code{select = 6}).}

\item{ntr}{see \code{msz}.}

\item{newbin}{binary predictors to evaluate the logic trees at.  If
\code{newbin} is omitted, the original (training) data is used.}

\item{newresp}{the response. If \code{newbin} is omitted, the
original (training) response is used. If \code{newbin} is specified
and \code{newresp} is omitted, the resulting
data frame will not have a
response column.}

\item{newsep}{separate (linear) predictors. If \code{newbin} is
omitted, the original (training) predictors are used, even if
\code{newsep} is specified.}

\item{newweight}{case weights. If \code{newbin} is omitted, the
original (training) weights are used. If \code{newbin} is specified
and \code{newweight} is omitted, the weights are taken to be 1.}

\item{newcens}{censoring indicator. For proportional hazards models
and exponential survival models
only.  If \code{newbin} is omitted, the original (training) censoring
indicators are used. If \code{newbin} is specified and \code{newcens}
is omitted, the censoring indicators are taken to be 1.}
}


\value{A data frame. The first column is the response, later columns
are weights, censoring indicator, separate predictors (all of which
are only provided if they are relevant) and all logic trees. Column
names should be transparent.}


\details{This function calls \code{eval.logreg}. }


\references{

Ruczinski I, Kooperberg C, LeBlanc ML (2003).  Logic Regression,
\emph{Journal of Computational and Graphical Statistics}, \bold{12}, 475-511.

Ruczinski I, Kooperberg C, LeBlanc ML (2002).  Logic Regression -
methods and software.  \emph{Proceedings of the MSRI workshop on
Nonlinear Estimation and Classification} (Eds: D. Denison, M. Hansen,
C. Holmes, B. Mallick, B. Yu), Springer: New York, 333-344.
}


\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fhcrc.org}}


\seealso{
\code{\link{logreg}},
\code{\link{eval.logreg}},
\code{\link{predict.logreg}},
\code{\link{logreg.testdat}}
}


\examples{
data(logreg.savefit1,logreg.savefit2,logreg.savefit6)
#
# fit a single mode
# myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 1000)
# logreg.savefit1 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21], 
#               type = 2, select = 1, ntrees = 2, anneal.control = myanneal)
frame1 <- frame.logreg(logreg.savefit1)
#
# a complete sequence
# myanneal2 <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 0)
# logreg.savefit2 <- logreg(select = 2, ntrees = c(1,2), nleaves =c(1,7), 
#               oldfit = logreg.savefit1, anneal.control = myanneal2)
frame2 <- frame.logreg(logreg.savefit2)
#
# a greedy sequence
# logreg.savefit6 <- logreg(select = 6, ntrees = 2, nleaves =c(1,12), oldfit = logreg.savefit1)
frame6  <- frame.logreg(logreg.savefit6, msz = 3:5) # restrict the size

}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}
