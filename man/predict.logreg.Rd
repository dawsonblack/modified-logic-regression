\name{predict.logreg}
\alias{predict.logreg}
\title{Predicted values Logic Regression}

\description{Computes predicted values for one or more Logic
Regression models that were fitted by a single call to \code{logreg}.}


\usage{\method{predict}{logreg}(object, msz, ntr, newbin, newsep, ...)}


\arguments{

\item{object}{Object of class \code{logreg}, that resulted from
applying the function \code{logreg} 
with
\code{select = 1} (single model fit),
\code{select = 2} (multiple model fit), or
\code{select = 6} (greedy stepwise fit).}


\item{msz}{if \code{predict.logreg} is executed on an object of class
\code{logreg}, that resulted from applying the function \code{logreg}
with \code{select = 2} (multiple model fit)
or \code{select = 6} (greedy stepwise fit) all logic trees for all
fitted models are returned. To restrict the model size and the number
of trees to some models, specify \code{msz} and \code{ntr}
(for \code{select = 2}) or just \code{msz} (for \code{select = 6}).}

\item{ntr}{see \code{msz}}

\item{newbin}{binary predictors to evaluate the logic trees at.  If
\code{newbin} is omitted, the original (training) data is used. }

\item{newsep}{separate (linear) predictors. If \code{newbin} is
omitted, the original (training) predictors are used, even if \code{
newsep} is specified. }

\item{...}{other options are ignored}
}


\value{If \code{object$select = 1}, a vector with fitted values,
otherwise a data frame with fitted values, where columns correspond to
models. }


\details{This function calls \code{frame.logreg}.}  


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
\code{\link{frame.logreg}},
\code{\link{logreg.testdat}}
}


\examples{
data(logreg.savefit1,logreg.savefit2,logreg.savefit6,logreg.testdat)
#
# myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 1000)
# logreg.savefit1 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21], type = 2,
#                select = 1, ntrees = 2, anneal.control = myanneal)
z1 <- predict(logreg.savefit1)
plot(z1, logreg.testdat[,1]-z1, xlab="fitted values", ylab="residuals")
# myanneal2 <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 0)
# logreg.savefit2 <- logreg(select = 2, nleaves =c(1,7), oldfit = logreg.savefit1,
#                anneal.control = myanneal2)
z2  <- predict(logreg.savefit2)
# logreg.savefit6 <- logreg(select = 6, ntrees = 2, nleaves =c(1,12), oldfit = logreg.savefit1)
z6 <- predict(logreg.savefit6, msz = 3:5)

}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}
