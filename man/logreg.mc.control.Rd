\name{logreg.mc.control}
\alias{logreg.mc.control}
\title{Control for Logic Regression}


\description{Control of MCMC annealing parameters needed in
\code{logreg}.}


\usage{logreg.mc.control(nburn=1000, niter=25000, hyperpars=0, update=0,
       output=4)}


\arguments{

\item{nburn}{number of burn in MCMC iterations that are ignored 
when computing summaries}

\item{niter}{number of MCMC iterations that are used to compute
summary statistics}

\item{hyperpars}{hyperparameters. The code allows up to 10 such
parameters, but currently only one is used. In particular,
\code{log(P(size=k)/P(size=k+1))} equals \code{hyperpars[1]}, where
P is the prior on model size. Since a maximum model size (specified
in \code{logreg} is being used, \code{hyperpars[1]} can even be 
smaller than 0.}

\item{update}{every how many iterations there should be an update of
the scores. I.e. if \code{update = 1000}, a score will get printed
every 1000 iterations. So if \code{iter = 100000} iterations, there
will be 100 updates on your screen.  If \code{update = 0}, a one
line summary for each fitted model is printed.  If \code{update = -1},
there is virtually no printed output.}

\item{output}{If \code{abs(output) > 1} bivariate statistics
are gathered, if \code{abs(output) > 2} trivariate  statistics
are also gathered, otherwise only univariate statistics are gathered. If
\code{output > 0}
all fitted models are saved in a text file ``slogiclisting.tmp'',
if \code{output < 0} this does not happen.}
}

\details{
Considerations for setting \code{nburn} and \code{niter} are as for any
MCMC problem. In our experience Logic Regression mixes quickly, and
a real small \code{nburn} (1000, for example) suffices. If there are
many trees and large models \code{niter} may need to be large.

A more detailed description of the output options can be found
in the helpfile of \code{logreg}.}


\value{A list with arguments \code{nburn}, \code{niter}, \code{hyperpars},
\code{update}, and \code{output}, that can be used as the value of
the argument \code{mc.control} of \code{logreg}.}


\references{

Ruczinski I, Kooperberg C, LeBlanc ML (2003).  Logic Regression,
\emph{Journal of Computational and Graphical Statistics}, \bold{12}, 475-511.

Ruczinski I, Kooperberg C, LeBlanc ML (2002).  Logic Regression -
methods and software.  \emph{Proceedings of the MSRI workshop on
Nonlinear Estimation and Classification} (Eds: D. Denison, M. Hansen,
C. Holmes, B. Mallick, B. Yu), Springer: New York, 333-344.

Kooperberg C, Ruczinki I (2005). Identifying interacting SNPs using
Monte Carlo Logic Regression, Genetic Epidemiology, in press.}



\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fhcrc.org}.
}


\seealso{
\code{\link{logreg}},
\code{\link{logreg.tree.control}},
\code{\link{logreg.anneal.control}}
}


\examples{
mymccontrol <- logreg.mc.control(nburn = 500, niter = 500000, update = 25000,
hyperpars = log(2), output = -2)
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}

