\name{plot.logreg}
\alias{plot.logreg}
\title{Plots for Logic Regression}


\description{Makes plots for objects fitted by \code{logreg}. }

\usage{plot.logreg(x, pscript=FALSE, title=TRUE, ...)}


\arguments{

\item{x}{object of class \code{logreg}, typically the result of the
function \code{logreg}. }

\item{pscript}{if \code{TRUE} all plots will be stored in postscript
files with distinct names.}

\item{title}{if \code{TRUE} this generates a title for some plots,
typically listing the number of trees and the model size. } 

\item{...}{graphical parameters can be given as arguments to plot.}
}


\value{
The type of the plots generated depends on the value of
\code{x$select}.

if \code{select = 1} the fitted trees for the results of
\emph{find the best scoring model of any size} are plotted;

if \code{select = 2} or \code{select = 6} the fitted trees are plotted, as well as a
graph of the scores for various model sizes versus model size for the
results of \emph{find the best scoring models for various sizes},
or \emph{fit a sequence of logic regression models using a stepwise
greedy algorithm};

if \code{select = 3} training and test set scores as a function
of model size for the results of \emph{carry out cross-validation for
model selection} are plotted;

if \code{select = 4} a histogram of the permutation scores with
various important values highlighted for the results of \emph{carry
out a permutation test to check for signal in the data} are plotted;

if \code{select = 5} a series of histograms of the permutation
scores with various important values highlighted for the results of
\emph{carry out a permutation test for model selection} are plotted;

if  \code{select = 7} a histogram of the size frequency of the fitted models,
a histogram of the frequency that predictors are marginally in the model,
and (if that information was collected) an image plot
for the observed frequency that predictors were jointly in the model, and
an image plot of the observed/expected ratio of that joint frequency.
See Kooperberg and Ruczinski (2004) for the definition of the expected frequency.
}


\references{

Ruczinski I, Kooperberg C, LeBlanc ML (2003).  Logic Regression,
\emph{Journal of Computational and Graphical Statistics}, \bold{12}, 475-511.

Ruczinski I, Kooperberg C, LeBlanc ML (2002).  Logic Regression -
methods and software.  \emph{Proceedings of the MSRI workshop on
Nonlinear Estimation and Classification} (Eds: D. Denison, M. Hansen,
C. Holmes, B. Mallick, B. Yu), Springer: New York, 333-344.

Kooperberg C, Ruczinki I (2005). Identifying interacting SNPs using
Monte Carlo Logic Regression, Genetic Epidemiology, in press.

Selected chapters from the dissertation of Ingo Ruczinski, available
from \url{http://bear.fhcrc.org/~ingor/logic/documents/myphd-logic.pdf}}


\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fhcrc.org}.
}


\seealso{
\code{\link{logreg}}, 
\code{\link{plot.logregmodel}}, 
\code{\link{plot.logregtree}}, 
\code{\link{logreg.testdat}} 
}


\examples{
data(logreg.savefit1,logreg.savefit2,logreg.savefit3,logreg.savefit4,
     logreg.savefit5,logreg.savefit6,logreg.savefit7)
#
# fit a single model
# myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 1000)
# logreg.savefit1 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21],
#                type = 2, select = 1, ntrees = 2, anneal.control = myanneal)
# the best score should be in the 0.96-0.98 range
plot(logreg.savefit1)
#
# fit multiple models
# myanneal2 <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 0)
# logreg.savefit2 <- logreg(select = 2, ntrees = c(1,2), nleaves =c(1,7), 
#                oldfit = logreg.savefit1, anneal.control = myanneal2)
plot(logreg.savefit2)
# After an initial steep decline, the scores only get slightly better
# for models with more than four leaves and two trees. 
#
# cross validation
# logreg.savefit3 <- logreg(select = 3, oldfit = logreg.savefit2)
plot(logreg.savefit3)
# 4 leaves, 2 trees should give the best test set score
#
# null model test
# logreg.savefit4 <- logreg(select = 4, anneal.control = myanneal2, oldfit = logreg.savefit1)
plot(logreg.savefit4)
# A histogram of the 25 scores obtained from the permutation test. Also shown
# are the scores for the best scoring model with one logic tree, and the null
# model (no tree). As the permutation scores are not even close to the score
# of the best model with one tree (fit on the original data), there is strong
# evidence against the null hypothesis that there was no signal in the data. 
#
# Permutation tests
# logreg.savefit5 <- logreg(select = 5, oldfit = logreg.savefit2)
plot(logreg.savefit5)
# The permutation scores improve until we condition on a model with two
# trees and four leaves, and then do not change very much anymore. This 
# indicates that the best model has indeed four leaves.
#
# a greedy sequence
# logreg.savefit6 <- logreg(select = 6, ntrees = 2, nleaves =c(1,12), oldfit = logreg.savefit1)
plot(logreg.savefit6)
# Monte Carlo Logic Regression
# logreg.savefit7 <- logreg(select = 7, oldfit = logreg.savefit1, mc.control=
#                logreg.mc.control(nburn=1000, niter=100000, hyperpars=log(2)))
plot(logreg.savefit7)
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}

