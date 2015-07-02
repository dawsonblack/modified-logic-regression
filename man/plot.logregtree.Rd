\name{plot.logregtree}
\alias{plot.logregtree}
\title{A plot of one Logic Regression tree.}


\description{Makes a plot of one Logic Regression tree, fitted by
\code{logreg}.}

\usage{\method{plot}{logregtree}(x, nms, full=TRUE, and.or.cx=1.0, leaf.sz=1.0, 
                leaf.txt.cx=1.0, coef.cx=1.0, indents=rep(0,4), coef=TRUE,
                coef.rd=4, ...)}


\arguments{

\item{x}{an object of class \code{logregtree}, or the \code{trees}
component of such an object.  Typically this object will be part of
the result of an object of class \code{logreg}, generated with
\code{select = 1} (single model fit) or \code{select = 2} (multiple
model fit). }

\item{nms}{names of variables. If nms is provided variable names will
be plotted, otherwise indices will be used. }

\item{full}{if \code{TRUE}, the tree occupies the entire window with
margins specified by \code{indents}. }

\item{and.or.cx}{character expansion (size) for the operators
and/or.}

\item{leaf.sz}{character expansion for the size of the leaves. }

\item{leaf.txt.cx}{character expansion for the text in the leaves. }

\item{coef.cx}{character expansion for the coefficient string. }

\item{indents}{indents for plot - bottom, left, top, right. }

\item{coef}{if \code{TRUE}, the coefficient of the tree is plotted.}

\item{coef.rd}{controls how many digits of the above coefficient are 
displayed.}

\item{...}{graphical parameters can be given as arguments to plot.}
}


\value{This function makes a plot of one logic tree. The character
expansion terms (\code{and.or.cx, leaf.sz, leaf.txt.cx, coef.cx}) defaults of
1.0 are chosen to generate a pretty plot of a single tree with up to
eight leaves (4 levels deep). To plot more than one tree, or trees of
different complexity, scale accordingly. }


\references{

Ruczinski I, Kooperberg C, LeBlanc ML (2003).  Logic Regression,
\emph{Journal of Computational and Graphical Statistics}, \bold{12}, 475-511.

Ruczinski I, Kooperberg C, LeBlanc ML (2002).  Logic Regression -
methods and software.  \emph{Proceedings of the MSRI workshop on
Nonlinear Estimation and Classification} (Eds: D. Denison, M. Hansen,
C. Holmes, B. Mallick, B. Yu), Springer: New York, 333-344.

Selected chapters from the dissertation of Ingo Ruczinski, available
from \url{http://kooperberg.fhcrc.org/logic/documents/ingophd-logic.pdf}}


\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fredhutch.org}.
}


\seealso{
\code{\link{logreg}},
\code{\link{frame.logreg}},
\code{\link{logreg.testdat}}
}


\examples{
data(logreg.savefit2)
# 
# myanneal2 <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 0)
# logreg.savefit2 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21],
#                type = 2, select = 2, ntrees = c(1,2), nleaves =c(1,7),
#                anneal.control = myanneal2)
for(i in 1:logreg.savefit2$nmodels) for(j in 1:logreg.savefit2$alltrees[[i]]$ntrees[1]){
   plot.logregtree(logreg.savefit2$alltrees[[i]]$trees[[j]])
   title(main=paste("model",i,"tree",j))
}
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}
