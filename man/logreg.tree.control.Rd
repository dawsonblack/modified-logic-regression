\name{logreg.tree.control}
\alias{logreg.tree.control}
\title{Control for logreg}


\description{Control of various secondary parameters of tree shape
needed in \code{logreg}.}


\usage{logreg.tree.control(treesize=8, opers=1, minmass=0, n1)}


\arguments{

\item{treesize}{specify the maximum number of allowed leaves per
logic tree.  Allowing one leave means that the tree is (at most) a
simple predictor, two leaves allows for trees such as (X1 or X2) or
(not X3 and X4). Four, eight and sixteen leaves allow for two, three
or four levels of operators.  To be able to interpret the results, do
not choose too many leaves. Since the model selection techniques
usually trim down the trees, it is recommend to allow at least four or
eight leaves per tree. }

\item{opers}{The default is to allow both "and" and "or" operators in
the logic trees. If the interest is in logic statements in disjunctive
normal form, use only one of the two operator types. Choose 1 for
both operators, 2 for only "and" and 3 for only "or". }

\item{minmass}{specify the minimum number of cases for which any tree
needs to be 1 and for which any tree needs to be 0 to be considered as
a logic tree in the model. This is to prevent that \code{logreg}, will
select trees with, for example, 999 1s and one 0 out of 1000 cases.
The default is to take 5\% of the cases or 15, whatever is less. }

\item{n1}{if you specify the sample size \code{n1}, it is checked
that \code{minmass} is smaller than \code{n1/4}.  This option is used
by \code{logreg}, but is likely not useful for direct use. }
}


\value{A list with components \code{treesize}, \code{opers}, and
\code{minmass}, that can be used as the value of the argument
\code{tree.control} of \code{logreg}. }


\details{Missing arguments take defaults. If the argument
\code{treesize} is a list with arguments \code{treesize},
\code{opers}, and \code{minmass}, those values take precedent of
directly specified values. }


\references{

Ruczinski I, Kooperberg C, LeBlanc ML (2003).  Logic Regression,
\emph{Journal of Computational and Graphical Statistics}, \bold{12}, 475-511.

Ruczinski I, Kooperberg C, LeBlanc ML (2002).  Logic Regression -
methods and software.  \emph{Proceedings of the MSRI workshop on
Nonlinear Estimation and Classification} (Eds: D. Denison, M. Hansen,
C. Holmes, B. Mallick, B. Yu), Springer: New York, 333-344.

Selected chapters from the dissertation of Ingo Ruczinski, available
from \url{http://bear.fhcrc.org/~ingor/logic/documents/myphd-logic.pdf}}



\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fhcrc.org}.
}


\seealso{
\code{\link{logreg}},
\code{\link{logreg.anneal.control}},
\code{\link{logreg.mc.control}}}


\examples{
mytreecontrol <- logreg.tree.control(treesize = 16, minmass = 10)
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}

