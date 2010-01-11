\name{logregtree}
\alias{logregtree}
\title{Format of class logregtree}

\description{This help file contains a description of the format of
class logregtree. }


\usage{logregtree()}


\value{
An object of class logregtree is typically a substructure of an object
of the class \code{logregmodel}.  It will typically be the result of
using the fitting function \code{logreg}.  An object of class
logictree has the following components:
\item{whichtree}{the sequence number of the current tree within the
model.}
\item{coef}{the coefficients of this tree.}
\item{trees}{a matrix (data.frame) with five columns; see below for
the format.}}


\details{
When storing trees, we number the location of the nodes using the
following scheme (this is an example for a tree with at most 8
\emph{terminal} nodes, but the generalization should be obvious):

\tabular{ccccccccccccccc}{
\tab \tab \tab \tab \tab \tab \tab 1\tab \tab \tab \tab \tab \tab \tab \cr
\tab \tab \tab 2\tab \tab \tab \tab \tab \tab \tab \tab 3\tab \tab \tab \cr
\tab 4\tab \tab \tab \tab 5\tab \tab \tab \tab 6\tab \tab \tab \tab 7\tab \cr
8\tab \tab 9\tab \tab 10\tab \tab 11\tab \tab 
12\tab \tab 13\tab \tab 14\tab \tab 15\cr
}

Each node may or may not be present in the current tree. If it is
present, it can contain an operator (``and'' or ``or''), in which case
it has to child nodes, or it can contain a variable, in which case the
node is a terminal node. It is also possible that the node does not
exist (as the user only specifies the maximum tree size, not the tree
size that is actually fitted).

Output files have one line for each node. Each line contains 5
numbers:
\enumerate{
\item
the node number.
\item
does this node contain an ``and'' (1), an ``or'' (2), a variable (3),
or is the node empty (0).
\item
if the node contains a variable, which one is it; e.g. if this number
is 3 the node contains X3.
\item
if the node contains a variable, does it contain the regular variable
(0) or its complement (1)
\item
is the node empty (0) or not (1) (this information is redundant with
the second number)}

\bold{Example}
\tabular{ccccccccccccccc}{
\tab \tab \tab \tab \tab \tab \tab AND\tab \tab \tab \tab \tab \tab \tab \cr
\tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \cr
\tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \cr
\tab \tab \tab OR\tab \tab \tab \tab \tab \tab \tab OR\tab \tab \tab \tab \cr
\tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \cr
\tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \cr
\tab OR\tab \tab \tab \tab OR\tab \tab \tab \tab X20\tab \tab \tab \tab OR\tab \cr
\tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab \cr
X17\tab \tab X12\tab \tab X3\tab \tab 
X13c\tab \tab \tab \tab \tab \tab X2\tab \tab X1\cr
}

is represented as

\tabular{rrrrr}{
1 \tab  1 \tab   0 \tab   0 \tab   1\cr
2 \tab  2 \tab   0 \tab   0 \tab   1\cr
3 \tab  2 \tab   0 \tab   0 \tab   1\cr
4 \tab  2 \tab   0 \tab   0 \tab   1\cr
5 \tab  2 \tab   0 \tab   0 \tab   1\cr
6 \tab  3 \tab  20 \tab   0 \tab   1\cr
7 \tab  2 \tab   0 \tab   0 \tab   1\cr
8 \tab  3 \tab  17 \tab   0 \tab   1\cr
9 \tab  3 \tab  12 \tab   0 \tab   1\cr
10 \tab  3 \tab   3 \tab   0 \tab   1\cr
11 \tab  3 \tab  13 \tab   1 \tab   1\cr
12 \tab  0 \tab   0 \tab   0 \tab   0\cr
13 \tab  0 \tab   0 \tab   0 \tab   0\cr
14 \tab  3 \tab   2 \tab   0 \tab   1\cr
15 \tab  3 \tab   1 \tab   0 \tab   1\cr
}
}


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
Charles Kooperberg \email{clk@fhcrc.org}.
}


\seealso{
\code{\link{logreg}},
\code{\link{plot.logregtree}},
\code{\link{print.logregtree}},
\code{\link{logregmodel}}
}


\examples{
logregtree()       # displays this help file
help(logregtree)   # equivalent
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}
