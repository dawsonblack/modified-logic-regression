\name{logreg.myown}
\alias{logreg.myown}
\title{Writing your own Logic Regression scoring function}


\description{Help file for writing your own scoring function for
\code{logreg}!}


\usage{logreg.myown()}

\details{

You can write your own scoring function for \code{logreg}!  This may
be useful if you have a model other than those which we already
programmed in.

Essentially you need to provide two routines in the file
\bold{My\_own\_scoring.f}:\cr
(i) A routine \emph{My\_own\_fitting} which fits your model: it
provides a coefficient (beta) for each of the logic trees and provides
a score of how good the model is. Low scores are good. (So add a minus
sign if your score is a log-likelihood.)\cr
(ii) A routine \emph{My\_own\_scoring} which - given the betas -
provides the score of your model. [If you don't use cross-validation,
this second routine is not needed, though some dummy routine to
satisfy the compiler should still be provided.]

After recompilation, you can fit your model using the option \code{
type = 0} in \code{logreg}.  Below we give an example for a version
of the My.own functions for conditional logistic regression which are
also provided as \bold{inst/condlogic.ff} when you downloaded the files.

\bold{PROGRAMMING DETAILS}

Below is a list of variables that are passed on. Most of them are as
you expect - response, predictors (binary ones and continuous ones),
number of cases, number of predictors.  In addition there are two
columns - \code{dcph} and \code{weight} - that can either be used to pass on an
auxiliary variable for each case (discrete for 
\code{dcph} and continuous for
\code{weight}), or even some overall auxiliary variables - as these numbers
are not used anywhere else.  If you do not need any of the variables -
just ignore them!

\code{prtr}:\cr 
the predictions of the logic trees in the current model: this is an
integer matrix of size \code{n1} times \code{ntr} - although only the
first \code{nop} columns contain useful information.
\cr
\code{rsp}:\cr
the response variable: this is a real (single precision) vector of
length \code{n1}.  
\cr
\code{dcph}: \cr
censor times: this is an integer vector of length \code{n1} this could
be used as an auxiliary (integer) vector - as it is just passed on.
(There is no check that this is a 0/1 variable, when you use your own
scoring function.)  For example, you could use this to pass on
something like cluster membership.
\cr
\code{weight}:\cr
weights for the cases this is a real vector of length \code{n1}.  this
could be used as an auxiliary (real) vector - as it is just passed on.
There is no check that these numbers are positive, when you choose
your own scoring function.
\cr
\code{ordrs}:\cr 
the order (by response size) of the cases This is an integer vector of
length \code{n1}. For the case with the smallest response this one is
1, for the second smallest 2, and so on. Ties are resolved arbitrary.
Always computed, although only used for proportional hazards
models. Use it as you wish.
\cr
\code{n1}:\cr 
the total number of cases in the data.
\cr
\code{ntr}:\cr 
the number of logic trees ALLOWED in the tree.
\cr
\code{nop}:\cr 
the number of logic trees in the CURRENT model.  The subroutines
should work if \code{nop} is 0.
\cr
\code{wh}:\cr 
the index of the tree that has been edited in the last move - i.e. the
column of \code{prtr} that has changes since the last call.
\cr
\code{nsep}:\cr 
number of variables that get fit a separate parameter The subroutines
should work if \code{nsep} is 0.
\cr
\code{seps}:\cr 
array of the above variables - this is a single precision matrix of
size \code{nsep} times \code{n1}. Note that \code{seps} and
\code{prtr} are stored in different directions.


For \emph{My\_own\_fitting} you should return:\cr
\code{betas}:\cr 
a vector of parameters of the model that you fit.  \code{betas(0)}
should be the parameter for the intercept \code{betas(1:nsep)} should
be the parameters for the continuous variables in seps
\code{betas((nsep+1):(nsep+nop))} should be the parameters for the
binary trees in prtr if you have more parameters, use \code{dcph}, or
\code{weight}; these variables will not be printed however.
\cr
\code{score}:\cr 
whatever score you assign to your model small should be good
(i.e. deviance or -log.likelihood).
\cr
\code{reject}:\cr 
an indicator whether or not to reject the proposed move *regardless*
of the score (for example when an iteration necessary to determine the
score failed to converge (0 = move is OK ; 1 = reject move) set this
one to 0 if there is no such condition.

You are allowed to change the values of dcph, and weight.

For \emph{My\_own\_scoring} additional input is:\cr
\code{betas}:
the coefficients

You should return:\cr
\code{score}:
whatever score you assign to your model small should be good
(i.e. deviance or -log.likelihood). If the model "crashes", you should
simply return a very large number.

While we try to prevent that models are singular, it is possible that
for your model a single or degenerate model is passed on for
evaluation. For My\_own\_fitting you can pass the model back with
\code{reject = 1}, for My\_own\_scoring you can pass it on with a very large
value for \code{score}.  Currently My\_own\_scoring.f contains empty frames for
the scoring functions; condlogic.ff contains an example with
conditional logistic regression.

The logic regression program is written in Fortran 77.

\bold{CONDITIONAL LOGISTIC REGRESSION}

A function for a conditional logistic regression score function is
attached as an example function on how to write your own scoring
function for Logic Regression.  Obviously you can also use it if you
have conditional logistic data.

Conditional logistic regression is common model fitting technique for
matched case-control studies, in which each case is matched with one
or more controls. (In conditional logistic regression several cases
could be matched to several controls, in the implementation provided
here only one case can be matched with each group of controls.)
Conditional logistic regression models are parameterized like regular
logistic regression models, except that the intercept is no longer
identifiable. See, for example, Breslow and Day - Volume 1 (1990,
Statistical Methods in Cancer Research, International Agency for
Research on Cancer, Lyon) for details.  Conditional logistic
regression models are most easily fit using a stratified proportional
hazards model (if there is one-to-one case-control matching it can
also be fit using logistic regression, but that method breaks down if
there is more than one control per case). Each group of a case and
controls is one stratum. All cases get an arbitrarily event time of
1.00, and all controls get a censoring time of 2.00.

In our implementation we use the response column to indicate the
matching. For all controls this column is 0, for a case it is k,
indicating that the next k records are the matched controls for the
current case. Thus, we order our cases so that each case is followed
by its controls. Cases with a negative response are put in a stratum
-1, which is not used in any computations.  This has implications for
cross-validation. See below.

In \emph{My\_own\_fitting} and \emph{My\_own\_scoring} we first allocate
various vectors (strata, index, censoring variable) that are local, as
well as some work arrays that are used by our fitting routines.  (We
need to set some of the parameters for that, see the help page of
\code{logreg} for details.)  We then define \code{idx(j)=j} for
\code{j=1,n1}, and we define the \code{strata} and \code{delta}
vectors.  We use slightly modified versions of the proportional
hazards routines that are already used otherwise in the Logic
Regression program, to include stratification. After the model is
fitted, we assign minus the partial likelihood to \code{score(1)} and
(for \emph{My\_own\_fitting}) we pass on the betas.

Recompile after replacing \emph{My\_own\_scoring.f} by
\emph{condlogic.ff}

The permutation and null model versions are not directly usable (we
could do some permutation tests, but they require more programming),
but we can use cross-validation. Obviously we should keep cases and
controls match. To that extend, we would run permutation with a
negative seed (see \code{logreg}) and we would take care ourselves
that case-control groups are in a random order, and that every block
has the same number of records. We achieve the later by adding some
records with response -1.  In particular, suppose that we have 19
pairs of case- (single) control data, and that we want to do 3-fold
cross validation. We would permute the sequence of the 19 pairs, and
add two records with response -1 after the 13th pair, and two records
with -1 at the end of the file, so that the total data file would have
42 records.
}


\references{

Ruczinski I, Kooperberg C, LeBlanc ML (2003).  Logic Regression,
\emph{Journal of Computational and Graphical Statistics}, \bold{12}, 475-511.

Ruczinski I, Kooperberg C, LeBlanc ML (2002).  Logic Regression -
methods and software.  \emph{Proceedings of the MSRI workshop on
Nonlinear Estimation and Classification} (Eds: D. Denison, M. Hansen,
C. Holmes, B. Mallick, B. Yu), Springer: New York, 333-344.

Kooperberg C, Ruczinski I, LeBlanc ML, Hsu L (2001). Sequence Analysis
using Logic Regression, \emph{Genetic Epidemiology}, \bold{21},
S626-S631.

Selected chapters from the dissertation of Ingo Ruczinski, available
from \url{http://bear.fhcrc.org/~ingor/logic/documents/myphd-logic.pdf}}



\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fhcrc.org}.
}


\seealso{
\code{\link{logreg}}
}


\examples{
logreg.myown()      # displays this help file
help(logreg.myown)  # equivalent
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}
