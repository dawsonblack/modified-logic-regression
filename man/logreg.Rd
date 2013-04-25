\name{logreg}
\alias{logreg}
\title{Logic Regression}


\description{Fit one or a series of Logic Regression models,  carry
out cross-validation or permutation tests for such models,
or fit Monte Carlo Logic Regression models.

Logic regression is a (generalized) regression methodology that is
primarily applied when most of the covariates in the data to be
analyzed are binary. The goal of logic regression is to find
predictors that are Boolean (logical) combinations of the original
predictors.  Currently the Logic Regression methodology has scoring
functions for linear regression (residual sum of squares), logistic
regression (deviance), classification (misclassification), 
proportional hazards models (partial likelihood),
and exponential survival models (log-likelihood). A feature of the
Logic Regression methodology is that it is easily possible to extend
the method to write ones own scoring function if you have a different
scoring function.  \code{logreg.myown} contains information on how to
do so.
}


\usage{logreg(resp, bin, sep, wgt, cens, type, select, ntrees, nleaves, 
       penalty, seed, kfold, nrep, oldfit, anneal.control, tree.control,
       mc.control)}


\arguments{

\item{resp}{vector with the response variables. Let \code{n1} be the
length of this column.}

\item{bin}{matrix or data frame with binary data.  Let \code{n2} be
the number of columns of this object.  \code{bin} should have
\code{n1} rows.}

\item{sep}{(optional) matrix or data frame that is fitted additively
in the logic regression model.  \code{sep} should have \code{n1}
rows. When exponential survival models (\code{type = 5}) are
used, the additive predictors have to be binary. When
logistic regression models (\code{type = 3}) are used \code{logreg}
is much faster when all additive predictors are binary.}

\item{wgt}{(optional) vector of length \code{n1} with case weights;
default is \code{rep(1,n1)}.}

\item{cens}{(optional) an indicator variable with censoring
indicators if \code{type} equals 4 (proportional hazards model)
or 5 (exponential survival model); default is \code{rep(1,n1)}.}

\item{type}{type of model to be fit: (1) classification, (2)
regression, (3) logistic regression, (4) proportional hazards model
(Cox regression), (5) exponential
survival model, or (0) your own scoring function. If \code{type = 0},
the code needs to be recompiled, uncompiled \code{type = 0}
results in a constant score of 0, which may be useful
to generate a sample from the prior when \code{select = 7} (Monte Carlo
Logic Regression).}

\item{select}{type of model selection to be carried out: (1) fit a
single model, (2) fit multiple models, (3) cross-validation, (4)
null-model permutation test, (5) conditional permutation test,
(6) a greedy stepwise algorithm, or (7) Monte Carlo Logic Regression
(using MCMC). See
details below.}

\item{ntrees}{number of logic trees to be fit. A single number if you
select to fit a single model (\code{select = 1}), carry out the
null-model permutation test (\code{select = 4}), carry out greedy
stepwise selection  (\code{select = 6}) or,
select using MCMC (\code{select = 7}), or a range
(e.g. \code{c(ntreeslow,ntreeshigh)}) for any of the other selection
options.  In our applications, we usually ended up with models having
between one and four trees. In general, fitting one and two trees in
the initial exploratory analysis is a good idea.}

\item{nleaves}{maximum number of leaves to be fit in all trees
combined.  A single number if you select to fit a single model
(\code{select = 1}) carry out the null-model permutation test
(\code{select = 4}), carry out greedy
stepwise selection  (\code{select = 6}) or,
select using MCMC (\code{select = 7}),
or a range (e.g. \code{c(nleaveslow,nleaveshigh)})
for any of the other selection options.  If \code{select} is 1, 4, 6, or 7,
the default is \code{-1}, which is expanded to become \code{ntrees * tree.control\$treesize
}.}

\item{penalty}{specifying the penalty parameter allows you to
penalize the score of larger models. The penalty takes the form
\code{penalty} times the number of leaves in the model.
(For some score functions, we compute
average scores: the penalty is naturally
adjusted.)  Thus \code{penalty = 2} is somewhat comparable to
AIC. Note, however, that there is no relation between the size of the
model and the number of parameters (dimension) of the model, as is
usual the case for AIC like penalties. 
\code{penalty} is only relevant when \code{select = 1}.}

\item{seed}{a seed for the random number generator.  The random seed
is taken to be \code{abs(seed)}. \ For
the cross-validation version, if  \code{seed < 0} the sequence of
the cases is not permuted for division in training-test sets. This
is useful if you already permuted the sequence, and wish to compare
results with other approaches, or if there is a relation between the
sequence of the cases, for example for a matched case-control study.}

\item{kfold}{the number of groups the cases are randomly assigned
to. In turn, the model is trained on \code{(kfold - 1)} of those
groups, and scored on the group left out. Common choices
are \code{kfold = 5} and \code{kfold = 10}. Only
relevant for cross-validation (\code{select = 3}).  }

\item{nrep}{the number of runs on permuted data for each model size.
We recommend first running this program with a small number of
repetitions (e.g. 10 or 25) before sending off a big job.  Only
relevant for the null-model test (\code{select = 4}) or the
permutation test (\code{select = 5}).}

\item{oldfit}{object of class \code{logreg}, typically the result of a
previous call to \code{logreg}.  All options that are not specified
default to the value used in \code{oldfit}.
For the permutation test (\code{select = 5}) an
\code{oldfit} object obtained with \code{select = 2} (fit multiple
models) is mandatory, as the best models of each size need to be in
place.}

\item{anneal.control}{simulated annealing parameters - best set using
the function \code{logreg.anneal.control}.}

\item{tree.control}{several secondary parameters - best set using the
function \code{logreg.tree.control}.}

\item{mc.control}{Markov chain Monte Carlo parameters - best set using the
function \code{logreg.mc.control}.}
}


\value{
An object of the class \code{logreg}. This contains virtually all
input parameters, and in addition\cr

If \code{select = 1:}

an object of class \code{logregmodel}: the Logic
Regression model. This model contains a list of \code{ntrees} objects
of class \code{logregtree}.\cr

If \code{select = 2} or \code{select = 6}:

\code{nmodels:} 
the number of models fitted.\cr
\code{allscores:} 
a matrix with 3 columns, containing the scores of all models. Column 1
contains the score, column 2 the number of leaves and column 3 the
number of trees.\cr
\code{alltrees:} 
a list with \code{nmodels} objects of class \code{logregmodel}.\cr

If \code{select = 3:}

\code{cvscores:}
a matrix with the results of cross validation. The \code{train.ave}
and \code{test.ave} columns for train and test contain running
averages of the scores for individual validation sets. As such these
scores are of most interest for the rows where \code{k=kfold}.\cr

If \code{select = 4:}

\code{nullscore:}
score of the null-model.\cr
\code{bestscore:} 
score of the best model.\cr
\code{randscores:}
scores of the permutations; vector of length \code{nrep}.\cr

If \code{select = 5:}

\code{bestscore:} 
score of the best model.\cr
\code{randscores:}
scores of the permutations; each column corresponds to one model
size.\cr

If \code{select = 7:}

\code{size:} a matrix with two columns, indicating which size models
were fit how often.\cr
\code{single:} a vector with as many elements as there are binary
predictors. \code{single[i]} shows how often predictor  \code{i} is
in any of the MCMC models. Note that when a predictor is twice in the
same model, it is only counted once, thus, in particular,
\code{sum(size[,1]*size[,2]} will typically be slightly larger
than \code{sum(single)}.\cr
\code{double:} square matrix with as size the number of binary predictors.
\code{double[i,j]} shows how often predictors \code{i} and \code{j} are
in the same tree of the same MCMC model if
\code{i>j}, if \code{i<=j} \code{double[i,j]} equals zero. Note that
for models
with several logic trees two predictors can both be in the model but not
be in the same tree.\cr
\code{triple:} square 3D array with as size the number of binary predictors.
See \code{double}, but here \code{triple[i,j,k]} shows how often 
three predictors are jointly in one logic tree.\cr
In addition, the file
\code{slogiclisting.tmp} 
in the current working directory can be
created. This file contains a compact listing of all models visited. Column 1:
proportional to the log posterior probability of the model; column 2: score
(log-likelihood); column 3: how often was this model visited, column 4 through 3 + maximum number
of leaves: summary of the first tree, if there are two trees,
column 4 + maximum number of leaves through 3 + twice the maximum number
of leaves contains the second tree, and so on.\cr
In this compact notation, leaves are in the same sequence as the rows in
a \code{logregtree} object; a zero means that the leave is empty, a
1000 means an ``and'' and a 2000 an ``or'', any other positive number
means a predictor and a negative number means ``not'' that predictor.\cr
The \code{mc.control} element \code{output} can be used to surppress the
creation of
\code{double}, \code{triple}, and/or \code{slogiclisting.tmp}. There doesn't
seem to be too much use in surppressing \code{double}. Surpressing
\code{triple} speeds up computations a bit (in particular on
machines with limited memory when there are many binary predictors), and reduces the size of
both the code and the object, surppressing \code{slogiclisting.tmp}
saves the creation of a possibly very large file, which can slow
down the code considerably. See \code{logreg.mc.control} for details.
}


\details{
Logic Regression is an adaptive regression methodology that attempts
to construct predictors as Boolean combinations of binary covariates.

In most regression problems a model is developed that only relates the
main effects (the predictors or transformations thereof) to the
response.  Although interactions between predictors are considered
sometimes as well, those interactions are usually kept simple
(two- to three-way interactions at most).  But often, especially when
all predictors are binary, the interaction between many predictors is
what causes the differences in response. This issue often arises in
the analysis of SNP microarray data or in data mining problems.  Given
a set of binary predictors X, we try to create new, better predictors
for the response by considering combinations of those binary
predictors. For example, if the response is binary as well (which is
not required in general), we attempt to find decision rules such as
``if X1, X2, X3 and X4 are true'', or ``X5 or X6 but not X7 are
true'', then the response is more likely to be in class 0. In other
words, we try to find Boolean statements involving the binary
predictors that enhance the prediction for the response.  In more
specific terms: Let X1,...,Xk be binary predictors, and let Y be a
response variable.  We try to fit regression models of the form
g(E[Y]) = b0 + b1 L1+ ...+ bn Ln, where Lj is a Boolean expression of
the predictors X, such as Lj=[(X2 or X4c) and X7].  The above
framework includes many forms of regression, such as linear regression
(g(E[Y])=E[Y]) and logistic regression (g(E[Y])=log(E[Y]/(1-E[Y]))).
For every model type, we define a score function that reflects the
``quality'' of the model under consideration. For example, for linear
regression the score could be the residual sum of squares and for
logistic regression the score could be the deviance. We try to find
the Boolean expressions in the regression model that minimize the
scoring function associated with this model type, estimating the
parameters bj simultaneously with the Boolean expressions Lj. In
general, any type of model can be considered, as long as a scoring
function can be defined. For example, we also implemented the Cox
proportional hazards model, using the partial likelihood as the score.

Since the number of possible Logic Models we can construct for a given
set of predictors is huge, we have to rely on some search algorithms
to help us find the best scoring models.  We define the move set by a
set of standard operations such as splitting and pruning the tree
(similar to the terminology introduced by Breiman et al (1984) for
CART).  We investigated two types of algorithms: a greedy and a
simulated annealing algorithm. While the greedy algorithm is very
fast, it does not always find a good scoring model. The simulated
annealing algorithm usually does, but computationally it is more
expensive.  Since we have to be certain to find good scoring models,
we usually carry out simulated annealing for our case studies.
However, as usual, the best scoring model generally over-fits the
data, and methods to separate signal and noise are needed.  To assess
the over-fitting of large models, we offer the option to fit a model
of a specific size. For the model selection itself we developed and
implemented permutation tests and tests using cross-validation. If
sufficient data is available, an analysis using a training and a test
set can also be carried out.  These tests are rather complicated, so
we will not go into detail here and refer you to Ruczinski I,
Kooperberg C, LeBlanc ML (2003), cited below. 

There are two alternatives to the simulated annealing algorithm. One
is a stepwise greedy selection of models. This is when setting
\code{select = 6}, and yields a sequence of models from size 1 through
a maximum size. At each time among all the models that are one larger
than the current model the best model is selected, yielding a sequence
of models of different sizes. Usually these models are not the best 
possible, and, if the simulated annealing chain is long enough, you
should expect that the models selected using \code{select = 2} are better.

The second alternative is to run a Markov Chain Monte Carlo (MCMC) algorithm.
This is what is done in Monte Carlo Logic Regression. The algorithm used
is a reversible jump MCMC algorithm, due to Green (1995). Other than
the length of the Markov chain, the only parameter that needs to be set
is a parameter for the geometric prior on model size. Other than in many
MCMC problems, the goal in Monte Carlo Logic Regression is not to yield
one single best predicting model, but rather to provide summaries of
all models. These are exactly the elements that are shown above as
the output when \code{select = 7}.


\bold{MONITORING}

The help file for \code{logreg.anneal.control}, contains more
information on how to monitor the simulated annealing optimization for
logreg. Here is some general information.

\bold{Find the best scoring model of any size} \code{(select = 1)}

During the iterations the following information is printed out:

\tabular{rrrrrrr}{
log-temp \tab current score \tab  best score \tab  acc / \tab rej / \tab sing  \tab current parameters  \cr
  -1.000 \tab  1.494 \tab  1.494 \tab   0(  0) \tab    0 \tab  0 \tab 2.88   -1.99    0.00\cr
  -1.120 \tab  1.150 \tab  1.043 \tab 655( 54) \tab  220 \tab 71 \tab 3.63    0.15   -1.82\cr
  -1.240 \tab  1.226 \tab  1.043 \tab 555( 49) \tab  316 \tab 80 \tab 3.83    0.05   -1.71\cr
   ... \tab  \tab  \tab  \tab  \tab  \tab \cr
  -2.320 \tab  0.988 \tab  0.980 \tab 147( 36) \tab  759 \tab 58 \tab 3.00   -2.11    1.11\cr
  -2.440 \tab  0.982 \tab  0.980 \tab  25( 31) \tab  884 \tab 60 \tab 2.89   -2.12    1.24\cr
  -2.560 \tab  0.988 \tab  0.979 \tab  35( 61) \tab  850 \tab 51 \tab 3.00   -2.11    1.11\cr
   ... \tab  \tab  \tab  \tab  \tab  \tab \cr
  -3.760 \tab  0.964 \tab  0.964 \tab   2( 22) \tab  961 \tab 15 \tab 2.57   -2.15    1.55\cr
  -3.880 \tab  0.964 \tab  0.964 \tab   0( 17) \tab  961 \tab 22 \tab 2.57   -2.15    1.55\cr
  -4.000 \tab  0.964 \tab  0.964 \tab   0( 13) \tab  970 \tab 17 \tab 2.57   -2.15    1.55\cr
}

\emph{log-temp:}\cr 
logarithm (base 10) of the temperature at the last iteration before the print out.\cr
\emph{current score:}\cr 
the score after the last iterations.\cr
\emph{best score:}\cr 
the single lowest score seen at any iteration.\cr
\emph{acc:}\cr 
the number of proposed moves that were accepted since the last print
out for which the model changed, within parenthesis, the number of
those that were identical in score to the move before acceptance.\cr
\emph{rej:}\cr 
the number of proposed moves that gave numerically acceptable results,
but were rejected by the simulated annealing algorithm since the last
print out.\cr
\emph{sing:}\cr 
the number of proposed moves that were rejected because they gave
numerically unacceptable results, for example because they yielded a
singular system.\cr
\emph{current parameters:}\cr 
the values of the coefficients (first for the intercept, then for the
linear (separate) components, then for the logic trees).

This information can be used to judge the convergence of the simulated
annealing algorithm, as described in the help file of
\code{logreg.anneal.control}.  Typically we want (i) the number of
acceptances to be high in the beginning, (ii) the number of
acceptances with different scores to be low at the end, and (iii) the
number of iterations when the fraction of acceptances is moderate to
be as large as possible.

\bold{Find the best scoring models for various sizes}
\code{(select = 2)}

During the iterations the same information as for \emph{find the best
scoring model of any size}, for each size model considered.

\bold{Carry out cross-validation for model selection}
\code{(select = 3)}

Information about the simulated annealing as described above can be
printed out. Otherwise, during the cross-validation process
information like

\code{
Step  5 of 10 [ 2 trees; 4 leaves] The CV score is   1.127  1.120  1.052  1.122
}

The first of the four scores is the \emph{training}-set score on the
current validation sample, the second score is the average of all the
\emph{training}-set scores that have been processed for this model
size, the third is the \emph{test}-set score on the current
validation sample, and the fourth score is the average of all the
\emph{test}-set scores that have been processed for this model size.
Typically we would prefer the model with the lowest test-set score
average over all cross-validation steps.

\bold{Carry out a permutation test to check for signal in the data}
\code{(select = 4)}

Information about the simulated annealing as described above can be
printed out. Otherwise, first the score of the model of size 0
(typically only fitting an intercept) and the score of the best model
are printed out. Then during the permutation lines like

\code{
Permutation number   21  out of   50  has score       1.47777
}

are printed. Each score is the result of fitting a logic tree model,
on data where the response has been permuted.  Typically we would
believe that there is signal in the data if most permutations have
worse (higher) scores than the best model, but we may believe that
there is substantial over-fitting going on if these permutation scores
are much better (lower) than the score of the model of size 0.

\bold{Carry out a permutation test for model selection}
\code{(select = 5)}

To be able to run this option, an object of class \code{logreg} that
was run with \code{(select = 2)} needs to be in place.  Information
about the simulated annealing as described above can be printed
out. Otherwise, lines like

\code{
Permutation number    8  out of   25  has score       1.00767    model size  3 with 2 tree(s)}

are printed.  We can compare these scores to the tree of the same size
and the best tree.  If the scores are about the same as the one for
the best tree, we think that the ``true'' model size may be the one
that is being tested, while if the scores are much worse, the true
model size may be larger. The comparison with the model of the same
size suggests us again how much over-fitting may be going on.
\code{plot.logreg} generates informative histograms.

\bold{Greedy stepwise selection of Logic Regression models}
\code{(select = 6)}

The scores of the best greedy models of each size are printed. 

\bold{Monte Carlo Logic Regression}
\code{(select = 7)}

A status line is printed every so many iterations. This information
is probably not very useful, other than that it helps you figure out how
far the code is.

\bold{PARAMETERS}

As Logic Regression is written in Fortran 77 some parameters had to be
hard coded in. The default values of these parameters are

maximum sample size: 20000\cr
maximum number of columns in the input file: 1000\cr
maximum number of leaves in a logic tree: 128\cr
maximum number of logic trees: 5\cr
maximum number of separate parameters: 50\cr
maximum number of total parameters(separate + trees): 55\cr

If these parameters are not large enough (an error message will let you
know this), you need to reset them in \bold{slogic.f} and recompile.
In particular, the statements defining these parameters are 

\code{PARAMETER (LGCn1MAX   = 20000)}\cr
\code{PARAMETER (LGCn2MAX   =  1000)}\cr
\code{PARAMETER (LGCnknMAX  =   128)}\cr 
\code{PARAMETER (LGCntrMAX  =     5)}\cr 
\code{PARAMETER (LGCnsepMAX =    50)}\cr 
\code{PARAMETER (LGCbetaMAX =    55)}\cr 

The unfortunate thing is that you will have to change these parameter
definitions several times in the code. So search until you have found
them all.}


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

Kooperberg C, Ruczinki I (2005). Identifying interacting SNPs using
Monte Carlo Logic Regression, \emph{Genetic Epidemiology}, \bold{28}, 157-170.

Selected chapters from the dissertation of Ingo Ruczinski, available
from \url{http://kooperberg.fhcrc.org/logic/documents/ingophd-logic.pdf}}


\author{
Ingo Ruczinski \email{ingo@jhu.edu} and
Charles Kooperberg \email{clk@fhcrc.org}.
}


\seealso{
\code{\link{eval.logreg}},
\code{\link{frame.logreg}},
\code{\link{plot.logreg}},
\code{\link{print.logreg}},
\code{\link{predict.logreg}},
\code{\link{logregtree}},
\code{\link{plot.logregtree}},
\code{\link{print.logregtree}},
\code{\link{logregmodel}},
\code{\link{plot.logregtree}},
\code{\link{print.logregtree}},
\code{\link{logreg.myown}},
\code{\link{logreg.anneal.control}},
\code{\link{logreg.tree.control}}, 
\code{\link{logreg.mc.control}},
\code{\link{logreg.testdat}}
}


\examples{
data(logreg.savefit1,logreg.savefit2,logreg.savefit3,logreg.savefit4,
     logreg.savefit5,logreg.savefit6,logreg.savefit7,logreg.testdat)

myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 500, update = 100)
# in practie we would use 25000 iterations or far more - the use of 500 is only
# to have the examples run fast
\dontrun{myanneal <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 500)}
fit1 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21], type = 2,
               select = 1, ntrees = 2, anneal.control = myanneal)
# the best score should be in the 0.97-0.98 range
plot(fit1)
# you'll probably see X1-X4 as well as a few noise predictors
# use logreg.savefit1 for the results with 25000 iterations
 plot(logreg.savefit1)
 print(logreg.savefit1)
 z <- predict(logreg.savefit1)
 plot(z, logreg.testdat[,1]-z, xlab="fitted values", ylab="residuals")
# there are some streaks, thanks to the very discrete predictions
#
# a bit less output
 myanneal2 <- logreg.anneal.control(start = -1, end = -4, iter = 500, update = 0)
# in practie we would use 25000 iterations or more - the use of 500 is only 
# to have the examples run fast
\dontrun{myanneal2 <- logreg.anneal.control(start = -1, end = -4, iter = 25000, update = 0)}
#
# fit multiple models
 fit2 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21], type = 2,
               select = 2, ntrees = c(1,2), nleaves =c(1,7), anneal.control = myanneal2)
# equivalent
 fit2 <- logreg(select = 2, ntrees = c(1,2), nleaves =c(1,7), oldfit = fit1,
                anneal.control = myanneal2)
 plot(fit2)
# use logreg.savefit2 for the results with 25000 iterations
 plot(logreg.savefit2)
 print(logreg.savefit2)
# After an initial steep decline, the scores only get slightly better
# for models with more than four leaves and two trees. 
#
# cross validation
 fit3 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21], type = 2,
                select = 3, ntrees = c(1,2), nleaves=c(1,7), anneal.control = myanneal2)
# equivalent
 fit3 <- logreg(select = 3, oldfit = fit2)
 plot(fit3)
# use logreg.savefit3 for the results with 25000 iterations
 plot(logreg.savefit3)
# 4 leaves, 2 trees should top
# null model test
 fit4 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21], type = 2,
                select = 4, ntrees = 2, anneal.control = myanneal2)
# equivalent
 fit4 <- logreg(select = 4, anneal.control = myanneal2, oldfit = fit1)
 plot(fit4)
# use logreg.savefit4 for the results with 25000 iterations
plot(logreg.savefit4)
# A histogram of the 25 scores obtained from the permutation test. Also shown
# are the scores for the best scoring model with one logic tree, and the null
# model (no tree). Since the permutation scores are not even close to the score
# of the best model with one tree (fit on the original data), there is overwhelming
# evidence against the null hypothesis that there was no signal in the data. 
 fit5 <- logreg(resp = logreg.testdat[,1], bin=logreg.testdat[, 2:21], type = 2,
            select = 5, ntrees = c(1,2), nleaves=c(1,7), anneal.control = myanneal2,
            nrep = 10, oldfit = fit2)
# equivalent
 fit5 <- logreg(select = 5, nrep = 10, oldfit = fit2)
 plot(fit5)
# use logreg.savefit5 for the results with 25000 iterations and 25 permutations
 plot(logreg.savefit5)
# The permutation scores improve until we condition on a model with two trees and
# four leaves, and then do not change very much anymore. This indicates that the
# best model has indeed four leaves.
#
# greedy selection
 fit6 <- logreg(select = 6, ntrees = 2, nleaves =c(1,12), oldfit = fit1)
 plot(fit6)
# use logreg.savefit6 for the results with 25000 iterations
 plot(logreg.savefit6)
#
# Monte Carlo Logic Regression
fit7 <- logreg(select = 7, oldfit = fit1, mc.control=
               logreg.mc.control(nburn=500, niter=2500, hyperpars=log(2)))
# we need many more iterations for reasonable results
\dontrun{logreg.savefit7 <- logreg(select = 7, oldfit = fit1, mc.control=
               logreg.mc.control(nburn=1000, niter=100000, hyperpars=log(2)))}
#
plot(fit7)
# use logreg.savefit7 for the results with 25000 iterations
plot(logreg.savefit7)

}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}
