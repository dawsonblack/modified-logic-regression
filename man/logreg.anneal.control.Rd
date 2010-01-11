\name{logreg.anneal.control}
\alias{logreg.anneal.control}
\title{Control for Logic Regression}


\description{Control of simulated annealing parameters needed in
\code{logreg}.}


\usage{logreg.anneal.control(start=0, end=0, iter=0, earlyout=0, update=0)}


\arguments{

\item{start}{the upper temperature (on a log10 scale) in the annealing
chain. I.e. if \code{start = 3}, the annealing chain starts at
temperature 1000. The acceptance function is the usual
\code{min(1,exp(-diff(scores)/temp))}, so any temperature larger than what
would be expected as possible differences between any two models
pretty much generates a random walk in the beginning, and means that
you need to wait longer on results. A too low starting temperature
means that the chain may end up in a local optimal (rather than global
optimal) solution.  If you select both \code{start} and \code{end} the
default of 0, the program will attempt to find reasonable numbers
itself (it is known to be only moderately successful in this though).}

\item{end}{the lower temperature (on a log10 scale) in the annealing
chain. I.e. if \code{end} is -2, the annealing chain ends at
temperature 0.01. If this temperature is very low one can use the
early out possibility listed below, as otherwise the chain may run
longer than desired!}

\item{iter}{the total number of iterations in the annealing chain.
This is the total over all annealing chains, not the number of
iterations of a chain at a given temperature. If this number is too
small the chain may not find a good (the best) solution, if the chain
is too long the program may take long...}

\item{earlyout}{if the \code{end} temperature is very low, the simulated
annealing algorithm may not move any more, but one still needs to wait
on all possible moves being evaluated (and rejected)!  An early out
possibility is offered. If during consecutive five blocks of
\code{earlyout} iterations, in each block 10 or fewer moves are
accepted (for which the score changes), the program terminates.  This
is a desirable option after one is convinced the program otherwise
runs fine: it can be dangerous on the first run.}

\item{update}{every how many iterations there should be an update of
the scores. I.e. if \code{update = 1000}, a score will get printed
every 1000 iterations. So if \code{iter = 100000} iterations, there
will be 100 updates on your screen.  If you \code{update = 0}, a one
line summary for each fitted model is printed.  If \code{update = -1},
there is virtually no printed output.}
}


\value{A list with arguments \code{start}, \code{end}, \code{iter},
\code{earlyout}, and \code{update}, that can be used as the value of
the argument \code{anneal.control} of\code{logreg}.}


\details{
Missing arguments take defaults. If the argument \code{start} is a
list with arguments \code{start}, \code{end}, \code{iter},
\code{earlyout}, and \code{update}, those values take precedent of
directly specified values.

This is a rough outline how the automated simulated annealing works:
The algorithm starts running at a very high temperature, and decreases
the temperature until the acceptance ratio of moves is below a certain
threshold (in the neighborhood of 95\%). At this point we run longer
chains at fixed temperatures, and stop the search when the last "n"
consecutive moves have been rejected. If you think that the search was
either not sufficiently long or excessively long (both of which can
very well happen since it is pretty much impossible to specify default
values that are appropriate for all sorts of data and models), you can
over-write the default values.

If you want more detailed information continue reading....

These are some more detailed suggestions on how to set the parameters
for the beginning temperature, end temperature and number of
iterations for the Logic Regression routine. Note that if \code{start}
temperature and \code{end} temperature are both zero, the routine uses
its default values. The number of iterations \code{iter} is irrelevant
in this case. In our opinion, the default values are OK, but not
great, and you can usually do better if you're willing to invest time
in learning how to set the parameters.

The starting temperature is the log(10) value of \code{start} -
i.e., if \code{start} is 2
it means iterations start at a temperature of 100. The
\code{end} temperature is again the log(10) value. The number of iterations
are equidistant on a log-scale.

Considerations in setting these parameters.....

1) \code{start} temperature. If this is too high you're "wasting time", as
the algorithm is effectively just making a random walk at high
temperatures.  If the starting temperature is too low, you may already
be in a (too) localized region of the search space, and never reach a
good solution.  Typically a starting temperature that gives you 90\%
or so acceptances (ignoring the rejected attempts, see below) is
good. Better a bit too high than too low.  But don't waste too much
time.

2) \code{end} temperature. By the time that you reach the 
\code{end} temperature the
number of accepted iterations should be only a few per 1000, and the
best score should no longer change. Even zero acceptances is fine. If
there are many more acceptances, lower \code{end}. If there
are zero acceptances for many cycles in a row, raise it a bit. You can
set a lower \code{end} temperature than needed using the \code{earlyout}
test: if
in 5 consecutive cycles of 1000 iterations there are fewer than a
specified number of acceptances per cycle, the program terminates.

3) number of iterations. What really counts is the number of
iterations in the "crunch time", when the number of acceptances is,
say, more than 5\% but fewer than 40\% of the iterations.  If you
print summary statistics in blocks of 1000, you want to see as many
blocks with such acceptance numbers as possible. Obviously within what
is reasonable.

Here are two examples, with my analysis....

(A) \code{logreg.anneal.control(start = 2, end = 1, iter = 50000, update = 1000)
}

The first few lines are (cutting of some of the last columns...)

\tabular{rrrrrrr}{
log-temp \tab current score \tab  best score \tab  acc / \tab rej / \tab sing  \tab current parameters  \cr
   2.000 \tab 1198.785 \tab 1198.785 \tab   0     \tab 0   \tab 0   \tab 0.508 -0.368  -0.144 \cr
   1.980 \tab 1197.962 \tab 1175.311 \tab 719(18) \tab 34  \tab 229 \tab 1.273 -0.275  -0.109 \cr
   1.960 \tab 1197.909 \tab 1168.159 \tab 722(11) \tab 38  \tab 229 \tab 0.416 -0.345  -0.173 \cr
   1.940 \tab 1181.545 \tab 1168.159 \tab 715(19) \tab 35  \tab 231 \tab 0.416 -0.345  -0.173 \cr
   ... \tab  \tab  \tab  \tab  \tab  \tab \cr
   1.020 \tab 1198.258 \tab 1167.578 \tab 663(16) \tab 128 \tab 193 \tab 1.685 -0.216  -0.024 \cr
   1.000 \tab 1198.756 \tab 1167.578 \tab 641(23) \tab 104 \tab 232 \tab 1.685 -0.216  -0.024 \cr
   1.000 \tab 1198.756 \tab 1167.578 \tab 1( 0)   \tab 0   \tab 0   \tab 1.685 -0.216  -0.024 \cr
}

Ignore the last line. This one is just showing a refitting of the best
model.  Otherwise, this suggests
(i) \code{end} is ***way*** too high, as there are still have
more than 600 acceptances in blocks of 1000. It is hard to judge what
\code{end} should be from this run.
(ii) The initial number of acceptances is really high
\code{(719+18)/(719+18+34))=95\%} - but when \code{1.00} is reached it's at about
85\%.  One could  change \code{start} to 1, or keep it at 2 and play it save.


(B) \code{logreg.anneal.control(start = 2, end = -2, iter = 50000, update = 1000)}
-  different dataset/problem

The first few lines are 
\tabular{rrrrrrr}{
log-temp \tab current score \tab  best score \tab  acc / \tab rej / \tab sing  \tab current parameters  \cr
   2.000 \tab   1198.785 \tab 1198.785 \tab    0( 0) \tab  0 \tab  0 \tab 0.50847  -0.36814 \cr  
   1.918 \tab   1189.951 \tab 1172.615 \tab  634(23) \tab 22 \tab   322 \tab 0.38163  -0.28031  \cr
   1.837 \tab   1191.542 \tab 1166.739 \tab  651(24) \tab 32 \tab   293 \tab 1.75646  -0.22451 \cr 
   1.755 \tab   1191.907 \tab 1162.902 \tab  613(30) \tab 20 \tab   337 \tab 1.80210  -0.32276  \cr
}
The last few are
\tabular{rrrrrrr}{
log-temp \tab current score \tab  best score \tab  acc / \tab rej / \tab sing  \tab current parameters  \cr
  -1.837 \tab   1132.731 \tab 1131.866 \tab    0(18) \tab  701 \tab  281 \tab 0.00513  -0.45994 \cr
  -1.918 \tab   1132.731 \tab 1131.866 \tab    0(25) \tab  676 \tab  299 \tab 0.00513  -0.45994 \cr
  -2.000 \tab   1132.731 \tab 1131.866 \tab    0(17) \tab  718 \tab  265 \tab 0.00513  -0.45994 \cr
  -2.000 \tab   1132.731 \tab 1131.866 \tab    0( 0) \tab  0 \tab  1 \tab 0.00513  -0.45994 \cr
}
But there really weren't any acceptances since
\tabular{rrrrrrr}{
log-temp \tab current score \tab  best score \tab  acc / \tab rej / \tab sing  \tab current parameters  \cr
  -0.449 \tab   1133.622 \tab 1131.866 \tab    4(21) \tab  875 \tab  100 \tab 0.00513  -0.45994 \cr
  -0.531 \tab   1133.622 \tab 1131.866 \tab    0(19) \tab  829 \tab  152 \tab 0.00513  -0.45994 \cr
  -0.612 \tab   1133.622 \tab 1131.866 \tab    0(33) \tab  808 \tab  159 \tab 0.00513  -0.45994 \cr
}
Going down from 400 to fewer than 10 acceptances went pretty fast....
\tabular{rrrrrrr}{
log-temp \tab current score \tab  best score \tab  acc / \tab rej / \tab sing  \tab current parameters  \cr
   0.776 \tab   1182.156 \tab 1156.354 \tab  464(31) \tab  258 \tab  247 \tab 1.00543  -0.26602 \cr  
   0.694 \tab   1168.504 \tab 1150.931 \tab  306(17) \tab  355 \tab  322 \tab 1.56695  -0.43351 \cr 
   0.612 \tab   1167.747 \tab 1150.931 \tab  230(38) \tab  383 \tab  349 \tab 1.56695  -0.43351 \cr  
   0.531 \tab   1162.085 \tab 1145.920 \tab  124(12) \tab  571 \tab  293 \tab 1.15376  -0.15223 \cr   
   0.449 \tab   1143.841 \tab 1142.321 \tab   63(20) \tab  590 \tab  327 \tab 2.20150  -0.43795 \cr   
   0.367 \tab   1176.152 \tab 1142.321 \tab  106(21) \tab  649 \tab  224 \tab 2.20150  -0.43795 \cr   
   0.286 \tab   1138.384 \tab 1131.866 \tab   62(18) \tab  731 \tab  189 \tab 0.00513  -0.45994 \cr   
   0.204 \tab   1138.224 \tab 1131.866 \tab   11(27) \tab  823 \tab  139 \tab 0.00513  -0.45994 \cr  
   0.122 \tab   1150.370 \tab 1131.866 \tab   15(12) \tab  722 \tab  251 \tab 0.00513  -0.45994 \cr   
   0.041 \tab   1144.536 \tab 1131.866 \tab   30(19) \tab  789 \tab  162 \tab 0.00513  -0.45994 \cr   
  -0.041 \tab   1137.898 \tab 1131.866 \tab   21(25) \tab  911 \tab 43 \tab 0.00513  -0.45994 \cr  
  -0.122 \tab   1139.403 \tab 1131.866 \tab   12(30) \tab  883 \tab 75 \tab 0.00513  -0.45994  \cr
}

What does this tell me -
(i) \code{start} was probably a bit high - no real harm
done,
(ii) \code{end} was lower than needed. Since there really
weren't any acceptances after 10log(T) was about (\code{-0.5}), an ending
log-temperature of (\code{-1}) would have been fine,
(iii) there were far too few runs. The crunch time didn't take more
than about 10 cycles (10000 iterations). You see that this is the time
the "best model" decreased quite a bit - from 1156 to 1131. I would
want to spend considerably more than 10000 iterations during this
period for a larger problem (how many depends very much on the size of
the problem). So, I'd pick (A)\code{logreg.anneal.control(start = 2,
end = -1, iter = 200000, update = 5000)}.  Since the total range is
reduced from \code{2-(-2)=4} to \code{2-(-1)=3}, over a range of 10log temperatures
of 1 there will be \code{200000/3=67000} rather than \code{50000/4=12500}
iterations.  I would repeat this run a couple of times.

In general I may sometimes run several models, and check the scores of
the best models.  If those are all the same, I'm very happy, if
they're similar but not identical, it's OK, though I may run one or
two longer chains. If they're very different, something is wrong.  For
the permutation test and cross-validation I am usually less picky on
convergence.
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
\code{\link{logreg.mc.control}},
\code{\link{logreg.tree.control}}
}


\examples{
myannealcontrol <- logreg.anneal.control(start = 2, end = -2, iter = 50000, update = 1000)
}


\keyword{logic}
\keyword{methods}
\keyword{nonparametric}
\keyword{tree}

