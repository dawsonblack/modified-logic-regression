\name{logreg.testdat}
\alias{logreg.testdat}
\title{Test data for Logic Regression}


\description{
logreg.testdat has 500 cases, and 21 columns. Column 1 is the response
Y, column k+1, k=1,...,20 is (binary) predictor Xk.  Each predictor Xk
is simulated as an independent Bernoulli(pk) random variables, with
success probabilities pk between 0.1 and 0.9. The response variable is
simulated from the model
 
Y = 3 + 1 L1 - 2 L2 + N(0,1),
 
where L1=(X1 or X2) and L2=(X3 or X4). So the task is to use the
linear model in the logic regression framework to find L1 and L2.
}

\examples{
data(logreg.testdat)
}

\keyword{datasets}


\seealso{
\code{\link{logreg}}
}
