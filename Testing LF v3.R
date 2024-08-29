library(LogicForest)

set.seed(10051988)
N_c <- 100
N_r <- 500
init <- data.frame(matrix(0, nrow = N_r, ncol = N_c))
colnames(init, prefix = 'X')

for(n in 1:N_c){
  p <- runif(1, min = 0.2, max = 0.6)
  init[,n] <- rbinom(N_r, 1, p)
}

y_p <- -1 + .5*init$X1 + .5*init$X2 + 2*init$X3*init$X4 + 2*init$X5*init$X6
p <- 1/(1 + exp(-y_p))
init$Y <- rbinom(N_r, 1, p)

summary(init$Y)

init$shape <- 1 - .05*init$X1 - .05*init$X2 - .2*init$X3*init$X4 - .2*init$X5*init$X6
init$scale <- 1.5 - .05*init$X1 - .05*init$X2 - .2*init$X3*init$X4 - .2*init$X5*init$X6
init$TIME_Y <- rgamma(N_r, shape = init$shape, scale = init$scale)

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

ptm <- proc.time()
LF.fit <- logforest(resp.type = "bin", resp = init$Y, Xs = init[,1:N_c], nBS=100, nleaves = 8, numout = 30)
print(LF.fit)
proc.time() - ptm
LF_INT_FREQ <- data.frame(LF.fit$PI.frequency)
LF_MAIN_FREQ <- data.frame(LF.fit$Predictor.frequency)

ptm <- proc.time()
LF.fit <- logforest(resp.type = "lin", resp = init$TIME_Y, Xs = init[,1:N_c], nBS=500, nleaves = 8, numout = 30)
print(LF.fit)
proc.time() - ptm
LF_INT_FREQ <- data.frame(LF.fit$PI.frequency)
LF_MAIN_FREQ <- data.frame(LF.fit$Predictor.frequency)

ptm <- proc.time()
LF.fit <- logforest(resp.type = "exp_surv", resp = init$Y, resp.time = init$TIME_Y, Xs = init[,1:N_c], nBS=100, nleaves = 8, numout = 30)
print(LF.fit)
proc.time() - ptm
LF_INT_FREQ <- data.frame(LF.fit$PI.frequency)
LF_MAIN_FREQ <- data.frame(LF.fit$Predictor.frequency)

fit <- logreg(resp = init$TIME_Y, cens = init$Y, bin = init[,1:N_c],
              type = 5, select = 1, ntrees = 1, nleaves = nleaves, anneal.control = anneal.params)
