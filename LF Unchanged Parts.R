library(LogicReg)
#library(data.table)

### INTERNAL FUNCTION TO EVALUATE IMPORTANCE OF PREDICTOR COMBINATIONS ###
prime.imp<-function (tree, data, Xs, mtype)
{
  mat<-TTab(data, tree, Xs=Xs, mtype=mtype)
  mat[mat == 0] <- -1
  n.var <- ncol(mat)
  if (is.null(colnames(mat)))
    colnames(mat) <- paste("X", 1:n.var, sep = "")
  ia <- p.combos(n.var)
  ia.rowS <- rowSums(ia)
  vec.primes <- character(0)
  list.pimps<-list()
  list.cover <- list()
  mat.in <- NULL
  name.paste <- function(x)
  {
    x <- x[x != ""]
    paste(x, collapse = " & ")
  }
  tmp.list<-vector("list",n.var)
  for (i in 1:n.var)
  {
    pairt <- p.combos(i, conj = -1)
    n.p <- nrow(pairt)  #number of combos for i vars
    ia2 <- matrix(ia[ia.rowS == i, ], ncol = n.var)
    tmp <- matrix(0, nrow(ia2) * n.p, n.var)
    for (j in 1:nrow(ia2)) tmp[((j - 1) * n.p + 1):(j * n.p), ia2[j, ] == 1] <- pairt
    if (length(vec.primes) > 0) {
      tmp9 <- tmp %*% t(mat.in)
      tmp10 <- diag(mat.in %*% t(mat.in))
      tmp11 <- t(tmp10) %x% rep(1, nrow(tmp)) == tmp9
      tmp.in <- which(rowSums(tmp11) == 0)
      tmp <- tmp[tmp.in, ]
    }
    tmp2 <- tmp %*% t(mat) == i
    ids <- which(rowSums(tmp2) == 2^(n.var - i))
    tmp <- tmp[ids, ]
    tmp.list[[i]]<-tmp
    if (length(ids) > 0)
    {
      mat.in <- rbind(mat.in, tmp)
      for (k in ids) list.cover[[length(list.cover) + 1]] <- which(tmp2[k, ])
      mat.names <- matrix(rep(colnames(mat), e = length(ids)), ncol = n.var)
      mat.names[tmp == 0] <- ""
      mat.pimps<-vector("list", nrow(mat.names))
      for (m in 1:nrow(mat.names))
      {
        mat.pimps[[m]]<-which(colnames(data)%in%mat.names[m,])
      }
      list.pimps<-append(list.pimps, mat.pimps)
      mat.names[tmp == -1] <- paste("!", mat.names[tmp == -1], sep = "")
      tmp.prime <- apply(mat.names, 1, name.paste)
      vec.primes <- c(vec.primes, tmp.prime)
    }
    cover <- unique(unlist(list.cover))
    if (length(cover) == nrow(mat))
      break
  }
  tmp.mat<-NULL
  for (i in 1:length(tmp.list)){
    tmp.mat<-rbind(tmp.mat, tmp.list[[i]])
  }
  vec.pimpvars<-sort(unique(unlist(list.pimps)))
  n.prime <- length(vec.primes)
  listPI <- list(vec.primes=vec.primes, tmp.mat=tmp.mat, vec.pimpvars=vec.pimpvars, list.pimps=list.pimps)
  class(listPI) <- "primeImp"
  listPI
}


### INTERNAL FUNCTION TO EVALUATE IMPORTANCE OF PREDICTOR COMBINATIONS ###
Perms<-function(n)
{
  mat<-matrix(0, nrow=2^n, ncol=n)
  for(i in 1:n)
  {
    mat[,i]<-rep(0:1, times=2^(i-1), each=2^(n-i))
  }
  mat
}


### INTERNAL FUNCTION TO EVALUATE IMPORTANCE OF PREDICTOR COMBINATIONS ###
p.combos<-function(n.pair,conj=0)
{
  mat<-matrix(0,2^n.pair,n.pair)
  for(i in 1:n.pair)
  {
    mat[, i]<-rep(rep(c(1,conj),e=2^(n.pair-i)),2^(i-1))
  }
  mat
}


### INTERNAL FUNCTION TO EVALUATE IMPORTANCE OF PREDICTOR COMBINATIONS ###
pimp.mat<-function(pimps.out, testdata)
{
  tmp.mat<-pimps.out$tmp.mat
  zero.ids<-c()
  for(i in 1:ncol(tmp.mat))
  {
    ids<-if(all(tmp.mat[,i]==0)) {ids<-i}
    zero.ids<-append(zero.ids, ids)
  }
  if (length(zero.ids) > 0) {tmp.mat<-tmp.mat[,-zero.ids]}
  pimp.ids<-pimps.out$vec.pimpvars
  subdata<-as.matrix(testdata[,pimp.ids])
  if (is.null(dim(tmp.mat))) {tmp.mat<-matrix(1,1,1)}
  if (nrow(tmp.mat)!=length(pimps.out$vec.primes))
  {tmp.mat<-t(tmp.mat)}
  if (is.matrix(tmp.mat)) {npimps<-nrow(tmp.mat)}
  if (is.vector(tmp.mat)) {npimps<-1}
  n<-nrow(subdata)
  pimp.datamat<-matrix(0, nrow=n, ncol=npimps)
  colnames(pimp.datamat)<-pimps.out$vec.primes
  for (i in 1:npimps)
  {
    if (is.matrix(tmp.mat)) {match.matrix<-matrix(0, nrow=n, ncol=ncol(tmp.mat))}
    if (is.vector(tmp.mat)) {match.matrix<-matrix(0, nrow=n, ncol=length(tmp.mat))}
    for (j in 1:n)
    {
      if (is.matrix(tmp.mat))
      {
        for (k in 1:ncol(tmp.mat))
        {
          if (tmp.mat[i,k]==1 & subdata[j,k]==1) {match.matrix[j,k]<-1}
          if (tmp.mat[i,k]==-1 & subdata[j,k]==0) {match.matrix[j,k]<-1}
          if (tmp.mat[i,k]==0 & subdata[j,k]==1|tmp.mat[i,k]==0 & subdata[j,k]==0) {match.matrix[j,k]<-1}
        }
      }
      if(is.vector(tmp.mat))
      {
        for (k in 1:length(tmp.mat))
        {
          if (tmp.mat[k]==1 & subdata[j,k]==1) {match.matrix[j,k]<-1}
          if (tmp.mat[k]==-1 & subdata[j,k]==0) {match.matrix[j,k]<-1}
          if (tmp.mat[k]==0 & subdata[j,k]==1|tmp.mat[k]==0 & subdata[j,k]==0) {match.matrix[j,k]<-1}
        }
      }
      pimp.datamat[j,i]<-ifelse(all(match.matrix[j,]==1), 1, 0)
    }
  }
  pimp.names<-pimps.out$vec.primes
  pimp.info<-list(pimp.names=pimp.names, pimp.datamat=pimp.datamat)
  pimp.info
}




### INTERNAL FUNCTION USED IN DETERRMINING PREDICTION FOR CLASSIFICATION MODELS ###
proportion.positive<-function(predictmatrix, cutoff)
{
  q<-nrow(predictmatrix)
  ntrees<-ncol(predictmatrix)
  status<-c()
  predict.pos<-c()
  for (a in 1:q)
  {
    number.diseasepositive<-sum(predictmatrix[a,])
    proportion.predictpositive<-number.diseasepositive/ntrees
    if (proportion.predictpositive >= cutoff)
      disease.status <- 1
    else if (proportion.predictpositive < cutoff)
      disease.status <- 0
    status<-append(status, disease.status)
    predict.pos<-append(predict.pos, proportion.predictpositive)
  }
  predmat<-cbind(predict.pos, status)
  ans<-list(predmat=predmat)
}



### PRINTS OUT LOGIC FOREST MODEL ###
print.logforest<-function(x, ...)
{
  if (class(x)!= "logforest")
    stop("n not of class logforest")
  num<-x$numout
  OOBerr<-x$OOBmiss
  pscore<-sort(x$Predictor.importance, decreasing=TRUE)
  pred.freq<-x$Predictor.frequency[names(pscore)[1:num]]
  piscore<-sort(x$PI.importance, decreasing=TRUE)
  vimp.freq<-x$PI.frequency[names(piscore)[1:num]]
  p.1st<-paste("Top ", num, " Predictors", sep="")
  pi.1st<-paste("Top ", num, " Interactions", sep="")
  if (x$norm==TRUE)
  {
    pred.score<-round(pscore[1:num]/pscore[1], digits=4)
    vimp.score<-round(piscore[1:num]/piscore[1], digits=4)
    p.cnames<-c(p.1st, "Normalized Predictor Importance","Frequency")
    pi.cnames<-c(pi.1st, "Normalized Interaction Importance","Frequency")
  }
  else {
    pred.score<-round(pscore[1:num], digits=2)
    vimp.score<-round(piscore[1:num], digits=2)
    p.cnames<-c(p.1st, "Predictor Importance","Frequency")
    pi.cnames<-c(pi.1st, "Interaction Importance","Frequency")
  }
  prds<-cbind(names(pred.score), pred.score, pred.freq)
  pis<-cbind(names(vimp.score), vimp.score, vimp.freq)
  colnames(prds)<-p.cnames
  colnames(pis)<-pi.cnames
  rownames(pis)<-rownames(prds)<-c(1:num)
  cat("Number of logic regression trees =", length(x$AllFits), sep=" ")
  cat("\n")
  if (x$model.type=="Classification") cat("Out of Bag Misclassification =", OOBerr, sep=" ")
  if (x$model.type=="Linear Regression") cat("Out of Bag MSE =", OOBerr, sep=" ")
  cat("\n")
  cat("\n")
  cat(num, " most important predictors \n", sep="")
  cat("\n")
  print.default(prds, quote=FALSE, print.gap=3)
  cat("\n")
  cat(num, " most important interactions \n", sep="")
  cat("\n")
  print.default(pis, quote=FALSE, print.gap=3)
}

### PRINT FUNCTION FOR PREDICTION FROM LOGIC FOREST ###
print.LFprediction<-function(x, ...)
{
  if(class(x)!="LFprediction")
    stop("x not of class LFprediction")
  prdct<-x$LFprediction
  if(model.type=="Classification")
  {
    prop<-x$proportion_one

    if(length(x)==3)
    {
      cat("OOB Predicted values\n")
      cat("\n")
      print.default(prdct, quote=FALSE)
      cat("\n")
      cat("Proportion of OOB trees that predict 1")
      cat("\n")
      print.default(prop, quote=FALSE)
    }
    if(length(x)==4)
    {
      cat("Predicted values\n")
      cat("\n")
      print.default(prdct, quote=FALSE)
      cat("\n")
      cat("Proportion of trees that predict 1")
      cat("\n")
      print.default(prop, quote=FALSE)
    }
  }

  if(model.type=="Linear Regression")
  {
    if(length(x)==4)
    {
      cat("OOB Predicted values\n")
      cat("\n")
      print.default(prdct, quote=FALSE)
      cat("\n")
      cat("OOB mean squared error")
      cat("\n")
      print.default(x$OOBmse, quote=FALSE)
    }
    if(length(x)==3)
    {
      cat("Predicted values\n")
      cat("\n")
      print.default(prdct, quote=FALSE)
    }
  }
}



