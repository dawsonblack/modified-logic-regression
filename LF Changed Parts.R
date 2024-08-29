#' Logic Forest
#'
#' Constructs an ensemble of logic regression models using bagging for classification and identification of important predictors and predictor interactions
#' @usage logforest(resp, Xs, nBSXVars, anneal.params, nBS=100, h=0.5, norm=TRUE, numout=5, nleaves)
#' @param resp.type string indicating regression type: 'bin' for classification, 'lin' for linear regression, 'exp_surv' for exponential time-to-event, and 'cph_surv' for Cox-PH time-to-event
#' @param resp numeric vector of binary (for time-to-event, indicates censoring) or continuous response values
#' @param resp.time numeric vector of time to event/censoring values
#' @param Xs matrix or dataframe of zeros and ones for all predictor variables
#' @param nBSXVars integer for the number of predictors used to construct each logic regression model.  The default value is all predictors in the data.
#' @param anneal.params a list containing the parameters for simulated annealing.  See the help file for the function \code{logreg.anneal.control} in the \code{LogicReg} package.  If missing, default annealing parameters are set at \code{start}=1, \code{end}=-2, and \code{iter}=50000.
#' @param nBS number of logic regression trees to be fit in the logic forest model.
#' @param h a number between 0 and 1 for the minimum proportion of trees in the logic forest that must predict a 1 for the prediction to be one.
#' @param norm logical.  If FALSE, predictor and interaction scores in model output are not normalized to range between zero and one.
#' @param numout number of predictors and interactions to be included in model output
#' @param nleaves the maximum number of end nodes generated for each tree
#'
#' @return An object of class \code{"logforest"} including a list of values

#' @export

logforest<-function(resp.type, resp, resp.time, Xs, nBSXVars, anneal.params, nBS=100, h=0.5, norm=TRUE, numout=5, nleaves)
{
  pred<-ncol(Xs)
  if (missing(anneal.params)) {anneal.params<-logreg.anneal.control(start=2, end=-1, iter=100000)}
  if (missing(nBSXVars)) {nBSXVars<-pred}
  if (missing(resp.time)) {resp.time <- matrix(NA, nrow=1, ncol=length(resp))}
  n<-nrow(Xs)
  if (is.null(colnames(Xs))) {x.nms<-paste("X", 1:pred, sep="")}
  else {x.nms<-colnames(Xs)}
  fitlist<-vector("list", nBS)
  IBdata<-vector("list", nBS)
  OOBdata<-vector("list", nBS)
  OOBpred<-matrix(nrow=n, ncol=nBS)
  single.predimport<-vector("list",nBS)
  vimp.import<-vector("list", nBS)
  treepreds.list<-vector("list", nBS)
  IPImatrix<-matrix(NA, nrow=0, ncol=(1+pred*2))
  colnames(IPImatrix)<-c("treenum", x.nms, paste("!",x.nms, sep=""))
  if (resp.type == "bin")  {mtype="Classification"}
  if (resp.type == "lin")  {mtype="Linear Regression"}
  if (resp.type == "exp_surv")  {mtype="Exp. Time-to-Event"}
  if (resp.type == "cph_surv")  {mtype="Cox-PH Time-to-Event"}
  for(b in 1:nBS)
  {
    if(missing(nleaves)) {nleaves<-sample(2:8, 1, replace=FALSE)}
    BSindices<-sample(1:n, n, replace=TRUE)
    OOBindices<-(1:n)[!is.element(1:n, BSindices)]
    BS<-Xs[BSindices, ] #Selects the bootstrap sample corresponding to the chosen indices
    OOB<-Xs[OOBindices, ]  #Selects the out-of-bag sample
    BSY<-resp[BSindices]  #Response variable for bootstrap sample
    OOBY<-resp[OOBindices]  #Response variable for out-of-bag sample
    BST<-resp.time[BSindices]  #Response variable for bootstrap sample
    OOBT<-resp.time[OOBindices]  #Response variable for out-of-bag sample
    XVarIndices<-sort(sample(1:pred, nBSXVars, replace=FALSE))
    rsBSX<-BS[ ,XVarIndices]  #Bootstrap with selected X-vars
    rsOOBX<-OOB[,XVarIndices]  #Out-of-bag with selected X-vars
    if (mtype=="Classification")
    {
      FinalBS<-cbind(rsBSX, BSY)   #Final Bootstrap with Y
      FinalOOB<-cbind(rsOOBX, OOBY)   #Final out-of-bag with Y
      colnames(FinalBS)<-colnames(FinalOOB)<-c(x.nms[XVarIndices], "Y")
      fit <- logreg(resp = BSY, bin = FinalBS[,1:nBSXVars],
                    type = 1, select = 1, ntrees = 1, nleaves = nleaves, anneal.control = anneal.params)
    }
    if (mtype=="Linear Regression")
    {
      FinalBS<-cbind(rsBSX, BSY)   #Final Bootstrap with Y
      FinalOOB<-cbind(rsOOBX, OOBY)   #Final out-of-bag with Y
      colnames(FinalBS)<-colnames(FinalOOB)<-c(x.nms[XVarIndices], "Y")
      fit <- logreg(resp = BSY, bin = FinalBS[,1:nBSXVars],
                    type = 2, select = 1, ntrees = 1, nleaves = nleaves, anneal.control = anneal.params)
    }
    if (mtype=="Cox-PH Time-to-Event")
    {
      FinalBS<-cbind(rsBSX, BSY, BST)   #Final Bootstrap with Y and Time
      FinalOOB<-cbind(rsOOBX, OOBY, OOBT)   #Final out-of-bag with Y and Time
      colnames(FinalBS)<-colnames(FinalOOB)<-c(x.nms[XVarIndices], "Y")
      fit <- logreg(resp = BST, cens = BSY, bin = FinalBS[,1:nBSXVars],
                    type = 4, select = 1, ntrees = 1, nleaves = nleaves, anneal.control = anneal.params)
    }
    if (mtype=="Exp. Time-to-Event")
    {
      FinalBS<-cbind(rsBSX, BSY, BST)   #Final Bootstrap with Y and Time
      FinalOOB<-cbind(rsOOBX, OOBY, OOBT)   #Final out-of-bag with Y and Time
      colnames(FinalBS)<-colnames(FinalOOB)<-c(x.nms[XVarIndices], "Y")
      fit <- logreg(resp = BST, cens = BSY, bin = FinalBS[,1:nBSXVars],
                    type = 5, select = 1, ntrees = 1, nleaves = nleaves, anneal.control = anneal.params)
    }

    OOBpred[as.vector(OOBindices),b]<-predict.logreg(fit, newbin=as.matrix(FinalOOB[,1:nBSXVars]))
    if (sum(fit$model$tree[[1]]$trees[,3])!=0)
    {
      pred.import<-pimp.import(fit=fit, data=FinalBS, testdata=FinalOOB, BSpred=length(XVarIndices),
                               pred=pred, Xs=XVarIndices, mtype=mtype)
      vimp.import[[b]]<-pred.import$pimp.vimp
      single.predimport[[b]]<-pred.import$single.vimp
      treepreds.list[[b]]<-pred.import$Xids
      ipimat<-cbind(rep(b, nrow(pred.import$Ipimat)), pred.import$Ipimat)
      IPImatrix<-rbind(IPImatrix, ipimat)
    }
    else
    {
      single.predimport[[b]]= 0
      vimp.import[[b]]= 0
      treepreds.list[[b]]=0
      IPImatrix<-rbind(IPImatrix, c(b, rep(0, pred*2)))
    }
    fitlist[[b]]<-fit   #list of all models
    IBdata[[b]]<-BSindices
    OOBdata[[b]]<-OOBindices
  }
  if (mtype=="Classification")
  {
    OOB.pred<-matrix(0, nrow=n, ncol=2)
    for (i in 1:n)
    {
      pred.ids<-which(OOBpred[i,]==1|OOBpred[i,]==0)
      pred.vec<-OOBpred[i,c(pred.ids)]
      OOBprop<-sum(pred.vec)/length(pred.vec)
      OOBi.pred<-ifelse(OOBprop>h, 1, 0)
      OOB.pred[i,]<-c(OOBi.pred, OOBprop)
    }
    colnames(OOB.pred)<-c("predicted_class","proportion")
    OOBmiss<-sum(abs(OOB.pred[,1]-resp))/n
  }
  if (mtype=="Linear Regression")
  {
    OOB.pred<-matrix(0, nrow=n, ncol=2)
    for (i in 1:n)
    {
      pred.ids<-which(is.na(OOBpred[i,])==F)
      pred.vec<-OOBpred[i,c(pred.ids)]
      OOBval<-sum(pred.vec)/length(pred.vec)
      OOB.pred[i,]<-c(i, OOBval)
    }
    colnames(OOB.pred)<-c("sample id","mean prediction")
    OOBmiss<-sum((OOB.pred[,2]-resp)^2)/n
  }
  if (mtype=="Cox-PH Time-to-Event")
  {
    OOB.pred<-matrix(0, nrow=n, ncol=2)
    for (i in 1:n)
    {
      pred.ids<-which(is.na(OOBpred[i,])==F)
      pred.vec<-OOBpred[i,c(pred.ids)]
      OOBval<-sum(pred.vec)/length(pred.vec)
      OOB.pred[i,]<-c(i, OOBval)
    }
    colnames(OOB.pred)<-c("sample id","mean prediction")
    OOBmiss<-sum((OOB.pred[,2]-resp)^2)/n
  }
  if (mtype=="Exp. Time-to-Event")
  {
    OOB.pred<-matrix(0, nrow=n, ncol=2)
    for (i in 1:n)
    {
      pred.ids<-which(is.na(OOBpred[i,])==F)
      pred.vec<-OOBpred[i,c(pred.ids)]
      OOBval<-sum(pred.vec)/length(pred.vec)
      OOB.pred[i,]<-c(i, OOBval)
    }
    colnames(OOB.pred)<-c("sample id","mean prediction")
    OOBmiss<-sum((OOB.pred[,2]-resp)^2)/n
  }
  pred.importmat<-matrix(0, nrow=nBS, ncol=pred)
  colnames(pred.importmat)<-x.nms
  for (i in 1:nBS)
  {
    pred.ind<-treepreds.list[[i]]
    m<-length(pred.ind)
    for (j in 1:m)
    {
      col<-pred.ind[j]
      pred.importmat[i,col]<-single.predimport[[i]][j]
    }
  }
  pred.imp<-colSums(pred.importmat)
  names(pred.imp)<-x.nms
  freq.table<-table(names(unlist(single.predimport)))
  all.pimps<-unique(names(unlist(vimp.import)))
  npimps<-length(unique(names(unlist(vimp.import))))
  pimp.freqtable<-table(names(unlist(vimp.import)))
  pimptable<-matrix(0, nrow=nBS, ncol=npimps)
  colnames(pimptable)<-all.pimps
  for (i in 1:nBS)
  {
    npimps.tree<-length(vimp.import[[i]])
    for (j in 1:npimps.tree)
    {
      cpimp<-vimp.import[[i]][j]
      col.id<-which(colnames(pimptable)%in%names(cpimp))
      pimptable[i,col.id]<-cpimp
    }
  }
  pimpsum<-colSums(pimptable)  #importance of each interaction
  t5PIs<-names(sort(pimpsum, decreasing=TRUE)[1:5])
  ans<-list(AllFits=fitlist, Top5.PI=t5PIs, Predictor.importance=pred.imp, Predictor.frequency=freq.table,
            PI.frequency=pimp.freqtable, PI.importance=pimpsum,  ModelPI.import=vimp.import,
            IPImatrix=IPImatrix, OOBmiss=OOBmiss, OOBprediction=OOB.pred, IBdata=IBdata, OOBdata=OOBdata, norm=norm,
            numout=numout, predictors=pred, Xs=Xs, model.type=mtype)
  class(ans)<-"logforest"
  ans
}


### INTERNAL FUNCTION TO EVALUATE IMPORTANCE OF PREDICTOR COMBINATIONS ###
pimp.import<-function(fit, data, testdata, BSpred, pred, Xs, mtype)
{
  if (mtype=="Classification") {mtyp=1}
  if (mtype=="Linear Regression") {mtyp=2}
  if (mtype=="Cox-PH Time-to-Event") {mtyp=3}
  if (mtype=="Exp. Time-to-Event") {mtyp=4}
  n<-nrow(testdata) #number obs in testdata (OOB)
  tree<-fit$model$trees[[1]]
  y<-testdata[,(BSpred+1)]
  x.names<-colnames(data[,1:BSpred])
  orig.pred<-predict.logreg(fit, newbin=testdata[,1:BSpred])
  if (mtyp==1) {orig.miss<-sum(abs(orig.pred-y))/n}
  if (mtyp==2) {orig.miss<-log2(sum((orig.pred-y)^2)/n)}
  if (mtyp==3) {orig.miss<-log2(sum((orig.pred-y)^2)/n)}
  if (mtyp==4) {orig.miss<-log2(sum((orig.pred-y)^2)/n)}
  pimpinfo<-prime.imp(tree=tree, data=data, Xs=Xs, mtype=mtyp)
  numorigpis<-nrow(pimpinfo$tmp.mat)
  Xvars <- Xs

  if (mtyp %in% 2:4)
  {
    cpimpinfo<-prime.imp(tree=find.ctree(tree), data=data, Xs=Xvars, mtype=mtyp)
    cloc<-which(rowSums(abs(cpimpinfo$tmp.mat))>1)
    cnewpimps<-cpimpinfo$vec.primes[cloc]
    loc<-which(rowSums(abs(pimpinfo$tmp.mat))>1)
    newpimps<-pimpinfo$vec.primes[loc]
    is.cmp<-rep(c(0, 1), times=c(length(newpimps), length(cnewpimps)))
    pimpinfo$vec.primes<-append(newpimps, cnewpimps)
    pimpinfo$tmp.mat<-rbind(pimpinfo$tmp.mat[loc,], cpimpinfo$tmp.mat[cloc,])
    pimpinfo$list.pimps<-append(pimpinfo$list.pimps, cpimpinfo$list.pimps)
    pimpinfo$cmp<-is.cmp

    findmtch<-fnd(pinm=strsplit(newvec.primes, split=" & "), is.cmp=is.cmp, mtchnm=mtchnm)
  }

  Ipimat<-matrix(0, nrow=nrow(pimpinfo$tmp.mat), ncol=pred*2)
  for (i in 1:nrow(pimpinfo$tmp.mat))
  {
    lc<-pimpinfo$tmp.mat[i,][which(pimpinfo$tmp.mat[i,]!=0)]*pimpinfo$list.pimps[[i]]
    lc<-ifelse(lc>0, lc, -1*lc+pred)
    Ipimat[i,lc]<-1
  }
  colnames(Ipimat)<-c(x.names, paste("!", x.names, sep=""))
  vec.Xvars<-pimpinfo$vec.pimpvars
  nxvars<-length(vec.Xvars)
  single.vimp<-c()
  single.vimp.nms<-c()
  Xids<-c()
  for (i in 1:nxvars)#checking single var importance for tree
  {
    id<-vec.Xvars[i]
    Xid<-Xs[id]
    permute.ind<-sample(1:n, n, replace=FALSE)
    permute.col<-testdata[permute.ind, id]
    pre.id<-if(id>1) as.matrix(testdata[,1:(id-1)])
    post.id<-if(id<pred) as.matrix(testdata[,(id+1):BSpred])
    permute.testdata<-cbind(pre.id, permute.col, post.id)
    perm.pred<-predict.logreg(fit, newbin=as.matrix(permute.testdata[,1:BSpred]))
    if (mtyp==1) {perm.misclass<-sum(abs(perm.pred-y))/n}   	
    if (mtyp==2) {perm.misclass<-log2(sum((perm.pred-y)^2)/n)}
    if (mtyp==3) {perm.misclass<-log2(sum((perm.pred-y)^2)/n)}
    if (mtyp==4) {perm.misclass<-log2(sum((perm.pred-y)^2)/n)}
    vimp<-perm.misclass-orig.miss
    single.vimp<-append(single.vimp, vimp)
    single.vimp.nms<-append(single.vimp.nms, x.names[id])
    Xids<-append(Xids, Xid)
  }
  names(single.vimp)<-single.vimp.nms
  pimpmat<-pimp.mat(pimps.out=pimpinfo, testdata=testdata)[[2]] #transforms data into predictors = pimps
  pimpnames<-pimp.mat(pimps.out=pimpinfo, testdata=testdata)[[1]] #vector of pimp names
  tmp.mat<-pimpinfo$tmp.mat
  zero.ids<-c()
  for(i in 1:ncol(tmp.mat))
  {
    ids<-if(all(tmp.mat[,i]==0)) {ids<-i}
    zero.ids<-append(zero.ids, ids)
  }
  if (length(zero.ids) > 0) {tmp.mat<-tmp.mat[,-zero.ids]}
  if (is.matrix(tmp.mat)) {npimps<-nrow(tmp.mat)}
  if (is.vector(tmp.mat)) {npimps<-1}
  pimp.vimp<-c()
  for (j in 1:npimps)
  {
    perm.ind<-sample(1:n, n, replace=FALSE)
    perm.col<-pimpmat[perm.ind, j]
    pre.j<-if(j>1) pimpmat[,1:(j-1)]
    post.j<-if(j<npimps) pimpmat[,(j+1):npimps]
    permute.pimpdata<-cbind(pre.j, perm.col, post.j)
    pimp.pred<-c()
    for (k in 1:n)
    {
      if(mtyp==1)
      {
        pred<-ifelse(any(permute.pimpdata[k,]==1), 1, 0)
        pimp.pred<-append(pimp.pred, pred)
      }
      if(mtyp %in% 2:4)
      {
        pred<-ifelse(any(permute.pimpdata[k,]==1), sum(fit$model$coef), fit$model$coef[1])
        pimp.pred<-append(pimp.pred, pred)
      }
    }
    if (mtyp==1) {permpimp.miss<-sum(abs(pimp.pred-y))/n}	
    if (mtyp==2) {permpimp.miss<-log2(sum((pimp.pred-y)^2)/n)}
    if (mtyp==3) {permpimp.miss<-log2(sum((pimp.pred-y)^2)/n)}
    if (mtyp==4) {permpimp.miss<-log2(sum((pimp.pred-y)^2)/n)}
    pvimp<-permpimp.miss-orig.miss #diff bt/ permutation and original
    pimp.vimp<-append(pimp.vimp, pvimp)
  }
  names(pimp.vimp)<-paste(pimpnames)
  out<-list(single.vimp=single.vimp, pimp.vimp=pimp.vimp, Ipimat=Ipimat, vec.Xvars=vec.Xvars, Xids=Xids)
}

### INTERNAL FUNCTION TO EVALUATE IMPORTANCE OF PREDICTOR COMBINATIONS ###
TTab<-function(data, tree, Xs, mtype)
{
  if(!is(tree,"logregtree"))
    stop("tree must be an object of class logregtree")
  mod.var<-tree$trees[,3]
  mod.var<-sort(mod.var[mod.var!=0])
  mod.var<-mod.var[!duplicated(mod.var)]
  model.var<-c()
  for (i in 1:length(mod.var))
  {
    modvar<-Xs[mod.var[i]]
    model.var<-append(model.var, modvar)
  }
  mat.perms<-Perms(length(mod.var))
  nms<-colnames(data)
  if (is.null(colnames(data))) {colnames(mat.perms)<-paste("X", model.var, sep="")}
  if (length(nms)>0) {colnames(mat.perms)<-nms[mod.var]}
  mat.bin<-matrix(0, nrow(mat.perms), max(mod.var))
  mat.bin[,mod.var]<-mat.perms
  pred.out<-eval.logreg(tree, mat.bin)
  mat.truth<-cbind(mat.perms, outcome=pred.out)
  if(mtype==1) {truth<-ifelse(tree$coef>0 | is.na(tree$coef),1,0)}
  if(mtype==2) {truth<-1}
  if(mtype==3) {truth<-1}
  if(mtype==4) {truth<-1}
  ids.truth<-mat.truth[,"outcome"]==truth
  mat.truth<-mat.truth[ids.truth,-ncol(mat.truth),drop=FALSE]
  mat.truth
}


### FUNCTION TO PREDICT OUTCOME FOR NEW OBSERVATIONS ###
predict.logforest<-function(object, newdata, cutoff,...)
{
  if (class(object)!= "logforest")
    stop("object not of class logforest")
  nBS<-length(object$AllFits)
  trees<- object$AllFits
  mtype<-object$model.type
  if(mtype=="Classification")
  {
    if(missing(cutoff)) cutoff<-0.5

    if (missing(newdata))
    {
      LFprediction<-object$OOBprediction[,1]
      proportion_one<-object$OOBprediction[,2]
      ans<-list(model.type=mtype, LFprediction=LFprediction, proportion_one=proportion_one)
    }
    if (!missing(newdata))
    {
      pred<-ncol(newdata)
      if (pred!=object$predictors)
        stop("the predictors in newdata do not match the original predictors")
      size<-nrow(newdata)
      predict.new<-matrix(0, nrow=size, ncol=nBS)
      for (i in 1:nBS)
      {
        newX<-newdata[,1:pred]
        newpredict<-predict.logreg(object=trees[[i]], newbin=as.matrix(newX))
        predict.new[,i]<- newpredict
      }
      predictions<-proportion.positive(predictmatrix=predict.new, cutoff=cutoff)
      predmatrix<-cbind(predict.new, predictions$predmat)
      predframe<-as.data.frame(predmatrix)
      names(predframe)[1:nBS]<-paste("tree", 1:nBS, sep="")
      names(predframe)[nBS+1]<-paste("proportion_one")
      names(predframe)[nBS+2]<-paste("PredictedValue")
      ans<-list(model.type=mtype, LFprediction=predframe$PredictedValue, proportion_one=predframe$proportion_one,
                AllTrees=predframe)
    }
  }

  if(mtype=="Linear Regression")
  {
    if (missing(newdata))
    {
      LFprediction<-object$OOBprediction[,2]
      OOBmse<-object$OOBmiss
      ans<-list(model.type=mtype, LFprediction=LFprediction, OOBmse=OOBmse, ptype="OOBprediction")
    }
    if (!missing(newdata))
    {
      pred<-ncol(newdata)
      if (pred!=object$predictors)
        stop("the predictors in newdata do not match the original predictors")
      size<-nrow(newdata)
      predict.new<-matrix(0, nrow=size, ncol=nBS)
      for (i in 1:nBS)
      {
        newX<-newdata[,1:pred]
        newpredict<-predict.logreg(object=trees[[i]], newbin=as.matrix(newX))
        predict.new[,i]<- newpredict
      }
      predictions<-rowMeans(predict.new)
      predmatrix<-cbind(predict.new, predictions)
      predframe<-as.data.frame(predmatrix)
      names(predframe)[1:nBS]<-paste("tree", 1:nBS, sep="")
      names(predframe)[nBS+1]<-paste("PredictedValue")
      ans<-list(model.type=mtype, LFprediction=predframe$PredictedValue, AllTrees=predframe)
    }
  }
  if(mtype=="Cox-PH Time-to-Event")
  {
    if (missing(newdata))
    {
      LFprediction<-object$OOBprediction[,2]
      OOBmse<-object$OOBmiss
      ans<-list(model.type=mtype, LFprediction=LFprediction, OOBmse=OOBmse, ptype="OOBprediction")
    }
    if (!missing(newdata))
    {
      pred<-ncol(newdata)
      if (pred!=object$predictors)
        stop("the predictors in newdata do not match the original predictors")
      size<-nrow(newdata)
      predict.new<-matrix(0, nrow=size, ncol=nBS)
      for (i in 1:nBS)
      {
        newX<-newdata[,1:pred]
        newpredict<-predict.logreg(object=trees[[i]], newbin=as.matrix(newX))
        predict.new[,i]<- newpredict
      }
      predictions<-rowMeans(predict.new)
      predmatrix<-cbind(predict.new, predictions)
      predframe<-as.data.frame(predmatrix)
      names(predframe)[1:nBS]<-paste("tree", 1:nBS, sep="")
      names(predframe)[nBS+1]<-paste("PredictedValue")
      ans<-list(model.type=mtype, LFprediction=predframe$PredictedValue, AllTrees=predframe)
    }
  }
  if(mtype=="Exp. Time-to-Event")
  {
    if (missing(newdata))
    {
      LFprediction<-object$OOBprediction[,2]
      OOBmse<-object$OOBmiss
      ans<-list(model.type=mtype, LFprediction=LFprediction, OOBmse=OOBmse, ptype="OOBprediction")
    }
    if (!missing(newdata))
    {
      pred<-ncol(newdata)
      if (pred!=object$predictors)
        stop("the predictors in newdata do not match the original predictors")
      size<-nrow(newdata)
      predict.new<-matrix(0, nrow=size, ncol=nBS)
      for (i in 1:nBS)
      {
        newX<-newdata[,1:pred]
        newpredict<-predict.logreg(object=trees[[i]], newbin=as.matrix(newX))
        predict.new[,i]<- newpredict
      }
      predictions<-rowMeans(predict.new)
      predmatrix<-cbind(predict.new, predictions)
      predframe<-as.data.frame(predmatrix)
      names(predframe)[1:nBS]<-paste("tree", 1:nBS, sep="")
      names(predframe)[nBS+1]<-paste("PredictedValue")
      ans<-list(model.type=mtype, LFprediction=predframe$PredictedValue, AllTrees=predframe)
    }
  }
  class(ans)<-"LFprediction"
  return(ans)
}

### Function to find complement of the origianl tree and find complement PIs ###
find.ctree<-function(tree)
{
  conc<-ifelse(tree$trees$conc==3, 3, ifelse(tree$trees$conc==0, 0, 3-tree$trees$conc))
  neg<-ifelse(tree$trees$conc==3, 1-tree$trees$neg, 0)
  ctrees<-cbind(tree$trees$number, conc, tree$trees$knot, neg, tree$trees$pick)
  colnames(ctrees)<-c("number","conc","knot","neg","pick")
  ctrees<-as.data.frame(ctrees)
  ctree<-list(whichtree=tree$whichtree, coef=tree$coef, trees=ctrees)
  class(ctree)<-"logregtree"
  return(ctree)
}
