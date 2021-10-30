###' The goal of this demo is show the steps of how MedFix and MedMix works.
rm(list=ls())
library(glmnet)
library(MASS)
library(rrBLUP)
library(parallel)
library(doParallel)
library(CompQuadForm)
#setwd('E:\\Dropbox\\Dropbox\\Research_active\\f2\\data_f2mouse\\causal_gene_screening\\code2019\\')
setwd('/mnt/myvol/f2mouse/hdmd/code2020')
source('highmed2019.r')
source('fromSKAT.R')




###' New wrapper functions for each method, they source the function script highmed2019.r

adlassoMixWrapper <- function(Z,X,X0,y,kernel='linear',ncores=1,pmax=length(y)-2){  #This function is a wrapper function that runs the adaptive lasso for linear mixed model for a sequence of penalty parameter lambda, and output the result that minimizes BIC
  if(kernel=='shrink_EJ'){
    K = A.mat(Z,shrink=list(method="EJ"))
  }
  if(kernel=='linear'){
    K = A.mat(Z,shrink=FALSE)
  }
  n = length(y)
  eigK = eigen(K+sqrt(n)*diag(n))
  Qeig = eigK$vectors
  thetaeig = eigK$values-sqrt(n)
  Xt = t(Qeig)%*%cbind(cbind(1,X0),X)
  yt = t(Qeig)%*%y
  p = ncol(X)
  lambda.seq = getLambda(Xt,yt,nlambda=100,intercept=F,penalty.factor = c(rep(0,1+ncol(X0)),rep(1,p)))
  #    getDoParWorkers()
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  results.all = foreach(lambda=lambda.seq) %dopar% {
    library(glmnet)
    library(MASS)
    source('highmed2019.r')
    try(adlassoMix(X,y,K,X0=X0,lambda=lambda,init=list(tau=var(y)/2,ve=var(y)/2),pmax = pmax,err.max=1e-5,iter.max=200,iter.lasso.max=1e4,method.var='REML')) ## 'MLE') ## method.var=MLE or REML
  }
  names(results.all) = lambda.seq
  negll.all = sapply(results.all, getElement,name='negll')
  s0.all = sapply(results.all, getElement,name='s0')
  lambda.all = sapply(results.all, getElement,name='lambda')
  bic.all = sapply(results.all, getElement,name='bic')
  mod.final = results.all[[which.min(bic.all)]]
  return(list(model=mod.final,bic=bic.all,negll=negll.all,s0=s0.all,lambda=lambda.all,results.all=results.all))
}

###'
###'
adlasso2typeWrapper <- function(y,X0,X,Z,pz=seq(0.01,0.99,length.out=20),pmax=length(y)-2,ncores=16){ ## wrapper of medfix with tuning parameter selection based on BIC
  ## the fixed model should NOT include the intercept in X0
  n = length(y)
  p0 = ncol(X0)
  p = ncol(X)
  q = ncol(Z)
  getDoParWorkers()
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  out.all = foreach(ppzz=pz) %dopar% {
    library(glmnet)
    library(MASS)
    source('highmed2019.r')
    adlasso.2type(y,X,Z,ppzz,X0,pmax=pmax)
  }
  stopCluster(cl)
  bic = sapply(out.all, getElement,name='bic.min')
  id.opt = which.min(bic)[1]
  out.final = out.all[[id.opt]]
  lambda.opt = out.final[['lambda.seq']][which.min(out.final[['bic']])[1]]
  mod.final = out.final[['model']]
  mod.final[['X0']] = X0
  mod.final[['X']] = X
  mod.final[['Z']] = Z
  mod.final[['y']] = y
  mod.eq= adlasso.2type(y,X,Z,pz=0.5,X0,pmax=pmax)[['model']]
  mod.eq[['X0']] = X0
  mod.eq[['X']] = X
  mod.eq[['Z']] = Z
  mod.eq[['y']] = y
  return(list(model=mod.final,mod.eq=mod.eq,pz.opt=pz[id.opt],lambda.opt=lambda.opt,bic.prop=bic,bic.lambda=out.final[['bic']],pz.seq=pz,lambda.seq=out.final[['lambda.seq']]))
}

e2mFixed <- function(mod,ncores=8){ ## calculate the exposure to mediator effect using the fixed model
  ### mod = mod.fixed
  X = mod[['X']]
  n = nrow(X)
  X0 = mod[['X0']]#[,-1]
  if(is.null(X0)){
    p0=1
  }else{
    p0=ncol(X0)+1
  }
  Z = mod[['Z']]
  b = as.matrix(mod[['b']])
  ixnz = (b!=0)
  ns = sum(ixnz)
  #    print(ns)
  if(ns>0){
    bnz = b[ixnz]
    Xnz = matrix(X[,ixnz],nrow=n)
    #    print('good so far')
    getDoParWorkers()
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    Ball = foreach(ii=1:ns) %dopar% {
      library(glmnet)
      library(MASS)
      source('highmed2019.r')
      ## ii = 1
      vs.lasso = adLasso(X=Z,y=Xnz[,ii],X0=X0)
      all.lasso = bic.glmnet(vs.lasso,cbind(1,cbind(X0,Z)),Xnz[,ii],p0=p0)
      mod.lasso = all.lasso$model
      list(coef=mod.lasso$b,cov=mod.lasso$cov.coef,dv=mod.lasso$dev,s0=mod.lasso$s0)#,bic=all.lasso$bic,df=all.lasso$df,lambda=all.lasso$lambda.seq)
    }
    stopCluster(cl)
  }else{
    Ball = NULL
  }
  return(Ball)
}



###' load the sample data
load('demoFullData.rdata')
ls()
ncores = 16

###' run the estimation step for MedFix
vs.fixed = adlasso2typeWrapper(y,X0,X,Z,ncores=ncores) ## run the variable selection step for the outcome model of MedFix
mod.fixed = vs.fixed[['model']]  # extract the model that minimizes BIC
mod.eq = vs.fixed[['mod.eq']] # extract the model that assign equal penalty on the two data types.
e2m.fixed= e2mFixed(mod.fixed,ncores=ncores) # for mod.fixed, run the mediator models for the mediators with non-zero estimated coefficient in the outcome model
e2m.eq= e2mFixed(mod.eq,ncores=ncores) # for mod.eq, run the mediator models for the mediators with non-zero estimated coefficient in the outcome model
###' run the estimation step for MedMix using linear kernel, and extract the model that minimizes BIC
vs.mix.linear = adlassoMixWrapper(Z,X,X0,y,kernel='linear',ncores=ncores)
mod.mix.linear = vs.mix.linear[['model']]
###' run the estimation step for MedMix using shrink_EJ kernel, and extract the model that minimizes BIC
vs.mix.shrink = adlassoMixWrapper(Z,X,X0,y,kernel='shrink_EJ',ncores=ncores)
mod.mix.shrink = vs.mix.shrink[['model']]

###' perform the hypotheses testing step that control fasle discovery proportion (FDP) below gamma
###' P(FDP>gamma)<alpha where gamma = 0.1, and alpha could be 0.05
p.adj.method = 'holmr' # I wrote an implementation of Holm test that controls FDP
pval.cut=0.05
tests.fixed=testMedFix(mod.fixed,e2m.fixed,p.adj.method=p.adj.method)
(res.fixed.pcut = medH.L2fixed(mod.fixed,e2m.fixed,tests.fixed,pval.cut=pval.cut))# report the results including the PVM calculation for MedFix_BIC

tests.eq=testMedFix(mod.eq,e2m.eq,p.adj.method=p.adj.method)
(res.eq.pcut = medH.L2fixed(mod.eq,e2m.eq,tests.eq,pval.cut=pval.cut))# report the results including the PVM calculation for MedFix_{0.5}

tests.mix.linear = testMedMix(mod.mix.linear,p.adj.method=p.adj.method)
(res.mix.linear.pcut = medH.L2(mod.mix.linear,tests.mix.linear,pval.cut=pval.cut))# report the results including the PVM calculation for MedMix_linear

tests.mix.shrink = testMedMix(mod.mix.shrink,p.adj.method=p.adj.method)
(res.mix.shrink.pcut = medH.L2(mod.mix.shrink,tests.mix.shrink,pval.cut=pval.cut))# report the results including the PVM calculation for MedMix_shrink


###' a function extracting the mediators making the cut

reportMediator <- function(mod,tst,pval.cut){
  coef.outcome = mod[['cov.coef']][['b.nz']]
  med = tst[['med.individual']]
  ix.output = (med$padj<=pval.cut)
  output = list(pval=med[ix.output,],coef.outcome=coef.outcome[ix.output])
  return(output)
}

###' Report the mediators selected by each method
mediators.fixed = reportMediator(mod.fixed,tests.fixed,pval.cut)
id.fixed = mediators.fixed[['pval']]$id

mediators.eq = reportMediator(mod.eq,tests.eq,pval.cut)
id.eq = mediators.eq[['pval']]$id

mediators.mix.linear = reportMediator(mod.mix.linear,tests.mix.linear,pval.cut)
id.mix.linear = mediators.mix.linear[['pval']]$id

mediators.mix.shrink = reportMediator(mod.mix.shrink,tests.mix.shrink,pval.cut)
id.mix.shrink = mediators.mix.shrink[['pval']]$id


###' performance comparison using a simulated dataset in demoFullData.rdata, gamma is the true coefficient of the mediators.
###' A function that compare the variable selection performance
vsCompare <- function(estimator,target){ ### inputs are ids of non-zero elements
  fp = length(setdiff(estimator,target))
  fn = length(setdiff(target,estimator))
  output = c(fp,fn,length(estimator),length(target))
  names(output) = c('fp','fn','discovery','truth')
  return(output)
}
###'  compare the mediator selection accuracy of these methods
id.gamma = which(gamma!=0)

vsCompare(id.fixed,id.gamma)

vsCompare(id.eq,id.gamma)

vsCompare(id.mix.linear,id.gamma)

vsCompare(id.mix.shrink,id.gamma)
