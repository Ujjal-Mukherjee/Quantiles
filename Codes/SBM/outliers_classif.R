## PQ/ WPQ-based outliers
## combine with classification scheme

## Algorithm:
# 1. Take a random test/train split of data, predict test data
# 2. Detect outliers in each cluster in train data, predict test data again
# 3. repeat for a number of random splits

setwd("C:/Study/My projects/Quantiles/Codes/SBM")
rm(list=ls());

library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(fda.usc)

sourceCpp("ProjQuantNew1.cpp", verbose=TRUE, rebuild=FALSE);

## function to compute outlier score
outlier.score = function(X, k=NULL, alpha){
  n = nrow(X)
  if(is.null(k)){
    #k = floor(.1*n)
    k = floor(sqrt(n))
  }
  depth.vec = rep(0,n)
  knn.vec = rep(0,n)
  
  # calculate distance matrix for full data
  dist.mat = as.matrix(dist(X))
  
  # get depth and knn average dist for all data
  for(i in 1:n){
    depth.vec[i] = ProjQuantileDepthMod(X[-i,], X[i,], .5)    
    ik = order(dist.mat[i,])[2:(k+1)]
    knn.vec[i] = mean(dist.mat[i,ik])
  }
  
  lknn = log(knn.vec)
  lhtped = -log(depth.vec)
  
  return(alpha*lknn + (1-alpha)*lhtped)
}

# DNA alteration data
set.seed(11252014)
require(FastHCS)
data(DnaAlteration)
nsplit = 1e4
n = nrow(DnaAlteration)
ntest = floor(n/10)

for(i in 1:nsplit){
  test = sample(1:n, ntest, replace=F)
  itest = as.matrix(DnaAlteration[test,])
  itrain = as.matrix(DnaAlteration[-test,])
  
  depth0 = rep(0,ntest); depth1 = depth0
  for(j in 1:ntest){
    depth0[j] = log(ProjQuantileDepthMod(itrain[which(itrain[,1]==0),-1], itest[j,-1], .5))
    depth1[j] = log(ProjQuantileDepthMod(itrain[which(itrain[,1]==1),-1], itest[j,-1], .5))
  }
  
  class.pred = (depth0 > depth1)
  pred.miss = sum(itest[,1] == class.pred)
}

