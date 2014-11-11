## PQ/ WPQ-based outliers
setwd("C:/Study/My projects/Quantiles/Codes/SBM")
rm(list=ls());

library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(fda.usc)

sourceCpp("ProjQuantNew.cpp", verbose=TRUE, rebuild=FALSE);

## function for empirical PQD
EPQD = function(X, xx){
  p = ncol(X)
  Fuxu.vec = rep(0,100)
  
  for(iu in 1:100){
    u = rnorm(p); u = u/sqrt(sum(u^2))
    uecdf = ecdf(X%*%u)
    Fuxu.vec[iu] = uecdf(xx%*%u)
  }
  1/(1+max(abs(Fuxu.vec-.5)))
}

## function to compute outlier score
outlier.score = function(X, k=NULL){
  n = nrow(X)
  if(is.null(k)){
    k = floor(.1*n)
  }
  depth.vec = rep(0,n)
  knn.vec = rep(0,n)
  
  # calculate distance matrix for full data
  dist.mat = as.matrix(dist(X))
  
  # get depth and knn average dist for all data
  for(i in 1:n){
    depth.vec[i] = ProjQuantileDepthMod(X[-i,], X[i,], .5)
    #depth.vec[i] = EPQD(X[-i,], X[i,])
    ik = order(dist.mat[i,])[1:k]
    knn.vec[i] = mean(dist.mat[i,ik])
  }
  
  return(log(knn.vec^2/depth.vec))
}

## colon data
require(cepp)
data(Colon)
n = length(Colon$Y)

score.vec = outlier.score(Colon$X)

cols = ifelse(Colon$Y==2, "red","green")
plot(score.vec, col=cols, pch=19, cex=.5)
abline(h=quantile(score.vec,.9), lty=2, lwd=2)
abline(h=quantile(score.vec,.1), lty=2, lwd=2)


t.test(depth.vec[which(Colon$Y==1)], depth.vec[which(Colon$Y==2)])
ks.test(depth.vec[which(Colon$Y==1)], depth.vec[which(Colon$Y==2)])

# depth.vec1 = rep(0,n)
# for(i in 1:n){
#   depth.vec1[i] = -log(ProjQuantileDepth(Colon$X, Colon$X[i,], 10,0.5,0.2,0.1))
# }
# depth.vec1
# 
# cols = ifelse(Colon$Y==2, "red","green")
# plot(depth.vec1, col=cols)
# t.test(depth.vec1[which(Colon$Y==1)], depth.vec1[which(Colon$Y==2)])
# ks.test(depth.vec1[which(Colon$Y==1)], depth.vec1[which(Colon$Y==2)])

# DNA alteration data
require(FastHCS)
data(DnaAlteration)

n = dim(DnaAlteration)[1]
depth.vec = rep(0,n)

DnaAlt.x = as.matrix(DnaAlteration[,-1])
DnaAlt.y = as.matrix(DnaAlteration[,1])

score.vec = outlier.score(DnaAlt.x)

cols = ifelse(DnaAlt.y==1, "red","green")
plot(score.vec, col=cols, pch=19, cex=.5)
abline(h=quantile(score.vec,.9), lty=2, lwd=2)
abline(h=quantile(score.vec,.1), lty=2, lwd=2)

t.test(depth.vec[which(DnaAlt.y==1)], depth.vec[which(DnaAlt.y==0)])
ks.test(depth.vec[which(DnaAlt.y==1)], depth.vec[which(DnaAlt.y==0)])
