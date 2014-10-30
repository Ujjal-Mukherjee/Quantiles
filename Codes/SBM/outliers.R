## PQ/ WPQ-based outliers
setwd("C:/Study/My projects/Quantiles/Codes/SBM")
rm(list=ls());

library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)

sourceCpp("ProjQuantNew.cpp", verbose=TRUE, rebuild=FALSE);

## colon data
require(cepp)
data(Colon)
n = length(Colon$Y)
depth.vec = rep(0,n)

for(i in 1:n){
  depth.vec[i] = ProjQuantileDepthMod(Colon$X, Colon$X[i,], .5)
}
depth.vec

cols = ifelse(Colon$Y==2, "red","green")
plot(log(depth.vec), col=cols)
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

for(i in 1:n){
  depth.vec[i] = ProjQuantileDepthMod(DnaAlt.x[-i,],DnaAlt.x[i,], .999)
}
depth.vec
cols = ifelse(DnaAlt.y==1, "red","green")
plot(depth.vec, col=cols)
t.test(depth.vec[which(DnaAlt.y==1)], depth.vec[which(DnaAlt.y==0)])
ks.test(depth.vec[which(DnaAlt.y==1)], depth.vec[which(DnaAlt.y==0)])

