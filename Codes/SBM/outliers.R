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
    u = runif(p,-1,1); u = u/sqrt(sum(u^2))
    uecdf = ecdf(X%*%u)
    Fuxu.vec[iu] = uecdf(xx%*%u)
  }
  1/(1+max(abs(Fuxu.vec-.5)))
}

## function to compute outlier score
outlier.score = function(X, type, k=NULL, alpha){
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
    if(type==1){
      depth.vec[i] = ProjQuantileDepthMod(X[-i,], X[i,], .5)
    }      
    else{
      depth.vec[i] = EPQD(X[-i,], X[i,])
    }
    
    ik = order(dist.mat[i,])[2:(k+1)]
    knn.vec[i] = mean(dist.mat[i,ik])
  }
  
  lknn = log(knn.vec)
  lhtped = -log(depth.vec)
  
  return(alpha*lknn + (1-alpha)*lhtped)
}

## simulated data
## Setup 1... bivariate normal.. 5% contamination far away
set.seed(11182014)
X1 = matrix(rnorm(950), ncol=2)
X2 = matrix(rnorm(50)+10, ncol=2)
X = rbind(X1,X2)

score.vec = outlier.score(X, type=1, k=5, alpha=.2) # type 1 = ray depth, type 2 = PQD
cols = c(rep("darkgreen",475), rep("darkred", 25))
label.vec = ifelse(score.vec>q9, 19, 1)

par(mfrow=c(1,2))
plot(score.vec, col=cols, pch=19, cex=.5)
abline(h=quantile(score.vec,.95), lty=2, lwd=1)
#abline(h=quantile(score.vec,.05), lty=2, lwd=2)

q9 = quantile(score.vec,.95)
plot(X, pch=label.vec, col=cols)
par(mfrow=c(1,1))

# writeup plots
score1 = outlier.score(X, type=1, alpha=.1)
score2 = outlier.score(X, type=1, alpha=.5)
score3 = outlier.score(X, type=1, alpha=.9)
cols = c(rep("green",475), rep("red", 25))

par(mfrow=c(1,3))
plot(score1, col=cols, pch=19, cex=.5, ylab="outlier score")
plot(score2, col=cols, pch=19, cex=.5, ylab="outlier score")
plot(score3, col=cols, pch=19, cex=.5, ylab="outlier score")
par(mfrow=c(1,1))

## Setup 2... bivariate normal.. 5% contamination far away, scattered
X1 = matrix(rnorm(950), ncol=2)
X2 = matrix(rnorm(50)+sample(c(5:8, -8:-5), 50, replace=T), ncol=2)
label.vec = c(rep("1",475), rep("2", 25))
X = rbind(X1,X2)

score.vec = outlier.score(X, type=1, alpha=0) # type 1 = ray depth, type 2 = PQD
cols = c(rep("darkgreen",475), rep("darkred", 25))
label.vec = ifelse(score.vec>q9, 19, 1)

par(mfrow=c(1,2))
plot(score.vec, col=cols, pch=19, cex=.5)
abline(h=quantile(score.vec,.95), lty=2, lwd=1)
#abline(h=quantile(score.vec,.05), lty=2, lwd=2)

q9 = quantile(score.vec,.95)
plot(X, pch=label.vec, col=cols)
par(mfrow=c(1,1))

# writeup plots
score1 = outlier.score(X, type=1, alpha=.1)
score2 = outlier.score(X, type=1, alpha=.5)
score3 = outlier.score(X, type=1, alpha=.9)
cols = c(rep("green",475), rep("red", 25))

par(mfrow=c(1,3))
plot(score1, col=cols, pch=19, cex=.5, ylab="outlier score")
plot(score2, col=cols, pch=19, cex=.5, ylab="outlier score")
plot(score3, col=cols, pch=19, cex=.5, ylab="outlier score")
par(mfrow=c(1,1))

## colon data
require(cepp)
data(Colon)
n = length(Colon$Y)

score.vec = outlier.score(Colon$X, type=1, alpha=.5)

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

score.vec = outlier.score(DnaAlt.x, type=1, alpha=.5)

cols = ifelse(DnaAlt.y==1, "red","green")
plot(score.vec, col=cols, pch=19, cex=.5)
abline(h=quantile(score.vec,.9), lty=2, lwd=2)

t.test(depth.vec[which(DnaAlt.y==1)], depth.vec[which(DnaAlt.y==0)])
ks.test(depth.vec[which(DnaAlt.y==1)], depth.vec[which(DnaAlt.y==0)])

# stackloss data
par(mfrow=c(2,2))
plot(lm(stack.loss~., data=stackloss))
par(mfrow=c(1,1))

score.vec = outlier.score(as.matrix(stackloss[,-4]),
                          type=1, alpha=.8)

plot(score.vec, pch=19, cex=.5)
abline(h=quantile(score.vec,.9), lty=2, lwd=2)

# writeup plots
X = as.matrix(stackloss[,-4])
score1 = outlier.score(X, type=1, alpha=0)
score2 = outlier.score(X, type=1, alpha=.5)
score3 = outlier.score(X, type=1, alpha=1)

par(mfrow=c(1,3))
plot(score1, pch=19, ylab="outlier score")
plot(score2, pch=19, ylab="outlier score")
plot(score3, pch=19, ylab="outlier score")
par(mfrow=c(1,1))

# hawkins bradu kass data
require(robustbase)

score.vec = outlier.score(as.matrix(hbk[,-4]),
                          type=1, alpha=.2)
plot(score.vec, pch=19, cex=.5)
#abline(h=quantile(score.vec,.9), lty=2, lwd=2)

# writeup plots
X = as.matrix(hbk[,-4])
score1 = outlier.score(X, type=1, alpha=.05)
score2 = outlier.score(X, type=1, alpha=.5)
score3 = outlier.score(X, type=1, alpha=.95)

par(mfrow=c(1,3))
plot(score1, pch=19, cex=.8, ylab="outlier score"); abline(v=14, lty=2)
plot(score2, pch=19, cex=.8, ylab="outlier score"); abline(v=14, lty=2)
plot(score3, pch=19, cex=.8, ylab="outlier score"); abline(v=14, lty=2)
par(mfrow=c(1,1))

# brain weight data
Animals1 = within(Animals, {lbodywt=log(body)
                            lbrainwt=log(brain)})

score.vec = outlier.score(as.matrix(Animals1[,3:4]), type=1)

par(mfrow=c(1,2))
with(Animals1, plot(lbodywt, lbrainwt, col="white"))
with(Animals1, text(lbodywt, lbrainwt, row.names(Animals), cex=.7))

plot(score.vec, col="white")
text(score.vec, row.names(Animals), cex=.7)
abline(h=quantile(score.vec,.9), lty=2, lwd=2)
par(mfrow=c(1,1))
