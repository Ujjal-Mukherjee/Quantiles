rm(list=ls())
set.seed(10092014)
n<-10000;

library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)

## cpp source file
sourceCpp("ProjQuantile.cpp", verbose=TRUE, rebuild=FALSE);

## function to generate from multivariate normal
my.mvrnorm = function(n, mu, Sigma){
  p = length(mu)
  # compute square root of covariance matrix
  eo=eigen(Sigma, symmetric=TRUE)
  sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
  
  # generate random normals from runif by box-muller transform
  rnorm.vec = sqrt(-2*log(runif(n*p)))*cos(2*pi*runif(n*p))
  
  # generate sample matrix
  sample.matrix = matrix(rep(mu, n), nrow=n, byrow=T) +
    matrix(rnorm.vec, nrow=n, ncol=p)%*%sigma.sqrt
  return(sample.matrix)
}

n = 1e4
## Simulating from Bivariate normal, mean (0,0), Sigma = ((1 .5), (.5 1))
sig = matrix(c(1,.5,.5,1), nrow=2)
X1 = my.mvrnorm(n, mu=c(0,0), Sigma=sig)
X1.pq = ProjQuant(X1, c(0,.9), 1000)
X1.wpq = WtProjQuantProfile(X1, c(0,.9), 1000, 100, 0.5, 0.2, 0.9)


plot(X1, pch=19, cex=.1, col=adjustcolor("black", alpha.f=.3))
# normal confidence ellipsoid
require(ellipse); lines(ellipse(sig, level=.9), col="red", lwd=2, lty=2)
lines(X1.pq[,1], X1.pq[,2], lwd=2, col="blue")
lines(X1.wpq[,1], X1.wpq[,2], lwd=2, col="darkgreen")


## Mixture normal simulation
sig2 = matrix(c(1,-.5,-.5,1), nrow=2)
X21 = my.mvrnorm(5e3, mu=c(2,0), Sigma=sig)
X22 = my.mvrnorm(5e3, mu=c(-2,0), Sigma=sig2)
X2 = rbind(X21,X22)
X2.pq = ProjQuant(X2, c(0,.9), 1000)
X2.wpq = WtProjQuantProfile(X2, c(0,.9), 1000, 100, 0.5, 0.2, 0.9)

plot(X2, pch=19, cex=.1)
lines(X2.pq[,1], X2.pq[,2], lwd=2, col="blue")
lines(X2.wpq[,1], X2.wpq[,2], lwd=2, col="darkgreen")

## Mixture normal, assymmetric
X31 = my.mvrnorm(5e3, mu=c(2,0), Sigma=sig)
X32 = my.mvrnorm(2e3, mu=c(-1,0), Sigma=.3*sig2)
X3 = rbind(X31,X32)
X3.pq = ProjQuant(X3, c(0,.9), 1000)

X3.wpq = WtProjQuantProfile(X3, c(0,.9), 1000, 100, 0.5, 0.2, .9)

plot(X3, pch=19, cex=.1)
lines(X3.pq[,1], X2.pq[,2], lwd=2, col="blue")
lines(X3.wpq[,1], X2.wpq[,2], lwd=2, col="darkgreen")

