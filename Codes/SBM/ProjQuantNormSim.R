rm(list=ls())
set.seed(10092014)
setwd("C:/Study/My projects/Quantiles/Codes/SBM")

library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(ellipse)

## cpp source file
sourceCpp("ProjQuantNew.cpp", verbose=TRUE, rebuild=FALSE);

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
X1 = my.mvrnorm(5e3, mu=c(0,0), Sigma=sig)

## Mixture normal simulation
sig2 = matrix(c(1,-.5,-.5,1), nrow=2)
X21 = my.mvrnorm(5e3, mu=c(1,0), Sigma=sig)
X22 = my.mvrnorm(5e3, mu=c(-1,0), Sigma=sig2)
X2 = rbind(X21,X22)

## Mixture normal, assymmetric
X31 = my.mvrnorm(5e3, mu=c(1,0), Sigma=sig)
X32 = my.mvrnorm(2e3, mu=c(-1,0), Sigma=.3*sig2)
X3 = rbind(X31,X32)

# BVN mixture
X41 = my.mvrnorm(5e3, mu=c(2,0), Sigma=sig)
X42 = my.mvrnorm(2e3, mu=c(-1,0), Sigma=.5*sig2)
X4 = rbind(X41,X42)

X4.pq = ProjQuant(X2, c(0,.5), 1000)
X4.pq = WtProjQuantProfile(X2, c(0,.9), 1000, 100, 0.5, 0.1, 0.9)
X4.wpqnorm = ProjQuantNorm(X2, c(0,.9), 2.3, 1000)

plot(X2, pch=19, cex=.1, col=adjustcolor("black", alpha.f=.3),
     main="Bivariate Normal mixture", xlab="X1", ylab="X2")

lines(X4.pq)
lines(X4.wpqnorm, col="red")
