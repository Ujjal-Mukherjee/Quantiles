rm(list=ls())
set.seed(10092014)
setwd("C:/Study/My projects/Quantiles/Codes/SBM")

library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(ellipse)

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

X = my.mvrnorm(1000, c(0,0), matrix(c(1,.3,.3,1), ncol=2))

u = c(.3,.5)
Xu = X%*%u
Xuperp = sqrt(apply(X^2,1,sum) - Xu^2)
mean(Xu*dnorm(Xuperp))

sourceCpp("ProjQuantNew.cpp", verbose=TRUE, rebuild=FALSE);

pq(X,X[1,])
wpq(X, X[1,], sigma=3)
