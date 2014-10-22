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

## function for PQ and WPQ plots for given data
PwPplot = function(dat, isEllipse=F, sig=NULL, t1, t2){
  
  par(mfrow=c(1,2))
  plot(dat, pch=19, cex=.1, col=adjustcolor("black", alpha.f=.1),
       main=t1, xlab="", ylab="")
  
  for(i in 5:9){
    iq = c(0,i/10)
    dat.pq = ProjQuant(dat, iq, 1000)
    lines(dat.pq[,1], dat.pq[,2], lty=2, col="red")
  }
  
  # WPQ
  plot(dat, pch=19, cex=.1, col=adjustcolor("black", alpha.f=.1),
       main=t2, xlab="", ylab="")
  for(i in 5:9){
    iq = c(0,i/10)
    dat.wpq = WtProjQuantProfile(dat, iq, 1000, 100, 0.5, 0.1, 0.9)
    if(isEllipse)
      lines(ellipse(sig, level=i/10), col="blue")
    lines(dat.wpq[,1], dat.wpq[,2], lty=2, col="red")
  }
  par(mfrow=c(1,1))
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

## bivariate exponential
require(VineCopula)
X1 = rexp(n)
X2 = rexp(n)
c12 = BiCopSim(n, 1, par=.5)
X4 = cbind(X1,X2)*c12


PwPplot(X1, isEllipse = T, sig=sig, t1="a.",t2="b.")
PwPplot(X2, t1="c.",t2="d.")
PwPplot(X3, t1="e.",t2="f.")
PwPplot(X4, t1="g.",t2="h.")

## Compare bivariate normal and mixture normal plot... AIStats paper
X1.pq = ProjQuant(X1, c(0,.9), 1000)
X1.wpq = WtProjQuantProfile(X1, c(0,.9), 1000, 100, 0.5, 0.1, 0.9)

plot(X1, pch=19, cex=.1, col=adjustcolor("black", alpha.f=.3),
     main="Bivariate Normal", xlab="X1", ylab="X2")
lines(X1.pq[,1], X1.pq[,2], lwd=2, col="blue")
lines(X1.wpq[,1], X1.wpq[,2], lwd=2, col="red")
lines(ellipse(sig, level=.9), lwd=2)

legend("bottomright", c("Normal ellipsoid","PQ","DCW"),
       lty=1, lwd=2, cex=.7,
       col=c("black","blue","red"))

# BVN mixture
X41 = my.mvrnorm(5e3, mu=c(2,0), Sigma=sig)
X42 = my.mvrnorm(2e3, mu=c(-1,0), Sigma=.5*sig2)
X4 = rbind(X41,X42)

X4.pq = ProjQuant(X4, c(0,.9), 1000)
X4.wpq = WtProjQuantProfile(X4, c(0,.9), 1000, 100, 0.5, 0.1, 0.9)

plot(X4, pch=19, cex=.1, col=adjustcolor("black", alpha.f=.3),
     main="Bivariate Normal mixture", xlab="X1", ylab="X2")
# PQ contours

lines(ellipse(cov(X4), centre=colMeans(X4), level=.9), lwd=2)
lines(X4.pq[,1], X4.pq[,2], lwd=2, col="blue")
lines(X4.wpq[,1], X4.wpq[,2], lwd=2, col="red")

legend("bottomright", c("Normal ellipsoid","PQ","DCW"),
       lty=1, lwd=2, cex=.7,
       col=c("black","blue","red"))

