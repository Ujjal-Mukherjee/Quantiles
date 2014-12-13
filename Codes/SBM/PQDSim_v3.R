## PQDSim: Initial simulations of Projection quantile depth and comparisions with projection depth
## v3: classification scenario for high-dimensional data

## Functions
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

ones = function(m,n){
  matrix(1, nrow=m, ncol=n)
}

## function to calculate weighted projection quantile depth
wEPQD1 = function(X, grid, sig, nu=1e3){
  
  p = ncol(X)
  b = apply(X, 2, median)
  X0 = X - ones(nrow(X),1) %*% b
  grid0 = grid - ones(nrow(grid),1) %*% b
    
  ## get matrix of weighted PQDs for all points
  npt = dim(grid)[1]
  Fuxu.mat = matrix(0, nrow=npt, ncol=nu)
  
  # loop over nu pts on unit circle then take max
  for(iu in 1:nu){
    u = as.matrix(rnorm(p)); u = u/sqrt(sum(u^2))
    I.minus.Pu = diag(p) - u%*%t(u)
    
    Xuperp = X0 %*% I.minus.Pu
    scaled.perp = sqrt(Xuperp^2 %*% ones(ncol(X),1))
    #w = ifelse(scaled.perp>sig, 0, 1)
    #w = sig*exp(-scaled.perp/sig)
    w = dnorm(scaled.perp, sd=sig)
    #w = dcauchy(Xuperp, scale=sig)
    uecdf = ecdf(w * (X0%*%u))
    
    gridperp = grid0 %*% I.minus.Pu
    scaled.gridperp = sqrt(gridperp^2 %*% ones(ncol(X),1))
    #wu = ifelse(scaled.gridperp>sig, 0, 1)
    #wu = sig*exp(-scaled.gridperp/sig)
    wu = dnorm(scaled.gridperp, sd=sig)
    #wu = dcauchy(sqrt(apply(xygrid^2,1,sum) - xygrid.u^2), scale=sig)
    Fuxu.mat[,iu] = uecdf(wu * (grid0%*%u))
  }
  EPQD.vec = 1/(1+apply(abs(Fuxu.mat-.5), 1, max))
  
  return(cbind(grid,EPQD.vec))
  
}

## Empirically calculate PQD with grid search
# Bivariate normal mixture
p = 2000
set.seed(120214)
X1 = my.mvrnorm(50, mu=rep(0,p), Sigma=diag(p))
X2 = my.mvrnorm(50, mu=rep(3,p), Sigma=diag(p))
Xa = rbind(X1,X2)
Xb = my.mvrnorm(50, mu=rep(6,p), Sigma=diag(p))
X = rbind(Xa,Xb)

## find out 2 depths
col.vec = c(rep("red",100), rep("green",50))

d1 = wEPQD1(Xa, X, sig=1)
plot(d1[,101], col=col.vec, pch=19)

system.time(d2 <- wEPQD1(Xb, X, sig=1))
plot(d2[,101], col=col.vec, pch=19)
