## Initial simulations of Projection quantile depth and comparisions with projection depth

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
wEPQD = function(X, sig, mingrid, maxgrid, res=.2, nu=1e3){
  
  # find scaling mean and cov
  b = colMeans(X)
  AtA = cov(X)
  AtA.inv = solve(cov(X))
  X0 = (X-b) %*% AtA.inv
    
  # make grid of points
  pts = seq(mingrid, maxgrid, by=res)
  lengrid = length(pts)
  xcoord = rep(pts, rep(lengrid,lengrid))
  ycoord = rep(pts, lengrid)
  xygrid = cbind(xcoord,ycoord)
  rm(xcoord,ycoord)
  grid0 = (xygrid - b) %*% AtA.inv
    
  ## get matrix of weighted PQDs for all points
  npt = dim(xygrid)[1]
  Fuxu.mat = matrix(0, nrow=npt, ncol=nu)
  
  normsq.X0 = apply(X0^2, 1, sum)
  normsq.grid0 = apply(grid0^2, 1, sum)
  
  # loop over nu pts on unit circle then take max
  for(iu in 1:nu){
    u = as.matrix(rnorm(2)); u = u/sqrt(sum(u^2))
    I.minus.Pu = diag(2) - u%*%t(u)
    perp.mean = I.minus.Pu %*% b
    perp.cov = I.minus.Pu %*% AtA %*% I.minus.Pu
    
    Xuperp = X %*% I.minus.Pu
    bperp = b %*% I.minus.Pu
    scaled.perp = sqrt(Xuperp^2 %*% ones(ncol(X),1))
    w = dnorm(scaled.perp, sd=sig)
    #w = dcauchy(Xuperp, scale=sig)
    uecdf = ecdf(w * (X%*%u))
    
    gridperp = xygrid %*% I.minus.Pu
    scaled.gridperp = sqrt(gridperp^2 %*% ones(ncol(X),1))
    wu = dnorm(scaled.gridperp, sd=sig)
    #wu = dcauchy(sqrt(apply(xygrid^2,1,sum) - xygrid.u^2), scale=sig)
    Fuxu.mat[,iu] = uecdf(wu * (xygrid%*%u))
  }
  EPQD.vec = 1/(1+apply(abs(Fuxu.mat-.5), 1, max))
  
  # check if contains origin... if so, assign NA to that depth
  which0 = which(xygrid[,1]==0 & xygrid[,2]==0)
  if(length(which0>0)){
    EPQD.vec[which0] = NA
  } 
  
  ## plot result
  par(mfrow=c(1,2))
  persp(pts, pts, matrix(EPQD.vec, nrow=lengrid, byrow=T),
        main="Projection Quantile Depth",
        xlab="x1", ylab="x2", zlab="PQD(x,F)",
        theta=-45, phi=45)
  
  # contour plot
  z = contour(pts, pts, matrix(EPQD.vec, nrow=lengrid, byrow=T),
              lwd=2, col="red", nlevels=20)
  points(X, pch=19, cex=.2)
  par(mfrow=c(1,1))
  
  return(EPQD.vec)
  
}

## Empirically calculate PQD with grid search
# Bivariate normal mixture
set.seed(120214)
sig = matrix(c(1,.9,.9,1), nrow=2)
sig2 = matrix(c(1,-.9,-.9,1), nrow=2)
X1 = my.mvrnorm(500, mu=c(2,7), Sigma=sig)
X2 = my.mvrnorm(500, mu=c(6,2), Sigma=sig)
X = rbind(X1,X2)

# z = WtProjQuantProfileMod(d,c(0,.9),.1, 1000)
# plot(d)
# lines(z)

k = wEPQD(X, sig=1, mingrid=-2, maxgrid=10)

# naive PQD
sig=.5

npt = dim(xygrid)[1]
EPQD.vec = rep(0, npt)
normsq.X = apply(X^2,1,sum)
normsq.grid = apply(xygrid^2,1,sum)

for(i in 1:npt){
  u = xygrid[i,]; u = u/sqrt(sum(u^2))
  Xu = X%*%u
  Xuperp = sqrt(normsq.X - Xu^2)
  w = dnorm(Xuperp, sd=sig)
  #w = dcauchy(Xuperp, scale=sig)
  uecdf = ecdf(w*Xu)
  
  igrid.u = xygrid[i,]%*%u
  iwu = dnorm(0, sd=sig)
  #wu = dcauchy(sqrt(apply(xygrid^2,1,sum) - xygrid.u^2), scale=sig)
  EPQD.vec[i] = 1/(1+abs(uecdf(iwu*igrid.u)-.5))
}

par(mfrow=c(1,2))
persp(pts, pts, matrix(EPQD.vec, nrow=lengrid, byrow=T),
      main="Naive Projection Quantile Depth",
      xlab="x1", ylab="x2", zlab="PQD(x,F)",
      theta=-45, phi=45)

# contour plot
z = contour(pts, pts, matrix(EPQD.vec, nrow=lengrid, byrow=T),
            lwd=2, col="red", nlevels=20)
points(X, pch=19, cex=.2)
par(mfrow=c(1,1))