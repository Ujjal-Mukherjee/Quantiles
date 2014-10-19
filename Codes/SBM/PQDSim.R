## Initial simulations of Projection quantile depth and comparisions with projection depth

## Functions

# define grid of points
pts = seq(-1,1,by=.1)
lengrid = length(pts)
xcoord = rep(pts, rep(lengrid,lengrid))
ycoord = rep(pts, lengrid)
xygrid = cbind(xcoord,ycoord)
rm(xcoord,ycoord)

# Bivariate standard normal
# projection depth
norm.vec = sqrt(xygrid[,1]^2+xygrid[,2]^2)
cN = qnorm(.75)
pd.vec = cN/(cN+norm.vec)
# projection quantile depth
pqd.vec = 2/(1+2*pnorm(norm.vec))

par(mfrow=c(1,2))
persp(pts, pts, matrix(pd.vec, nrow=lengrid, byrow=T),
      main="Projection Depth",
      xlab="x1", ylab="x2", zlab="PD(x,F)",
      theta=-45, phi=45)
persp(pts, pts, matrix(pqd.vec, nrow=lengrid, byrow=T),
      main="Projection Quantile Depth",
      xlab="x1", ylab="x2", zlab="PQD(x,F)",
      theta=-45, phi=45)
par(mfrow=c(1,1))

# comparing PD and PQD
plot(pd.vec, pqd.vec, type="l", lwd=2,
     main="Comparison of PD and PQD", xlab="PD", ylab="PQD")
