rm(list=ls())
set.seed(123)
n<-10000;

library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)

sourceCpp("ProjQuant.cpp", verbose=TRUE, rebuild=FALSE);

x1<-rnorm(n,0,10);
x1<-x1+min(x1)*0.5;
x2<-rnorm(n,0,6);

th<-pi/4;

r=function(th)return(cbind(c(cos(th),sin(th)),c(-sin(th),cos(th))));

dat<-cbind(x1,x2);

dat45<-dat%*%r(th);

th<-3*pi/4;

dat135<-dat%*%r(th);

d<-rbind(dat45,dat135);


d<-apply(d,2,FUN=function(x){return((x-mean(x))/max(x))});


plot(d[,1],d[,2],pch=".",)


u=c(0,0.9);

x<-ProjQuant(d,u,1000);


lines(x[,1],x[,2],col=3);

y<-WtProjQuantProfile(d,u,1000,100,0.5,0.2,0.9);

lines(y[,1],y[,2],col=4);

r<-apply(y,1,FUN=function(x){return(sqrt(x[1]^2+x[2]^2))})

theta<-seq(0,2*pi,length.out=dim(y)[1]);

m<-ksmooth(theta,r,kernel="normal", bandwidth=(0.1*1.06*sd(r)^(-0.2)));

xx<-m$y*cos(m$x);
yy<-m$y*sin(m$x);


lines(xx,yy, col=2, pch=".");

x1<-rnorm(n,0,1);
x2<-rnorm(n,0,5);

d<-cbind(x1,x2);

z<-vector(mode="numeric", length=dim(d)[1]);

for(i in 1:dim(d)[1]){

	z[i]<-ProjQuantileDepth(d,d[i,],5,0.5,0.2,0.8);	

}



library(scatterplot3d)

library(rgl)


plot3d(d[,1], d[,2], z, col="red", size=3);


