rm(list=ls())
set.seed(123)
n<-1000;


library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)

sourceCpp("DataDepth.cpp", verbose=TRUE, rebuild=FALSE);

#x1<-rnorm(n,0,1);
#x2<-rnorm(n,0,5);

#d<-cbind(x1,x2);


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




x<-c(0,0);
e<-c(0,0);

p<-QuadrantFrac(d,x,e);

print(p);

q<-ecdf(d,x);

print(q);

k<-100;

x3<-seq(min(d[,1]), max(d[,1]),length.out=k);
x4<-seq(min(d[,2]), max(d[,2]),length.out=k);

xx<-expand.grid(x=x3,y=x4);
xx<-as.matrix(xx);

z<-ecdfGrid(d,xx);

y<-matrix(nrow=k, ncol=k);

for(i in 1:k){

	for(j in 1:k){

		x<-c(x3[i],x4[j]);

		y[i,j]<-ecdf(d,x);
	}
}

res<-persp(x3, x4, y, theta=-45, phi=45, shade=TRUE, col=5);

p<-c(0,0);

dep<-depthx(d,p);

med<-medianp(d);

library(scatterplot3d)

library(rgl)

dep1<-DepthGrid(d);

scatterplot3d(d[,1],d[,2],dep1, pch=".", highlight.3d=TRUE)

plot3d(d[,1], d[,2], dep1, col="red", size=3);


y1<-matrix(nrow=k, ncol=k);

for(i in 1:k){

	for(j in 1:k){

		x<-c(x3[i],x4[j]);

		y1[i,j]<-depthx(d,x);
	}
}

res<-persp(x3, x4, y1, theta=-145, phi=95, shade=TRUE, col=5);


d<-read.table("FisherIris.csv", sep=",", header=TRUE, as.is=TRUE, fill=TRUE);

library(car)
library(MASS)

scatterplotMatrix(~Sepal.length+Sepal.width+Petal.length+Petal.width|Species, data=d);



