rm(list=ls())
set.seed(123)

n<-100;

#x1<-rnorm(n,0,10);
#x1<-x1+min(x1)*0.5;
#x2<-rnorm(n,0,6);

x1<-runif(n,-1,1)
x2<-runif(n,-1,1)

r<-x1^2+x2^2

r_index<-which(r>1)

f<-0.4;
x1<-x1[-r_index]+0.7;
x2<-x2[-r_index]*f;

th<-0*pi/4;

r=function(th)return(cbind(c(cos(th),sin(th)),c(-sin(th),cos(th))));

dat<-cbind(x1,x2);

dat45<-dat%*%r(th);

th<-pi/2;

x2<-x2;
dat<-cbind(x1,x2);

dat135<-dat%*%r(th);

d<-rbind(dat45,dat135);


d<-apply(d,2,FUN=function(x){return((x-mean(x))/max(x))});


plot(d[,1],d[,2],pch=".",col=8)

u<-c(0.0,0.9)

beta=(1+sqrt(t(u)%*%u))/2;
normu=sqrt(t(u)%*%u);

library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)

sourceCpp("ProjQuantNew1.cpp", verbose=TRUE, rebuild=FALSE);


p<-d%*%u;
a<-WECDF(p,seq(0,1,length.out=length(p)),0.2)

h_val=0.25

x<-WtProjQuantProfileMod(d,u,h_val,1000);
lines(x[,1],x[,2],col=2);

step=0.01;

dep.vec<-vector(mode="numeric", length=dim(x)[1]);

for(i in 1:dim(x)[1]){

	dep.vec[i]<-DataDepthECDF(d, x[i,]);

}

#sd_depth = sd(dep.vec);

sd_depth = max(dep.vec)-min(dep.vec);

j=0;

DepVec=function(h,d,u){

	x<-WtProjQuantProfileMod(d,u,h*sd(d),1000);
	lines(x[,1],x[,2],col=1);

	for(i in 1:dim(x)[1]){

		dep.vec[i]<-DataDepthECDF(d, x[i,]);

	}
  
  dep.vec=dep.vec-(1-normu)/4;
  
  norm.dep.sq=(t(dep.vec)%*%dep.vec);

  return(norm.dep.sq);
	
}


p=optimize(f=DepVec, interval=c(0.05,0.3), d=d, u=u)




##################################################################
x<-WtProjQuantProfileMod(d,u,0.2,1000);
lines(x[,1],x[,2],col=3);

x<-WtProjQuantProfileMod(d,u,0.3,1000);
lines(x[,1],x[,2],col=4);

x<-WtProjQuantProfileMod(d,u,0.4,1000);
lines(x[,1],x[,2],col=5);

x<-WtProjQuantProfileMod(d,u,0.6,1000);
lines(x[,1],x[,2],col=6);

x<-WtProjQuantProfileMod(d,u,0.8,1000);
lines(x[,1],x[,2],col=7);

x<-ProjQuant(d,u,1000);
lines(x[,1],x[,2],col=1);

legend("bottomright", legend=c("h=0.1", "h=0.2", "h=0.3", "h=0.4", "h=0.6", "h=0.8", "h=1.0"), col=c(2,3,4,5,6,7,1), lty=c(1,1,1,1,1,1,1));

#=================



plot(d[,1],d[,2],pch=".",col=8, xlim=c(-1.5,1.5))

u<-c(0.0,0.95)


x<-KernelProjQuantProfile(d,u,0.1,1000);
lines(x[,1],x[,2],col=2);

x<-KernelProjQuantProfile(d,u,0.05,1000);
lines(x[,1],x[,2],col=3);

u<-c(0,0.95)
x<-KernelProjQuantProfile(d,u,0.001,1000);
lines(x[,1],x[,2],col=4);

x<-ProjQuant(d,u,1000);
lines(x[,1],x[,2],col=1);

legend("bottomright", legend=c("h=0.1", "h=0.05", "h=0.01", "h=inf"), col=c(2,3,4,1), lty=c(1,1,1,1,1,1,1));

#============================
pdf("QuantileGraphics.pdf")
u<-c(0,0.9);

B<-1000;

x1<-WtProjQuantProfileMod(d,u,0.3,1000);

U<-t(apply(x1,1,FUN=function(x){return(x/sqrt(t(x)%*%x))}));

R<-matrix(nrow=dim(U)[1], ncol=B);

for(i in 1:B){

	d1<-d[sample.int(dim(d)[1],size=ceiling(0.8*dim(d)[1]), replace=TRUE),];
	x1<-WtProjQuantProfileMod(d1,u,0.3,1000);

	R[,i]<-apply(x1,1,FUN=function(x){return(sqrt(t(x)%*%x))});

}

UL<-apply(R,1,FUN=function(x){return(sort(x)[ceiling(length(x)*0.975)])});
LL<-apply(R,1,FUN=function(x){return(sort(x)[floor(length(x)*0.025)])});
ML<-apply(R,1,FUN=function(x){return(mean(x))});


UL<-UL*U;
LL<-LL*U;
ML<-ML*U;

plot(d[,1],d[,2],pch=".",col=8)

lines(UL[,1],UL[,2],col=2);
lines(LL[,1],LL[,2],col=2);
lines(ML[,1],ML[,2],col=4);

#===

u<-c(0,0.9);

B<-1000;

x1<-KernelProjQuantProfile(d,u,0.01,1000);

U<-t(apply(x1,1,FUN=function(x){return(x/sqrt(t(x)%*%x))}));

R<-matrix(nrow=dim(U)[1], ncol=B);

for(i in 1:B){

	d1<-d[sample.int(dim(d)[1],size=ceiling(0.8*dim(d)[1]), replace=TRUE),];
	x1<-KernelProjQuantProfile(d1,u,0.01,1000);

	R[,i]<-apply(x1,1,FUN=function(x){return(sqrt(t(x)%*%x))});

}

UL<-apply(R,1,FUN=function(x){return(sort(x)[ceiling(length(x)*0.975)])});
LL<-apply(R,1,FUN=function(x){return(sort(x)[floor(length(x)*0.025)])});
ML<-apply(R,1,FUN=function(x){return(mean(x))});


UL<-UL*U;
LL<-LL*U;
ML<-ML*U;

plot(d[,1],d[,2],pch=".",col=8)

lines(UL[,1],UL[,2],col=2);
lines(LL[,1],LL[,2],col=2);
lines(ML[,1],ML[,2],col=4);

dev.off();


