rm(list=ls());

#############################################################################
#Vector Norm  
#############################################################################

Vector.Norm=function(u){

	return(sqrt(t(u)%*%u)[1]);

			}
#############################################################################
#Generalized Quantile function
#############################################################################

GenQuant=function(X,q,u,L){

	k=dim(X)[2];
	n=dim(X)[1];

	U=u/Vector.Norm(u);
	
	Xu=(X%*%U);
	qu=(q%*%U);

	XuU=(X%*%U)%*%t(U);
	quU=(q%*%U)%*%t(U);

	Xu.p=X-XuU;
	qu.p=q-quU;

	b=Vector.Norm(u);

	Xu_qu=as.vector(Xu-rep(qu,n));
	
	Xu.p_qu=diag((Xu.p-rep(qu,n)%*%t(rep(1,k)))%*%t(Xu.p-rep(qu,n)%*%t(rep(1,k))));
				
	Chi_i<-abs(Xu_qu)*(rep(1,n)+L*((Xu_qu)^(-2))*(Xu.p_qu))+b*(Xu_qu);

	Val<-t(rep(1,n))%*%Chi_i;

	return(Val)

	}

#############################################################################
#Optimization Routine - gradient descent
#############################################################################

QuantOptim1<-function(X,q.init,u,L,Del=0.01){

	diff=1	

	k<-dim(X)[2];	

	p=1;

	while(diff>(0.01^k)){

	fni<-GenQuant(X,q.init,u,L);

		for(i in 1:k){

			q.temp.1<-q.init;

			q.temp.2<-q.init;
	
			q.temp.1[i]<-q.temp.1[i]+Del;

			q.temp.2[i]<-q.temp.2[i]-Del;			

			fn.temp.1<-GenQuant(X,q.temp.1,u,L);

			fn.temp.2<-GenQuant(X,q.temp.2,u,L);

			if(!is.na(fn.temp.1) && !is.na(fni)){ if(fn.temp.1<=fni) q.init[i]<-q.temp.1[i]};

			if(!is.na(fn.temp.2) && !is.na(fni)){if(fn.temp.2<=fni) q.init[i]<-q.temp.1[i]};

			}

	fnf<-GenQuant(X,q.init,u,L);

	if(!is.na(fni) && !is.na(fnf))diff<-abs(fni-fnf);

	p=p+1;

	if(p>10000) break;

	}
	
return(q.init);

}

#############################################################################
#Gradient vector
#############################################################################

GradVect<-function(X,q,u,L,Del=0.001){

	k<-dim(X)[2];

	GradVect<-vector(mode="numeric", length=k);

	for(i in 1:k){

		tempL<-q;
		tempH<-q;

		tempL[i]<-tempL[i]-Del;
		tempH[i]<-tempH[i]+Del;

		GradVect[i]<-(GenQuant(X,tempH,u,L)-GenQuant(X,tempL,u,L))/(2*Del);

		}

	return(GradVect);

}
#############################################################################

QuantOptim2<-function(X,u,L,Del=0.01){

	X.norm.vector<-diag(X%*%t(X));

	R.Max<-max(X.norm.vector);

	N<-1000;

	Vector.Grid<-seq(-R.Max, R.Max, length.out=N);

	Grid.Mat<-Vector.Grid%*%t(u/Vector.Norm(u));

	Quant.Value.Vec<-vector(mode="numeric", length=N);

	for(i in 1:N){

	Quant.Value.Vec[i]<-GenQuant(X,Grid.Mat[i,],u,L);

	}

	return(Grid.Mat[grep(min(Quant.Value.Vec), Quant.Value.Vec)[1],]);
	
}

#############################################################################

#############################################################################
#Projection Kernel  
#############################################################################

Projection.Kernel.Quantile=function(X,u){

	#Generate the unit vector in the direction of u

	U=u/Vector.Norm(u);

	Xu=(X%*%U);

	qu=sort(Xu)[floor(length(Xu)*((1+Vector.Norm(u)[1])/2))];

	Qu.Proj=qu*U;

	return(Qu.Proj);

	}


#############################################################################
#Generate Random Vector
#############################################################################

Rand.Vector=function(p,N,u){

	b=Vector.Norm(u);

	Zp=matrix(nrow=N, ncol=p);

	for(i in 1:p){

	Zp[,i]<-rnorm(N,0,1);

	}

	for(i in 1:N){

	Zp[i,]<-Zp[i,]/Vector.Norm(Zp[i,]);

	}

	u.P=b[1]*Zp;

	return(u.P);

	}

#############################################################################
X1<-rnorm(1000,0,1);
X2<-rnorm(1000,10,10);
X3<-rnorm(1000,0,5);

X<-as.matrix(cbind(X1,X2,X3));

#X<-as.matrix(cbind(X1,X2));

q<-c(1,10,5);

u<-c(0.2,0.3,0.2);

L=1;	

N=1000;


p=dim(X)[2];

u.P=Rand.Vector(p,N,u);

Quant.Mat<-matrix(nrow=N, ncol=p);

for(i in 1:N){

	Quant.Mat[i,]=QuantOptim2(X,u.P[i,],L);

		}

library(scatterplot3d);

scatterplot3d(X1,X2,X3,type="p",highlight.3d=FALSE,xlab="X1",ylab="X2",zlab="X3",main="Quantile Surface",color=ceiling(4), pch="*",
			xlim=c(min(X1),max(X1)),
			ylim=c(min(X2),max(X2)),
			zlim=c(min(X3),max(X3)));

par(new=TRUE);


scatterplot3d(Quant.Mat[,1],Quant.Mat[,2],Quant.Mat[,3],type="p",highlight.3d=FALSE,xlab="X1",ylab="X2",zlab="X3",main="Quantile Surface",color=ceiling(2), pch="*",
			xlim=c(min(X1),max(X1)),
			ylim=c(min(X2),max(X2)),
			zlim=c(min(X3),max(X3)));



scatterplot3d(Quant.Mat[,1],Quant.Mat[,2],Quant.Mat[,3],type="p",highlight.3d=FALSE,xlab="X1",ylab="X2",zlab="X3",main="Quantile Surface",color=ceiling(2), pch=".");




