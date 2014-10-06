setwd("C:/Study/My projects/Quantiles/Codes")
rm(list=ls());

library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)

sourceCpp("DataDepth.cpp", verbose=TRUE, rebuild=FALSE);
sourceCpp("ProjQuantile.cpp", verbose=TRUE, rebuild=FALSE);

d1<-read.table("FisherIris.csv", sep=",", header=TRUE, as.is=TRUE, fill=TRUE);

library(car)
library(MASS)

scatterplotMatrix(~Sepal.length+Sepal.width+Petal.length+Petal.width|Species, data=d1);
d<-d1[,-c(1)];
CenterScale<-function(d){

	d<-apply(d,2,FUN=function(x){return(x-mean(x))});
	d<-d/max(abs(d));

	return(d);
}

DataDepthClassifier<-function(d,x.new,classCol){

	ClassVec<-d[,which(colnames(d)==classCol)]
	UniqueClass<-unique(ClassVec);
	
	N_Class<-length(UniqueClass);
	
	Class_Depth<-as.data.frame(matrix(nrow=N_Class,ncol=2));
	colnames(Class_Depth)<-c("Class","Depth_Value");

	for(i in 1:N_Class){

		d1<-d[which(ClassVec==UniqueClass[i]),-which(colnames(d)==classCol)];
		d1<-as.matrix(d1);
		d1mean<-apply(d1,2,FUN=function(x){return(median(x))});
		d1max<-max(d1);
		d1<-(d1-rep(1,dim(d1)[1])%*%t(d1mean))/d1max;
		x.new<-(x.new-d1mean)/d1max;
		
		#DepthT = depthx(d1,x.new);
		DepthT = ProjQuantileDepth(d1,x.new,10,0.5,0.2,0.1);
		Class_Depth[i,1] = as.character(UniqueClass[i]);
		Class_Depth[i,2] = DepthT	
	}
	
	return(Class_Depth);

}

for(j in 1:100){

Sel<-sample.int(dim(d)[1], floor(0.8*dim(d)[1]));
dtrain<-d[Sel,];
dtest<-d[-Sel,];

Compare<-as.data.frame(matrix(nrow=dim(dtest)[1],ncol=6));
colnames(Compare)<-c("Original","Predicted",unique(d$Species)[1],unique(d$Species)[2],unique(d$Species)[3],"result");


for(i in 1:dim(dtest)[1]){

n.vec<-dtest[i,-which(colnames(d)=="Species")];
n.vec<-as.numeric(n.vec);

Out<-DataDepthClassifier(dtrain,n.vec,"Species");

Compare[i,1]<-dtest[i,which(colnames(d)=="Species")];
Compare[i,2]<-Out[which(Out[,2]==max(Out[,2])),1];
Compare[i,3]<-Out[which(Out[,1]==unique(d$Species)[1]),2];
Compare[i,4]<-Out[which(Out[,1]==unique(d$Species)[2]),2];
Compare[i,5]<-Out[which(Out[,1]==unique(d$Species)[3]),2];
if(Compare[i,1]==Compare[i,2]){
	Compare[i,6]=1;
}else{
	Compare[i,6]=0;
}
}

Frac_Correct<-sum(Compare$result)/length(Compare$result);

if(j==1)Accuracy<-Frac_Correct else Accuracy<-c(Accuracy, Frac_Correct);

}

plot(density(Accuracy, bw=0.03));

print(mean(Accuracy));


#######################################################



d<-read.table("Wines.csv", sep=",", header=TRUE, as.is=TRUE, fill=TRUE);

library(car)
library(MASS)

#scatterplotMatrix(~Alcohol+Malic_Acid+Ash+Alcalinity+Magnesium
#			+Phenols+Flavanoids+Non_Phenol+Proanthocianins+
#			Color+Hue+Dilution+Proline|Class, data=d);

#scatterplotMatrix(~d[,1:34]|d[,35])

for(j in 1:100){

Sel<-sample.int(dim(d)[1], floor(0.8*dim(d)[1]));
dtrain<-d[Sel,];
dtest<-d[-Sel,];

Compare<-as.data.frame(matrix(nrow=dim(dtest)[1],ncol=6));
colnames(Compare)<-c("Original","Predicted",unique(d$Class)[1],unique(d$Class)[2],unique(d$Class)[3],"result");


for(i in 1:dim(dtest)[1]){

n.vec<-dtest[i,-which(colnames(d)=="Class")];
n.vec<-as.numeric(n.vec);

Out<-DataDepthClassifier(dtrain,n.vec,"Class");

Compare[i,1]<-dtest[i,which(colnames(d)=="Class")];
Compare[i,2]<-Out[which(Out[,2]==max(Out[,2])),1];
Compare[i,3]<-Out[which(Out[,1]==unique(d$Class)[1]),2];
Compare[i,4]<-Out[which(Out[,1]==unique(d$Class)[2]),2];
Compare[i,5]<-Out[which(Out[,1]==unique(d$Class)[3]),2];
if(Compare[i,1]==Compare[i,2]){
	Compare[i,6]=1;
}else{
	Compare[i,6]=0;
}
}

Frac_Correct<-sum(Compare$result)/length(Compare$result);

if(j==1)Accuracy<-Frac_Correct else Accuracy<-c(Accuracy, Frac_Correct);

}

plot(density(Accuracy, bw=0.03));

print(mean(Accuracy));
