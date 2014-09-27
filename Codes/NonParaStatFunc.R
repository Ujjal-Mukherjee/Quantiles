#JOINT PROBABILITY DENSITY SURFACE
#=================================
#STEP 1: CREATE DATA
#=================================

rm(list=ls());
#f1<-read.table("SeaTemp.csv",sep=",",header=TRUE);
#Temp<-f1$pTemperature; #Sea Temperature;
#Time<-f1$T; #Time in a linear scale;


set.seed(4145211);
#Num<-ceiling(sqrt(45211));
Num<-256;
Time<-runif(Num,-1.25,1.25);

Error<-rnorm(Num,0,4);

Temp<-exp(3*Time)+5*sin((30*Time)/(2*pi))+Error;

#=================================
#STEP 2: NORMALIZE THE DATA
#=================================

N<-length(Temp); #Data length

ii<-rep(1,N); #Vector of 1's N times
Temp_N<-(1/sd(Temp))*(Temp-(1/N)*t(ii)%*%Temp); #Normalized Temp
Time_N<-(1/sd(Time))*(Time-(1/N)*t(ii)%*%Time); #Normalized Time

Temp_Min<-floor(min(Temp_N));
Time_Min<-floor(min(Time_N));
Temp_Max<-ceiling(max(Temp_N));
Time_Max<-ceiling(max(Time_N));

#=================================
#STEP 3:X RESOLUTION (STEP LENGTH)
#=================================

Resolution<-0.1;
X_N<-max(((Temp_Max-Temp_Min)/Resolution),((Time_Max-Time_Min)/Resolution));

Temp_Resolution<-(Temp_Max-Temp_Min)/(X_N-1);

Time_Resolution<-(Time_Max-Time_Min)/(X_N-1);

X_Temp_N<-as.vector(seq(Temp_Min,Temp_Max,by=Temp_Resolution));

X_Time_N<-as.vector(seq(Time_Min,Time_Max,by=Time_Resolution));

X_Density_N<-mat.or.vec(X_N,X_N);

#=================================
#STEP 3:KERNEL DENSITY ESTIMATION
#=================================

h<-1.06*(N^(1/5)); #Kernel Resolution

iii<-rep(1,N);

Gaussian_Const=1/sqrt(2*pi); 

for(i in 1:X_N){

	for(j in 1:X_N){

		X_Density_N[i,j]<-t(iii)%*%((1/(N*(h^2)))*(Gaussian_Const^2)*exp(((X_Temp_N[i]-Temp_N)^2)*(0.5/h))*exp(((X_Time_N[j]-Time_N)^2)*(0.5/h)));

	}
}

#=================================
#STEP 3:KERNEL DENSITY PLOT
#=================================

library(rgl);

persp3d(X_Temp_N,X_Time_N,X_Density_N,color=3,alpha=0.9,back="lines");

print("HIT <ENTER> TO CONTINUE ...");
readline();
print("WAIT....PROCESSING");

#=================================
#****    LINEAR REGRESSION    ****
#=================================
#     SELECTION OF BEST LAMBDA
#=================================

Y<-Temp+abs(min(Temp))*1.03; #Response Variable with axis shifted by 2 units 
           # for making all responses positive
#We will transform the Y by a simple power transform 
#from lambda=-3 to +3 in jumps of 0.1. 
#Min RSS would be used as a criterion for model selection

Lambda<-seq(-3,3,by=0.1);

Y_Mat<-mat.or.vec(N,length(Lambda));

for(i in 1:length(Lambda)){

if(Lambda[i]==0)
	Y_Mat[,i]<-log(Y) 
else 
	Y_Mat[,i]<-Y^Lambda[i];

}

X<-as.matrix(cbind(rep(1,N),Time)); #Predictor Variable

Y_Hat<-X%*%solve(t(X)%*%X)%*%t(X)%*%Y_Mat; #Predicted Y values

Residual_Mat<-Y_Mat-Y_Hat; #Residual Matrix

RSS_Mat<-t(iii)%*%(Residual_Mat^2);

SE_Mat<-sqrt((1/N)*RSS_Mat);

for(i in 1:length(Lambda)){

if(SE_Mat[i]==min(SE_Mat))
	Min_Lambda=Lambda[i];
}

plot(Lambda,SE_Mat,type="l",col="blue",main="Standard Error plot against Lambda values of Power Transformation of response",xlab="Lambda",ylab="Standard Error");

print("HIT <ENTER> TO CONTINUE ...");
readline();
print("WAIT....PROCESSING");

#=================================
#CREATION OF SUMMARY MATRIX 
#=================================
#CREATION OF ANOVA MATRIX
#=================================

#Power transformation of response with best Lambda

if(Min_Lambda==0)	Y_T<-log(Y) else 	Y_T<-Y^Min_Lambda; 

B<-solve(t(X)%*%X)%*%t(X)%*%Y_T; #Parameter vector

Standard_Error<-sqrt((1/(N-2))*t(Y_T-X%*%B)%*%(Y_T-X%*%B));

SE_B<-Standard_Error*sqrt(diag(solve(t(X)%*%X)));

t_Val<-B/SE_B;

p_Val<-2*(1-pt(t_Val,(N-2)));

Signif<-mat.or.vec(2,1);

for(i in 1:2){

if(p_Val[i]>=0.1)
	Signif[i]<-"NS "
else if(p_Val[i]>=0.05)
	Signif[i]<-"   "
else if(p_Val[i]>=0.01)
	Signif[i]<-"*  "
else if(p_Val[i]>=0.001)
	Signif[i]<-"** "
else
	Signif[i]<-"***"

}

R_Sq<-1-(((N-2)*(Standard_Error^2))/(t(iii)%*%((Y_T-(1/N)*t(iii)%*%Y_T)*(Y_T-(1/N)*t(iii)%*%Y_T))));

Linear_Reg_Summary<-data.frame(B,SE_B,t_Val,p_Val,Signif);

Residual_SS<-(N-2)*(Standard_Error^2);
Total_SS<-t(iii)%*%((Y_T-(1/N)*t(iii)%*%Y_T)*(Y_T-(1/N)*t(iii)%*%Y_T));
Reg_SS<-Total_SS-Residual_SS;

Reg_df<-1;
Res_df<-N-2;
Total_df<-N-1;

SS<-c(Reg_SS,Residual_SS,Total_SS);
df<-c(Reg_df,Res_df,Total_df);
MSS<-SS/df;

F_Val1<-MSS[1]/MSS[2];
p_F_Val1<-1-pf(F_Val1,1,N-2);

F_Val<-c(round(F_Val1,digits=8)," "," ");
p_F_Val<-c(round(p_F_Val1,digits=8)," "," ");


if(p_F_Val1>=0.1) S_F<-"NS " else if(p_F_Val1>=0.05) S_F<-"   " else if(p_F_Val1>=0.01) S_F<-"*  " else if(p_F_Val1>=0.001) S_F<-"** " else S_F<-"***";

Signif_F<-c(S_F," "," ");

Anova_Par<-c("Regression","Residual","Total");

Anova_Reg<-data.frame(Anova_Par,df,SS,MSS,F_Val,p_F_Val,Signif_F);

#=============================
#REGRESSION OUTPUT
#=============================
print("====================================================================");
print("LINEAR REGRESSION OUTPUT: SUMMARY");
print("--------------------------------------------------------------------");
cat("Linear_Reg_Equation :"," E(Temp^",Min_Lambda,"|Time)=",B[[1]],"+",B[[2]],"*Time\n");
print("--------------------------------------------------------------------");
print(Linear_Reg_Summary);
print("--------------------------------------------------------------------");
print("Signif: '***'<0.001, '**'<0.01, '*'<0.05");
cat("Standard Error is ",Standard_Error," with df = ",N-2,"\n");
cat("R Square = ",R_Sq,"\n");
print("--------------------------------------------------------------------");
print("ANOVA TABLE FOR THE REGRESSION");
print("--------------------------------------------------------------------");
print(Anova_Reg);
print("====================================================================");

Temp_Hat<-((as.matrix(cbind(rep(1,N),sort(Time))))%*%B)^(1/Min_Lambda)-2;

plot(Time,Temp,main="Scatter Plot of the response and predictor",xlab="Time",ylab="Temperature");
lines(sort(Time),Temp_Hat,col="blue");

print("HIT <ENTER> TO CONTINUE ...");
readline();
print("WAIT....PROCESSING");

#===============================
#  NON PARAMETRIC REGRESSION
#===============================
#  NADARAYA WATSON ESTIMATOR
#===============================
#      GAUSSIAN KERNEL
#===============================

X_Max=ceiling(max(Time));

X_Min=floor(min(Time));

X_Span=500;

X_Interval=(X_Max-X_Min)/(X_Span-1);

X_Seq=seq(X_Min,X_Max,by=X_Interval);

One_1<-rep(1,N);
One_2<-rep(1,X_Span);

H_NW<-1.06*sd(Time)*(N^(-0.2)); #Nadaraya Watson Resolution

NW_Estimation<-((exp(-1*((X_Seq%*%t(One_1)-One_2%*%t(Time))/H_NW)^2))%*%Temp)/((exp(-1*((X_Seq%*%t(One_1)-One_2%*%t(Time))/H_NW)^2))%*%(One_1));

plot(Time,Temp,main="Temp vs. Time", xlab="Time", ylab="Temp"); #Point plot of time and temp
abline(lm(Temp~Time),col="red",lty=2); #Simple Linear Plot
lines(X_Seq,NW_Estimation,col="black"); #Non Parametric regression plot

print("HIT <ENTER> TO CONTINUE ...");
readline();
print("WAIT....PROCESSING");

#==============================
#   LOCAL CURVE FITTING
#==============================

LocalCurveFitting=function(P){

#print("How many degrees of Polynomial should be fitted (Min: 1, Max: 4)???");

#P<-floor(as.numeric(readline()));

Local_Pol_Estimate<-mat.or.vec(X_Span,1);

Xx<-mat.or.vec(N,P+1);



for(i in 1:X_Span){

	for(j in 1:(P+1)){

		Xx[,j]<-(Time-X_Seq[i])^(j-1)/factorial(j-1);

			}

Wx<-diag(Gaussian_Const*exp((-1)*(Time-X_Seq[i])^2/(H_NW^2)));

Rnx<-solve(t(Xx)%*%Wx%*%Xx)%*%t(Xx)%*%Wx%*%Temp;

Local_Pol_Estimate[i]<-Rnx[1];

}

lines(X_Seq,Local_Pol_Estimate,col=P+1); #Local Polynomial regression plot

rm(P,Xx,Wx,Rnx,Local_Pol_Estimate);
}

#print("How many Local Polynomial curves you would like to fit??? (Max 4)");

#N_LPC<-floor(as.numeric(readline()));

#if(N_LPC<=4) for(k in 1:N_LPC){LocalCurveFitting()} else print("Wrong: Value should be less than 5");

LocalCurveFitting(1);
LocalCurveFitting(2);
LocalCurveFitting(3);
LocalCurveFitting(4);

legend(0,40,c("LINEAR REGRESSION","KERNEL REGRESSION (NW)","LOCAL POLYNOMIAL DEGREE 1","LOCAL POLYNOMIAL DEGREE 2","LOCAL POLYNOMIAL DEGREE 3","LOCAL POLYNOMIAL DEGREE 4"),col=c(2,1,2,3,4,5),lty=c(2,1,1,1,1,1),cex=0.7);

print("HIT <ENTER> TO CONTINUE ...");
readline();
print("WAIT....PROCESSING");

#=================================
#        SPLINE SMOOTHING
#=================================

Spline_Smoothing=function(Time,Temp,smooth){

Time1<-Time/max(abs(Time));
Temp1<-Temp/max(abs(Temp));

Data_Mat<-as.matrix(cbind(Time1,Temp1));
Sorted_Data_Mat<-Data_Mat[order(Data_Mat[,1]),];

X<-Sorted_Data_Mat[,1];
Y<-Sorted_Data_Mat[,2];

N<-length(X);

A<-vector(mode="numeric",length=N);
B<-vector(mode="numeric",length=N);
C<-vector(mode="numeric",length=N);
E<-vector(mode="numeric",length=N);

dy<-1.06*sd(Y)*N^(-0.2);

Diff_Matrix<-matrix(0,N-1,N);

for(i in 1:(N-1)){

	Diff_Matrix[i,i]=-1;
	Diff_Matrix[i,i+1]=1;

}

H<-Diff_Matrix%*%X;

small_number<-0.0005;

for(i in 1:(N-1)){

if(H[i]==0) H[i]<-small_number;

}

D<-diag(dy,nrow=N,ncol=N);

T<-matrix(data=0,nrow=(N-2),ncol=(N-2));

Q<-matrix(data=0,nrow=N,ncol=(N-2));

for(i in 2:(N-3)){

	T[i,i]<-2*(H[(i-1)]+H[i])/3;
	T[i,(i+1)]<-H[i]/3;
	T[(i+1),i]<-H[i]/3;
}

for(i in 2:(N-2)){

	Q[(i-1),i]<-1/H[(i-1)];
	Q[i,i]<--1/H[(i-1)]-1/H[i];
	Q[(i+1),i]<-1/H[i];
}

#smooth<-0.5;

library(MASS);

C1<-ginv(t(Q)%*%(D^2)%*%Q+smooth*T)%*%t(Q)%*%Y*smooth;

for(i in 1:(N-2)){

C[(i+1)]<-C1[i];
}

A<-Y-(1/smooth)*(D^2)%*%Q%*%C1;

for(i in 1:(N-1)){

E[i]<-(C[(i+1)]-C[i])/(3*H[i]);
B[i]<-(A[(i+1)]-A[i])/H[i]-C[i]*H[i]-D[i]*H[i]*H[i];

}

X1<-vector(mode="numeric",length=N);

X1[1]<-X[1];
X1[N]<-X[N];
for(i in 2:N){

X1[i]<-X1[(i-1)]+H[(i-1)];

}


ITER<-500;

STEP<-(floor(max(X))-ceiling(min(X)))/(ITER-1);

XX<-seq(ceiling(min(X)),floor(max(X)),by=STEP);

XX1<-sort(union(XX,X1));

NN<-length(XX1);

FX<-matrix(data=0,nrow=NN,ncol=7);

for(i in 1:NN){

	for(j in 1:N){

	if(XX1[i]<X1[j]){

		FX[i,1]<-XX1[i];
		FX[i,2]<-X1[(j-1)];
		FX[i,3]<-A[(j-1)];
		FX[i,4]<-B[(j-1)];
		FX[i,5]<-C[(j-1)];
		FX[i,6]<-E[(j-1)];
		break;
			}
		}

	}
FX[NN,1]<-XX1[NN];
FX[NN,2]<-X1[N];
FX[NN,3]<-A[N];
FX[NN,4]<-B[N];
FX[NN,5]<-C[N];
FX[NN,6]<-E[N];

for(i in 1:NN){

FX[i,7]<-FX[i,3]+FX[i,4]*(FX[i,1]-FX[i,2])+FX[i,5]*((FX[i,1]-FX[i,2])^2)+FX[i,6]*((FX[i,1]-FX[i,2])^3);

}

#plot(X*max(abs(Time)),Y*max(abs(Temp)),xlab="TIME",ylab="TEMP",main="SPLINE SMOOTHING");
lines(FX[,1]*max(abs(Time)),FX[,7]*max(abs(Temp)),col=floor(smooth*10));

}

plot(Time,Temp,xlab="TIME",ylab="TEMP",main="SPLINE SMOOTHING AT DIFFERENT LEVELS OF SMOOTHING");

Spline_Smoothing(Time,Temp,0.1);

Spline_Smoothing(Time,Temp,0.5);

Spline_Smoothing(Time,Temp,1);

Spline_Smoothing(Time,Temp,1.5);

Spline_Smoothing(Time,Temp,2);

Spline_Smoothing(Time,Temp,2.5);

Spline_Smoothing(Time,Temp,3);

Spline_Smoothing(Time,Temp,3.5);

Spline_Smoothing(Time,Temp,4);

lines(smooth.spline(Time,Temp,spar=0.8),lty=2);

print("HIT <ENTER> TO CONTINUE ...");
readline();
print("WAIT....PROCESSING");

##################################
#      WAVELET SIMULATION
##################################

Wavelet_Approximation=function(Time,Temp){

#Haar_Wavelet=function(Time,Temp)

#Daubacious_D4=function(Time,Temp)

#AUGMENTATION OF THE DATA SET
#------------------------------------

Max_Time=max(Time);

N<-length(Time);

K<-ceiling(log(N,2));

N1<-2^K;

Add_N<-N1-N;

Data_Mat<-as.matrix(cbind(Time,Temp));
Sorted_Data_Mat<-Data_Mat[order(Data_Mat[,1]),];


Temp_Aug1<-vector(mode="numeric",length=N1);

Time_Aug1<-vector(mode="numeric",length=N1);

for(i in 1:N){

Temp_Aug1[i]<-Sorted_Data_Mat[i,2];
Time_Aug1[i]<-Sorted_Data_Mat[i,1];

}

for(i in 1:Add_N){

Temp_Aug1[(N+i)]<-0;
Time_Aug1[(N+i)]<-Sorted_Data_Mat[N,1]+0.001*i;

}

Temp_Aug<-Temp_Aug1/(max(abs(Temp_Aug1)));
Time_Aug<-Time_Aug1/(max(abs(Time_Aug1)));

#FORWARD TRANSFORM
#-----------------------------------------

#Notes
#There are N1 = 256 (2^8) data points

#All additional redundant data have been added after the max.time and =temp@Max.Time

Haar0<-0.5;
Haar1<-0.5;
Wave_Haar0<-0.5;
Wave_Haar1<--0.5

h0<-(1+sqrt(3))/(4*sqrt(2));
h1<-(3+sqrt(3))/(4*sqrt(2));
h2<-(3-sqrt(3))/(4*sqrt(2));
h3<-(1-sqrt(3))/(4*sqrt(2));
g0<-h3;
g1<--h2;
g2<-h1;
g3<--h0;

hm<-1/4;

HAAR_Forward_0<-matrix(data=0,nrow=N1,ncol=N1);
HAAR_Forward_1<-matrix(data=0,nrow=N1/2,ncol=N1/2);
HAAR_Forward_2<-matrix(data=0,nrow=N1/4,ncol=N1/4);
HAAR_Forward_3<-matrix(data=0,nrow=N1/8,ncol=N1/8);
HAAR_Forward_4<-matrix(data=0,nrow=N1/16,ncol=N1/16);
HAAR_Forward_5<-matrix(data=0,nrow=N1/32,ncol=N1/32);
HAAR_Forward_6<-matrix(data=0,nrow=N1/64,ncol=N1/64);

HAAR_Backward_6<-matrix(data=0,nrow=N1/64,ncol=N1/64);
HAAR_Backward_5<-matrix(data=0,nrow=N1/32,ncol=N1/32);
HAAR_Backward_4<-matrix(data=0,nrow=N1/16,ncol=N1/16);
HAAR_Backward_3<-matrix(data=0,nrow=N1/8,ncol=N1/8);
HAAR_Backward_2<-matrix(data=0,nrow=N1/4,ncol=N1/4);
HAAR_Backward_1<-matrix(data=0,nrow=N1/2,ncol=N1/2);
HAAR_Backward_0<-matrix(data=0,nrow=N1,ncol=N1);

D4_Forward_0<-matrix(data=0,nrow=N1,ncol=N1);
D4_Forward_1<-matrix(data=0,nrow=N1/2,ncol=N1/2);
D4_Forward_2<-matrix(data=0,nrow=N1/4,ncol=N1/4);
D4_Forward_3<-matrix(data=0,nrow=N1/8,ncol=N1/8);
D4_Forward_4<-matrix(data=0,nrow=N1/16,ncol=N1/16);
D4_Forward_5<-matrix(data=0,nrow=N1/32,ncol=N1/32);
D4_Forward_6<-matrix(data=0,nrow=N1/64,ncol=N1/64);

D4_Backward_4<-matrix(data=0,nrow=N1/64,ncol=N1/64);
D4_Backward_4<-matrix(data=0,nrow=N1/32,ncol=N1/32);
D4_Backward_4<-matrix(data=0,nrow=N1/16,ncol=N1/16);
D4_Backward_3<-matrix(data=0,nrow=N1/8,ncol=N1/8);
D4_Backward_2<-matrix(data=0,nrow=N1/4,ncol=N1/4);
D4_Backward_1<-matrix(data=0,nrow=N1/2,ncol=N1/2);
D4_Backward_0<-matrix(data=0,nrow=N1,ncol=N1);

Mean_Filter_Basis_HAAR<-matrix(data=0,nrow=N1,ncol=N1);
Mean_Filter_Basis_D4<-matrix(data=0,nrow=N1,ncol=N1);


#HAAR AND D4 BASIS FUNCTION MATRICES
#-----------------------------------

i<-1;

while(i<=(N1-3)){

HAAR_Forward_0[i,i]<-Haar0;
HAAR_Forward_0[i,(i+1)]<-Haar1;
HAAR_Forward_0[(i+1),i]<-Wave_Haar0;
HAAR_Forward_0[(i+1),(i+1)]<-Wave_Haar1;

HAAR_Backward_0[i,i]<-1;
HAAR_Backward_0[i,(i+1)]<-1;
HAAR_Backward_0[(i+1),i]<-1;
HAAR_Backward_0[(i+1),(i+1)]<--1;

D4_Forward_0[i,i]<-h0;
D4_Forward_0[i,(i+1)]<-h1;
D4_Forward_0[i,(i+2)]<-h2;
D4_Forward_0[i,(i+3)]<-h3;

D4_Forward_0[(i+1),i]<-g0;
D4_Forward_0[(i+1),(i+1)]<-g1;
D4_Forward_0[(i+1),(i+2)]<-g2;
D4_Forward_0[(i+1),(i+3)]<-g3;

D4_Backward_0[i,i]<-h2;
D4_Backward_0[i,(i+1)]<-g2;
D4_Backward_0[i,(i+2)]<-h0;
D4_Backward_0[i,(i+3)]<-g0;

D4_Backward_0[(i+1),i]<-h3;
D4_Backward_0[(i+1),(i+1)]<-g3;
D4_Backward_0[(i+1),(i+2)]<-h1;
D4_Backward_0[(i+1),(i+3)]<-g1;

Mean_Filter_Basis_HAAR[i,i]<-Haar0;
Mean_Filter_Basis_HAAR[i,(i+1)]<-Haar1;
Mean_Filter_Basis_HAAR[(i+1),i]<-Haar0;
Mean_Filter_Basis_HAAR[(i+1),(i+2)]<-Haar1;

Mean_Filter_Basis_D4[i,i]<-hm;
Mean_Filter_Basis_D4[i,(i+1)]<-hm;
Mean_Filter_Basis_D4[i,(i+2)]<-hm;
Mean_Filter_Basis_D4[i,(i+3)]<-hm;

Mean_Filter_Basis_D4[(i+1),(i+1)]<-hm;
Mean_Filter_Basis_D4[(i+1),(i+2)]<-hm;
Mean_Filter_Basis_D4[(i+1),(i+3)]<-hm;
if(i<(N1-3)) Mean_Filter_Basis_D4[(i+1),(i+4)]<-hm;

i=i+2;
}

HAAR_Forward_0[(N1-1),(N1-1)]<-Haar0;
HAAR_Forward_0[(N1-1),((N1-1)+1)]<-Haar1;
HAAR_Forward_0[((N1-1)+1),(N1-1)]<-Wave_Haar0;
HAAR_Forward_0[((N1-1)+1),((N1-1)+1)]<-Wave_Haar1;

Mean_Filter_Basis_HAAR[(N1-1),(N1-1)]<-Haar0;
Mean_Filter_Basis_HAAR[(N1-1),(N1)]<-Haar1;
Mean_Filter_Basis_HAAR[(N1),(N1)]<-Haar1;

Mean_Filter_Basis_HAAR[(N1-1),(N1-1)]<-hm;
Mean_Filter_Basis_HAAR[(N1-1),(N1)]<-hm;
Mean_Filter_Basis_HAAR[(N1),(N1)]<-hm;

HAAR_Backward_0[(N1-1),(N1-1)]<-1;
HAAR_Backward_0[(N1-1),((N1-1)+1)]<-1;
HAAR_Backward_0[((N1-1)+1),(N1-1)]<-1;
HAAR_Backward_0[((N1-1)+1),((N1-1)+1)]<--1;

HAAR_Forward_1<-HAAR_Forward_0[1:(N1/2),1:(N1/2)];
HAAR_Forward_2<-HAAR_Forward_0[1:(N1/4),1:(N1/4)];
HAAR_Forward_3<-HAAR_Forward_0[1:(N1/8),1:(N1/8)];
HAAR_Forward_4<-HAAR_Forward_0[1:(N1/16),1:(N1/16)];
HAAR_Forward_5<-HAAR_Forward_0[1:(N1/32),1:(N1/32)];
HAAR_Forward_6<-HAAR_Forward_0[1:(N1/64),1:(N1/64)];

HAAR_Backward_6<-HAAR_Backward_0[1:(N1/64),1:(N1/64)];
HAAR_Backward_5<-HAAR_Backward_0[1:(N1/32),1:(N1/32)];
HAAR_Backward_4<-HAAR_Backward_0[1:(N1/16),1:(N1/16)];
HAAR_Backward_3<-HAAR_Backward_0[1:(N1/8),1:(N1/8)];
HAAR_Backward_2<-HAAR_Backward_0[1:(N1/4),1:(N1/4)];
HAAR_Backward_1<-HAAR_Backward_0[1:(N1/2),1:(N1/2)];

D4_Forward_1<-D4_Forward_0[1:(N1/2),1:(N1/2)];
D4_Forward_2<-D4_Forward_0[1:(N1/4),1:(N1/4)];
D4_Forward_3<-D4_Forward_0[1:(N1/8),1:(N1/8)];
D4_Forward_4<-D4_Forward_0[1:(N1/16),1:(N1/16)];
D4_Forward_5<-D4_Forward_0[1:(N1/32),1:(N1/32)];
D4_Forward_6<-D4_Forward_0[1:(N1/64),1:(N1/64)];

D4_Backward_6<-D4_Forward_0[1:(N1/64),1:(N1/64)];
D4_Backward_5<-D4_Forward_0[1:(N1/32),1:(N1/32)];
D4_Backward_4<-D4_Forward_0[1:(N1/16),1:(N1/16)];
D4_Backward_3<-D4_Forward_0[1:(N1/8),1:(N1/8)];
D4_Backward_2<-D4_Forward_0[1:(N1/4),1:(N1/4)];
D4_Backward_1<-D4_Forward_0[1:(N1/2),1:(N1/2)];

#PASS FILTER SELECTION MATRICES & STAGE VALUE MATRICES
#-----------------------------------------------------

PASS_1<-matrix(data=0,nrow=N1/2,ncol=N1);
PASS_2<-matrix(data=0,nrow=N1/4,ncol=N1/2);
PASS_3<-matrix(data=0,nrow=N1/8,ncol=N1/4);
PASS_4<-matrix(data=0,nrow=N1/16,ncol=N1/8);
PASS_5<-matrix(data=0,nrow=N1/32,ncol=N1/16);
PASS_6<-matrix(data=0,nrow=N1/64,ncol=N1/32);

HAAR_Temp_Forward_0<-vector(mode="numeric",length=N1);
HAAR_Temp_Forward_1<-vector(mode="numeric",length=N1/2);
HAAR_Temp_Forward_2<-vector(mode="numeric",length=N1/4);
HAAR_Temp_Forward_3<-vector(mode="numeric",length=N1/8);
HAAR_Temp_Forward_4<-vector(mode="numeric",length=N1/16);
HAAR_Temp_Forward_5<-vector(mode="numeric",length=N1/32);
HAAR_Temp_Forward_6<-vector(mode="numeric",length=N1/64);

D4_Temp_Forward_0<-vector(mode="numeric",length=N1);
D4_Temp_Forward_1<-vector(mode="numeric",length=N1/2);
D4_Temp_Forward_2<-vector(mode="numeric",length=N1/4);
D4_Temp_Forward_3<-vector(mode="numeric",length=N1/8);
D4_Temp_Forward_4<-vector(mode="numeric",length=N1/16);
D4_Temp_Forward_5<-vector(mode="numeric",length=N1/32);
D4_Temp_Forward_6<-vector(mode="numeric",length=N1/64);

HAAR_Temp_Backward_0<-vector(mode="numeric",length=N1);
HAAR_Temp_Backward_1<-vector(mode="numeric",length=N1/2);
HAAR_Temp_Backward_2<-vector(mode="numeric",length=N1/4);
HAAR_Temp_Backward_3<-vector(mode="numeric",length=N1/8);
HAAR_Temp_Backward_4<-vector(mode="numeric",length=N1/16);
HAAR_Temp_Backward_5<-vector(mode="numeric",length=N1/32);
HAAR_Temp_Backward_6<-vector(mode="numeric",length=N1/64);

D4_Temp_Backward_0<-vector(mode="numeric",length=N1);
D4_Temp_Backward_1<-vector(mode="numeric",length=N1/2);
D4_Temp_Backward_2<-vector(mode="numeric",length=N1/4);
D4_Temp_Backward_3<-vector(mode="numeric",length=N1/8);
D4_Temp_Backward_4<-vector(mode="numeric",length=N1/16);
D4_Temp_Backward_5<-vector(mode="numeric",length=N1/32);
D4_Temp_Backward_6<-vector(mode="numeric",length=N1/64);

i=1;

while(i<=(N1/2)){

PASS_1[i,(2*i-1)]<-1;

i=i+1;
}

PASS_2<-PASS_1[1:(N1/4),1:(N1/2)];
PASS_3<-PASS_1[1:(N1/8),1:(N1/4)];
PASS_4<-PASS_1[1:(N1/16),1:(N1/8)];
PASS_5<-PASS_1[1:(N1/32),1:(N1/16)];
PASS_6<-PASS_1[1:(N1/64),1:(N1/32)];

#--------------------------------
#HAAR SIMULATION 
#--------------------------------
#FORWARD FILTERING PASS
#--------------------------------

HAAR_Temp_Forward_0<-HAAR_Forward_0%*%Temp_Aug;
HAAR_Temp_Forward_1<-HAAR_Forward_1%*%PASS_1%*%HAAR_Temp_Forward_0;
HAAR_Temp_Forward_2<-HAAR_Forward_2%*%PASS_2%*%HAAR_Temp_Forward_1;
HAAR_Temp_Forward_3<-HAAR_Forward_3%*%PASS_3%*%HAAR_Temp_Forward_2;
HAAR_Temp_Forward_4<-HAAR_Forward_4%*%PASS_4%*%HAAR_Temp_Forward_3;
HAAR_Temp_Forward_5<-HAAR_Forward_5%*%PASS_5%*%HAAR_Temp_Forward_4;
HAAR_Temp_Forward_6<-HAAR_Forward_6%*%PASS_6%*%HAAR_Temp_Forward_5;

#-------------------------------
#BACKWARD RECONSTRUCTION PASS
#--------------------------------

HAAR_Temp_Backward_6<-HAAR_Temp_Forward_6;

Temporary_Reconstruction<-HAAR_Backward_6%*%HAAR_Temp_Backward_6;

i=1;
while(i<=((N1/32)-1)){

HAAR_Temp_Backward_5[i]<-Temporary_Reconstruction[i];
HAAR_Temp_Backward_5[i]<-HAAR_Temp_Forward_5[(i+1)];

i=i+2;
}

Temporary_Reconstruction<-HAAR_Backward_5%*%HAAR_Temp_Backward_5;

i=1;
while(i<=((N1/16)-1)){

HAAR_Temp_Backward_4[i]<-Temporary_Reconstruction[i];
HAAR_Temp_Backward_4[i]<-HAAR_Temp_Forward_4[(i+1)];

i=i+2;
}

Temporary_Reconstruction<-HAAR_Backward_4%*%HAAR_Temp_Backward_4;

i=1;
while(i<=((N1/8)-1)){

HAAR_Temp_Backward_3[i]<-Temporary_Reconstruction[i];
HAAR_Temp_Backward_3[i]<-HAAR_Temp_Forward_3[(i+1)];

i=i+2;
}

Temporary_Reconstruction<-HAAR_Backward_3%*%HAAR_Temp_Backward_3;

i=1;
while(i<=((N1/4)-1)){

HAAR_Temp_Backward_2[i]<-Temporary_Reconstruction[i];
HAAR_Temp_Backward_2[i]<-HAAR_Temp_Forward_2[(i+1)];

i=i+2;
}

Temporary_Reconstruction<-HAAR_Backward_2%*%HAAR_Temp_Backward_2;

i=1;
while(i<=((N1/2)-1)){

HAAR_Temp_Backward_1[i]<-Temporary_Reconstruction[i];
HAAR_Temp_Backward_1[i]<-HAAR_Temp_Forward_1[(i+1)];

i=i+2;
}

Temporary_Reconstruction<-HAAR_Backward_1%*%HAAR_Temp_Backward_1;

i=1;
while(i<=((N1)-1)){

HAAR_Temp_Backward_0[i]<-Temporary_Reconstruction[i];
HAAR_Temp_Backward_0[i]<-HAAR_Temp_Forward_0[(i+1)];

i=i+2;
}

Wavelet_Filter_Haar<-(Mean_Filter_Basis_HAAR%*%Temp_Aug)*max(abs(Temp_Aug1));

plot(Time_Aug*max(abs(Time_Aug1)),Temp_Aug*max(abs(Temp_Aug1)),xlab="Time",ylab="Temp",main="WAVELET APPROXIMATIONS");
lines(Time_Aug*max(abs(Time_Aug1)),Wavelet_Filter_Haar+HAAR_Temp_Backward_0*max(abs(Temp_Aug1)),col=2);

print("HIT <ENTER> TO CONTINUE ...");
readline();
print("WAIT....PROCESSING");

#--------------------------------
#DAUBECHIES 4 SIMULATION 
#--------------------------------
#FORWARD FILTERING PASS
#--------------------------------

D4_Temp_Forward_0<-D4_Forward_0%*%Temp_Aug;
D4_Temp_Forward_1<-D4_Forward_1%*%PASS_1%*%D4_Temp_Forward_0;
D4_Temp_Forward_2<-D4_Forward_2%*%PASS_2%*%D4_Temp_Forward_1;
D4_Temp_Forward_3<-D4_Forward_3%*%PASS_3%*%D4_Temp_Forward_2;
D4_Temp_Forward_4<-D4_Forward_4%*%PASS_4%*%D4_Temp_Forward_3;
D4_Temp_Forward_5<-D4_Forward_5%*%PASS_5%*%D4_Temp_Forward_4;
D4_Temp_Forward_6<-D4_Forward_6%*%PASS_6%*%D4_Temp_Forward_5;

#-------------------------------
#BACKWARD RECONSTRUCTION PASS
#--------------------------------

D4_Temp_Backward_6<-D4_Temp_Forward_6;

Temporary_Reconstruction<-D4_Backward_6%*%D4_Temp_Backward_6;

i=1;
while(i<=((N1/32)-1)){

D4_Temp_Backward_5[i]<-Temporary_Reconstruction[i];
D4_Temp_Backward_5[i]<-D4_Temp_Forward_5[(i+1)];

i=i+2;
}

Temporary_Reconstruction<-D4_Backward_5%*%D4_Temp_Backward_5;

i=1;
while(i<=((N1/16)-1)){

D4_Temp_Backward_4[i]<-Temporary_Reconstruction[i];
D4_Temp_Backward_4[i]<-D4_Temp_Forward_4[(i+1)];

i=i+2;
}

Temporary_Reconstruction<-D4_Backward_4%*%D4_Temp_Backward_4;

i=1;
while(i<=((N1/8)-1)){

D4_Temp_Backward_3[i]<-Temporary_Reconstruction[i];
D4_Temp_Backward_3[i]<-D4_Temp_Forward_3[(i+1)];

i=i+2;
}

Temporary_Reconstruction<-D4_Backward_3%*%D4_Temp_Backward_3;

i=1;
while(i<=((N1/4)-1)){

D4_Temp_Backward_2[i]<-Temporary_Reconstruction[i];
D4_Temp_Backward_2[i]<-D4_Temp_Forward_2[(i+1)];

i=i+2;
}

Temporary_Reconstruction<-D4_Backward_2%*%D4_Temp_Backward_2;

i=1;
while(i<=((N1/2)-1)){

D4_Temp_Backward_1[i]<-Temporary_Reconstruction[i];
D4_Temp_Backward_1[i]<-D4_Temp_Forward_1[(i+1)];

i=i+2;
}

Temporary_Reconstruction<-D4_Backward_1%*%D4_Temp_Backward_1;

i=1;
while(i<=((N1)-1)){

D4_Temp_Backward_0[i]<-Temporary_Reconstruction[i];
D4_Temp_Backward_0[i]<-D4_Temp_Forward_0[(i+1)];

i=i+2;
}

Wavelet_Filter_D4<-(Mean_Filter_Basis_D4%*%Temp_Aug)*max(abs(Temp_Aug1));

lines(Time_Aug*max(abs(Time_Aug1)),Wavelet_Filter_Haar+D4_Temp_Backward_0*max(abs(Temp_Aug1)),col=4,lty=2);

legend(0,40,c("HAAR WAVELETS","DAUBECHIES 4 WAVELET","WAVE-THRESH"),col=c(2,4,3),lty=c(1,2,3),cex=0.7);

library(wavethresh);
Temp.Decompose<-wd(Temp_Aug1);
Temp.Reconstruct<-wr(threshold(Temp.Decompose));

lines(Time_Aug1,Temp.Reconstruct,col=3,lty=3);

}

Wavelet_Approximation(Time,Temp);

print("HIT <ENTER> TO CONTINUE ...");
readline();

################################################
##########     BOOTSTRAPPING     ###############
################################################

################################################
# BOOTSTRAP CONFIDENCE INTERVAL FOR SLOPE COEFF 
################################################

#NOTE: USED A SIMPLE LINEAR MODEL WITHOUT ANY TRANSFORMATION OF THE DATA
#-----------------------------------------------------------------------

# Data Set used, X=TIME, Y=TEMP
#--------------------------------
#Method Used: Paired Boot-strap, Residual Boot Strap and External Boot-Strap
#---------------------------------------------------------------------------

Data.Original<-as.matrix(cbind(Temp,Time));

#PAIRED BOOTSTRAP FUNCTION
#-------------------------

PairedBootstrap<-function(Data){

N=nrow(Data);

Index<-seq(1,N);

Index_Boot<-sample(Index,size=N,replace=TRUE);

Data_Boot<-Data[Index_Boot,];

return(Data_Boot);

}

#RESIDUAL BOOTSTRAP FUNCTION
#----------------------------

ResidualBootstrap<-function(Data){

N=nrow(Data);

Linear_Regression<-lm(Data[,1]~Data[,2]);

Response_Hat<-predict(Linear_Regression)

Residuals<-Data[,1]-Response_Hat;

Boot_Residuals<-sample(Residuals,size=N,replace=TRUE);

Response_Boot<-Response_Hat+Boot_Residuals;

Data_Boot<-as.matrix(cbind(Response_Boot,Data[,2]));

return(Data_Boot);

}

#EXTERNAL BOOTSTRAP
#------------------

ExternalBootstrap<-function(Data){

N=nrow(Data);

Linear_Regression<-lm(Data[,1]~Data[,2]);

Response_Hat<-predict(Linear_Regression)

Error_SD<-sd(residuals(Linear_Regression));

External_Error<-rnorm(N,0,Error_SD);

Response_Boot<-Response_Hat+External_Error;

Data_Boot<-as.matrix(cbind(Response_Boot,Data[,2]));

return(Data_Boot);

}

B<-1000;

R1<-lm(Data.Original[,1]~Data.Original[,2]);

Beta_Coeff<-coefficients(R1)[[2]];

Beta_X<-sort(runif(B,(Beta_Coeff-4*coef(summary(R1))[2,2]),(Beta_Coeff+4*coef(summary(R1))[2,2])));

plot(Beta_X,dnorm(Beta_X,Beta_Coeff,coef(summary(R1))[2,2]),xlab="Slope Values",ylab="Probability Distribution",main="PDF Plot of Slope Values of Regression",type="l",ylim=c(0,0.7));

R1_Conf_Int<-confint(R1)[2,];

Boot_Strap_Slope<-matrix(data=0,nrow=B,ncol=3);

PN<-500;

Predicted_Paired<-matrix(data=0,nrow=B,ncol=PN);

Predicted_Residual<-matrix(data=0,nrow=B,ncol=PN);

Predicted_External<-matrix(data=0,nrow=B,ncol=PN);

Prediction_Points<-seq(floor(min(Time)),ceiling(max(Time)),((-floor(min(Time))+ceiling(max(Time)))/(PN-1)));

for(i in 1:B){

D1<-PairedBootstrap(Data.Original);

D2<-ResidualBootstrap(Data.Original);

D3<-ExternalBootstrap(Data.Original);

R1<-lm(D1[,1]~D1[,2]);

R2<-lm(D2[,1]~D2[,2]);

R3<-lm(D3[,1]~D3[,2]);

Boot_Strap_Slope[i,1]<-coefficients(R1)[[2]];

Boot_Strap_Slope[i,2]<-coefficients(R2)[[2]];

Boot_Strap_Slope[i,3]<-coefficients(R3)[[2]];

for(j in 1:PN){

Predicted_Paired[i,j]<-coefficients(R1)[[1]]+coefficients(R1)[[2]]*Prediction_Points[j];

Predicted_Residual[i,j]<-coefficients(R2)[[1]]+coefficients(R2)[[2]]*Prediction_Points[j];

Predicted_External[i,j]<-coefficients(R3)[[1]]+coefficients(R3)[[2]]*Prediction_Points[j];

}

}

lines(density(Boot_Strap_Slope[,1]),col=2);
lines(density(Boot_Strap_Slope[,2]),col=3);
lines(density(Boot_Strap_Slope[,3]),col=4);

legend(9.5,0.7,c("Gaussian Model","Paired Bootstrap","Residual Bootstrap","External Bootstrap"),col=c(1,2,3,4),lty=c(1,1,1,1));

print("HIT <ENTER> TO CONTINUE ...");
readline();

Alpha=0.05;

LL_D1<-sort(Boot_Strap_Slope[,1])[floor(0.5*B*Alpha)];
UL_D1<-sort(Boot_Strap_Slope[,1])[floor(B*(1-0.5*Alpha))];

LL_D2<-sort(Boot_Strap_Slope[,2])[floor(0.5*B*Alpha)];
UL_D2<-sort(Boot_Strap_Slope[,2])[floor(B*(1-0.5*Alpha))];

LL_D3<-sort(Boot_Strap_Slope[,3])[floor(0.5*B*Alpha)];
UL_D3<-sort(Boot_Strap_Slope[,3])[floor(B*(1-0.5*Alpha))];

Header_Conf_Table<-c("Slope","2.5%","97.5%","Interval");
LINEAR_GAUSSIAN_MODEL<-c(format(Beta_Coeff,digits=5),format(R1_Conf_Int[[1]],digits=5),format(R1_Conf_Int[[2]],digits=5),format((abs(R1_Conf_Int[[1]]-R1_Conf_Int[[2]])),digits=5));
PAIRED_BOOTSTRAP<-c(format(mean(Boot_Strap_Slope[,1]),digits=5),format(LL_D1,digits=5),format(UL_D1,digits=5),format((abs(LL_D1-UL_D1)),digits=5));
RESIDUAL_BOOTSTRAP<-c(format(mean(Boot_Strap_Slope[,2]),digits=5),format(LL_D2,digits=5),format(UL_D3,digits=5),format((abs(LL_D2-UL_D2)),digits=5));
EXTERNAL_BOOTSTRAP<-c(format(mean(Boot_Strap_Slope[,3]),digits=5),format(LL_D3,digits=5),format(UL_D3,digits=5),format((abs(LL_D3-UL_D3)),digits=5));

Confidence_Table<-data.frame(rbind(Header_Conf_Table,LINEAR_GAUSSIAN_MODEL,PAIRED_BOOTSTRAP,RESIDUAL_BOOTSTRAP,EXTERNAL_BOOTSTRAP));

print("The Confidence Interval Table for The Linear Regression Model:");
print("---------------------------------------------------------------");
print(Confidence_Table);
print("---------------------------------------------------------------");

print("HIT <ENTER> TO CONTINUE ...");
readline();

plot(Time,Temp,xlab="Time",ylab="Temp",main="Confidence Band Plot for Linear Models");
abline(lm(Temp~Time));

Confidence_Int_Boot<-matrix(data=0,nrow=6,ncol=PN);

for(i in 1:PN){

Order_Temp_1<-sort(Predicted_Paired[,i]);

Order_Temp_2<-sort(Predicted_Residual[,i]);

Order_Temp_3<-sort(Predicted_External[,i]);


Confidence_Int_Boot[1,i]<-Order_Temp_1[floor(0.5*B*Alpha)];

Confidence_Int_Boot[2,i]<-Order_Temp_1[ceiling(B*(1-0.5*Alpha))];

Confidence_Int_Boot[3,i]<-Order_Temp_2[floor(0.5*B*Alpha)];

Confidence_Int_Boot[4,i]<-Order_Temp_2[ceiling(B*(1-0.5*Alpha))];

Confidence_Int_Boot[5,i]<-Order_Temp_3[floor(0.5*B*Alpha)];

Confidence_Int_Boot[6,i]<-Order_Temp_3[ceiling(B*(1-0.5*Alpha))];

}

lines(Prediction_Points,Confidence_Int_Boot[1,],col=2)
lines(Prediction_Points,Confidence_Int_Boot[2,],col=2)
lines(Prediction_Points,Confidence_Int_Boot[3,],col=3)
lines(Prediction_Points,Confidence_Int_Boot[4,],col=3)
lines(Prediction_Points,Confidence_Int_Boot[5,],col=4)
lines(Prediction_Points,Confidence_Int_Boot[6,],col=4)

legend(0,40,c("Paired Bootstrap","Residual Bootstrap","External Bootstrap"),col=c(2,3,4),lty=c(1,1,1),cex=0.7);



############################################
##   SPLINE BOOTSTRAP
############################################

N<-nrow(Data.Original);

Spline_Boot_Strap<-matrix(data=0,nrow=B,ncol=N);

Ordered_Data<-Data.Original[order(Data.Original[,2]),];

for(i in 1:B){

S1<-PairedBootstrap(Ordered_Data);

P1<-smooth.spline(S1[,2],S1[,1],spar=0.8);

for(j in 1:N){

Response_Hat<-predict(P1,Ordered_Data[[j,2]]);

Spline_Boot_Strap[i,j]<-Response_Hat[[2]];

}

}

Parameter_Spline_Boot<-matrix(data=0,nrow=4,ncol=N)

for(i in 1:N){

Order_Temp<-sort(Spline_Boot_Strap[,i]);

Parameter_Spline_Boot[1,i]<-mean(Order_Temp);

Parameter_Spline_Boot[2,i]<-median(Order_Temp);

Parameter_Spline_Boot[3,i]<-Order_Temp[floor(0.5*B*Alpha)];

Parameter_Spline_Boot[4,i]<-Order_Temp[floor(B*(1-0.5*Alpha))];

}

plot(Time,Temp,xlab="Time",ylab="Temp",main="Confidence Band Plot for Splines");
lines(Ordered_Data[,2],Parameter_Spline_Boot[1,],type="l",col=2);
lines(Ordered_Data[,2],Parameter_Spline_Boot[3,],type="l",col=3,lwd=3);
lines(Ordered_Data[,2],Parameter_Spline_Boot[4,],type="l",col=4,lwd=3);
lines(Ordered_Data[,2],Parameter_Spline_Boot[2,],type="l",col=5);

legend(0,40,c("MEAN SPLINE","LOWER CONF LIM","UPPER CONF LIM","MEDIAN SPLINE"),col=c(2,3,4,5),lty=c(1,1,1,1),cex=0.7);


print("HIT <ENTER> TO CONTINUE ...");
readline();

###########################################
#     LOCAL POLYNOMIAL BOOTSTRAP
###########################################

N<-nrow(Data.Original);

library(locfit);

Locpol_Boot_Strap<-matrix(data=0,nrow=B,ncol=N);

Ordered_Data<-Data.Original[order(Data.Original[,2]),];

for(i in 1:B){

S1<-PairedBootstrap(Ordered_Data);

P1<-locfit(S1[,1]~S1[,2]);

for(j in 1:N){

Response_Hat<-predict(P1,Ordered_Data[[j,2]]);

Locpol_Boot_Strap[i,j]<-Response_Hat[[1]];

}

}

Parameter_Locpol_Boot<-matrix(data=0,nrow=4,ncol=N)

for(i in 1:N){

Order_Temp<-sort(Locpol_Boot_Strap[,i]);

Parameter_Locpol_Boot[1,i]<-mean(Order_Temp);

Parameter_Locpol_Boot[2,i]<-median(Order_Temp);

Parameter_Locpol_Boot[3,i]<-Order_Temp[floor(0.5*B*Alpha)];

Parameter_Locpol_Boot[4,i]<-Order_Temp[floor(B*(1-0.5*Alpha))];

}

plot(Time,Temp,xlab="Time",ylab="Temp",main="Confidence Band Plot for Local Polynomial");
lines(Ordered_Data[,2],Parameter_Locpol_Boot[1,],type="l",col=2);
lines(Ordered_Data[,2],Parameter_Locpol_Boot[3,],type="l",col=3,lwd=3);
lines(Ordered_Data[,2],Parameter_Locpol_Boot[4,],type="l",col=4,lwd=3);
lines(Ordered_Data[,2],Parameter_Locpol_Boot[2,],type="l",col=5);

legend(0,40,c("MEAN POLYNOMIAL","LOWER CONF LIM","UPPER CONF LIM","MEDIAN POLYNOMIAL"),col=c(2,3,4,5),lty=c(1,1,1,1),cex=0.7);


print("HIT <ENTER> TO CONTINUE ...");
readline();

###########################################
#     WAVELET BOOTSTRAP
###########################################

#AUGMENTATION OF THE DATA SET
#------------------------------------

Max_Time=max(Time);

N<-length(Time);

K<-ceiling(log(N,2));

N1<-2^K;

Add_N<-N1-N;

Data_Mat<-as.matrix(cbind(Time,Temp));
Sorted_Data_Mat<-Data_Mat[order(Data_Mat[,1]),];


Temp_Aug1<-vector(mode="numeric",length=N1);

Time_Aug1<-vector(mode="numeric",length=N1);

for(i in 1:N){

Temp_Aug1[i]<-Sorted_Data_Mat[i,2];
Time_Aug1[i]<-Sorted_Data_Mat[i,1];

}

for(i in 1:Add_N){

Temp_Aug1[(N+i)]<-0;
Time_Aug1[(N+i)]<-Sorted_Data_Mat[N,1]+0.001*i;

}

#Bootstrap Sequence
#------------------

Data.Original1<-as.matrix(cbind(Temp_Aug1,Time_Aug1));

N1<-nrow(Data.Original1);

library(wavethresh);

Wave_Boot_Strap<-matrix(data=0,nrow=B,ncol=N1);

Ordered_Data<-Data.Original1;

for(i in 1:B){

S1<-PairedBootstrap(Ordered_Data);

S2<-S1[order(S1[,2]),];

Wave_Forward<-wd(S2[,1]);

Wave_Reconstruct<-wr(threshold(Wave_Forward));

Wave_Boot_Strap[i,]<-Wave_Reconstruct;


}

Parameter_Wave_Boot<-matrix(data=0,nrow=4,ncol=N1)

for(i in 1:N1){

Order_Temp<-sort(Wave_Boot_Strap[,i]);

Parameter_Wave_Boot[1,i]<-mean(Order_Temp);

Parameter_Wave_Boot[2,i]<-median(Order_Temp);

Parameter_Wave_Boot[3,i]<-Order_Temp[floor(0.5*B*Alpha)];

Parameter_Wave_Boot[4,i]<-Order_Temp[floor(B*(1-0.5*Alpha))];

}

plot(Time_Aug1,Temp_Aug1,xlab="Time",ylab="Temp",main="Confidence Band Plot for Wavelet");
lines(Ordered_Data[,2],Parameter_Wave_Boot[1,],type="l",col=2);
lines(Ordered_Data[,2],Parameter_Wave_Boot[3,],type="l",col=3,lwd=5);
lines(Ordered_Data[,2],Parameter_Wave_Boot[4,],type="l",col=4,lwd=5);
lines(Ordered_Data[,2],Parameter_Wave_Boot[2,],type="l",col=5);

legend(0,40,c("MEAN WAVELET","LOWER CONF LIM","UPPER CONF LIM","MEDIAN WAVELET"),col=c(2,3,4,5),lty=c(1,1,1,1),cex=0.7);


print("HIT <ENTER> TO CONTINUE ...");
readline();



