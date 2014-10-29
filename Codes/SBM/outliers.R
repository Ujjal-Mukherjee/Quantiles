
library(FastHCS)
library(cepp)

## colon data
## FastHCS implementation
data(Colo)
system.time(Fit <- FastHCS(x=Colon$X, k=20))

colvec<-rep("orange",length(Colon$Y))
colvec[Colon$Y==1]<-"blue"
SDIND<-Fit$rew.fit$sd/Fit$rew.fit$cutoff.sd
ODIND<-Fit$rew.fit$od/Fit$rew.fit$cutoff.od
plot(SDIND,ODIND,col=colvec,pch=16)
abline(h=1,col="red",lty=2)
abline(v=1,col="red",lty=2)
