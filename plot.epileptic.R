library(MASS)
rm(list=ls())
library(survival)
library(joineR)
library(nlme)
library(mvtnorm)
library(mnormt)
library(JM)
library(splines)
library(grid)
library(parallel)
library(pracma)
library(ggplot2)
library(tidyverse)

#epileptic <- read.table("epileptic.txt", header = TRUE)
head(epileptic)

epileptic$time <- epileptic$time / 365.25
epileptic$with.time <- epileptic$with.time / 365.25
epileptic$status <- rep("Censored", nrow(epileptic))
epileptic$status[epileptic$with.status.uae == 1] <- "UAE"
epileptic$status[epileptic$with.status.isc == 1] <- "ISC"
epileptic$treat<-as.numeric(epileptic$treat=="LTG")

epileptic$cause[epileptic$status == "Censored"]<-0
epileptic$cause[epileptic$status == "ISC"] <- 1
epileptic$cause[epileptic$status =="UAE"] <- 2

head(epileptic,70)


################
#Spaghetti plots of calib dose versus revers time 
#stratified by event time mechanism and drug type

par(op)
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,3))
attach(epileptic)

plot( time[id==1]-with.time[id==1],dose[id==1],type="l",col="red",ylim=c(0,7),xlim=c(-7,0) ,main="(Censored)",
     ,xlab="",ylab="")
index=1:length(id)
index0=index[cause==0]
index1=index[cause==1]
index2=index[cause==2]
i0=as.numeric(attributes(table(id[index0]))$dimnames[[1]])
i1=as.numeric(attributes(table(id[index1]))$dimnames[[1]])
i2=as.numeric(attributes(table(id[index2]))$dimnames[[1]])


# treat: 0=CBZ
for(i in i0){ 
  if(treat[i]==0){ 
lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
    }

plot( time[id==5]-with.time[id==5],dose[id==5],type="l",col="red",ylim=c(0,7), xlim=c(-7,0),main="(ISC)",
      xlab="",ylab="")

for(i in i1){ 
  if(treat[i]==0){ 
    lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
  
}

plot( time[id==7]-with.time[id==7],dose[id==7],type="l",col="red",ylim=c(0,7), xlim=c(-7,0) ,main="(UAE)",
      xlab="",ylab="")
for(i in i2){ 
  if(treat[i]==0){ 
    lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
}

mtext(text="Time before treatment failure or censoring (years)",side=1,line=0,outer=TRUE)
mtext(text="Calibrated dose",side=2,line=0,outer=TRUE)

mtext("CBZ", side=4, line=0.7)

##############

plot( time[id==12]-with.time[id==12],dose[id==12],type="l",ylim=c(0,7),xlim=c(-7,0) ,
      xlab="",ylab="",col="red")
# treat: 0=CBZ
for(i in i0){ 
 
  if(treat[i]==1){ 
    lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
  }


plot( time[id==42]-with.time[id==42],dose[id==42],type="l",ylim=c(0,7), xlim=c(-7,0) ,
      xlab="",ylab="",col="red")

for(i in i1){ 
 
  if(treat[i]==1){ 
    lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
}

plot( time[id==52]-with.time[id==52],dose[id==52],type="l",ylim=c(0,7), xlim=c(-7,0) ,
      xlab="",ylab="",col="red")
for(i in i2){ 
 
  if(treat[i]==1){ 
    lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
}

mtext("LTG", side=4, line=0.7)



#######################     RANDOMLY SELECTED
#Spaghetti plots of calib dose versus revers time 
#stratified by event time mechanism and drug type

par(op)
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(2,3))


plot( time[id==1]-with.time[id==1],dose[id==1],type="l",col="red",ylim=c(0,7),xlim=c(-7,0) ,main="(Censored)",
     ,xlab="",ylab="")


# treat: 0=CBZ
for(i in sample(i0,25)){ 
  if(treat[i]==0){ 
lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
    }

plot( time[id==5]-with.time[id==5],dose[id==5],type="l",col="red",ylim=c(0,7), xlim=c(-7,0),main="(ISC)",
      xlab="",ylab="")

for(i in sample(i1,25)){ 
  if(treat[i]==0){ 
    lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
  
}

plot( time[id==7]-with.time[id==7],dose[id==7],type="l",col="red",ylim=c(0,7), xlim=c(-7,0) ,main="(UAE)",
      xlab="",ylab="")
for(i in sample(i2,25)){ 
  if(treat[i]==0){ 
    lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
}

mtext(text="Time before treatment failure or censoring (years)",side=1,line=0,outer=TRUE)
mtext(text="Calibrated dose",side=2,line=0,outer=TRUE)

mtext("CBZ", side=4, line=0.7)

##############

plot( time[id==12]-with.time[id==12],dose[id==12],type="l",ylim=c(0,7),xlim=c(-7,0) ,
      xlab="",ylab="",col="red")
# treat: 0=CBZ
for(i in sample(i0,25)){ 
 
  if(treat[i]==1){ 
    lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
  }


plot( time[id==42]-with.time[id==42],dose[id==42],type="l",ylim=c(0,7), xlim=c(-7,0) ,
      xlab="",ylab="",col="red")

for(i in sample(i1,25)){ 
 
  if(treat[i]==1){ 
    lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
}

plot( time[id==52]-with.time[id==52],dose[id==52],type="l",ylim=c(0,7), xlim=c(-7,0) ,
      xlab="",ylab="",col="red")
for(i in sample(i2,25)){ 
 
  if(treat[i]==1){ 
    lines(time[id==i]-with.time[id==i],dose[id==i],col=i,type="l")
  }
}

mtext("LTG", side=4, line=0.7)




##########Spaghetti plot calib dose versus times separated by drug type

index3=index[treat==0]
index4=index[treat==1]
i3=as.numeric(attributes(table(id[index3]))$dimnames[[1]])
i4=as.numeric(attributes(table(id[index4]))$dimnames[[1]])

par(mfrow=c(1,2))
par(op)
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(1,2))
plot( time[id==1],dose[id==1],type="l",col="red",ylim=c(0,7),xlim=c(0,10) ,main="(CBZ)",
     ,xlab="",ylab="")

for(i in i3){ 
   
lines(time[id==i],dose[id==i],col=i,type="l")
 
    }


plot( time[id==2],dose[id==2],type="l",col="red",ylim=c(0,7),xlim=c(0,10) ,main="(LTG)",
     ,xlab="",ylab="")

for(i in i4){ 
   
lines(time[id==i],dose[id==i],col=i,type="l")
 
    }


mtext(text="Time(years)",side=1,line=0,outer=TRUE)
mtext(text="Calibrated dose",side=2,line=0,outer=TRUE)


########################Randomly selected 
par(op)
par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(1,2))
plot( time[id==1],dose[id==1],type="l",col="red",ylim=c(0,7),xlim=c(0,10) ,main="(CBZ)",
     ,xlab="",ylab="")

for(i in sample(i3,20)){ 
   
lines(time[id==i],dose[id==i],col=i,type="l")
 
    }


plot( time[id==2],dose[id==2],type="l",col="red",ylim=c(0,7),xlim=c(0,10) ,main="(LTG)",
     ,xlab="",ylab="")

for(i in sample(i4,20)){ 
   
lines(time[id==i],dose[id==i],col=i,type="l")
 
    }


mtext(text="Time(years)",side=1,line=0,outer=TRUE)
mtext(text="Calibrated dose",side=2,line=0,outer=TRUE)



########################Kaplan-meier
par(mfrow=c(1,1))

cause1=cause
cause1[cause1==2]=0

library(survival)
ss=Surv(with.time,event=cause1 , type="right")

plot(survfit(ss~1),ylim=c(0.6,1),lwd=1.3,xlab="Time (year)",ylab="Survival Probability",
main="Kaplan-Meier Estimate") 

cause2=cause
cause2[cause2==1]=0
cause2[cause2==2]=1
ss=Surv(with.time,event=cause2 , type="right")

lines(survfit(ss~1),col=2,lwd=1.3) 


legend("bottomleft", c("UAE", "ISC"), lty =1,lwd=1.3,
       col=c("red","black")) 