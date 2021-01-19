library(MASS)
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


data(epileptic)
epileptic
write.table(epileptic, file="epileptic",row.names=TRUE,col.names=TRUE)
epileptic$time<-epileptic$time/365.25
epileptic$with.time<-epileptic$with.time/365.25
epileptic$treat<-as.numeric(epileptic$treat=="LTG")
Data=epileptic
Data
#rm(time)
attach(Data)

#################proposed two-satge approach for JMLC
################################### Longitudinal sub-model

#ctrl<-lmeControl(opt='optim')
#lmeFit<-lme(dose~treat * time,random=~time|id,control=ctrl,data=epileptic)
lmeFit <- lme(dose ~ treat * time,
              random = ~ time | id,
              data = epileptic)
summary(lmeFit)
b<-ranef(lmeFit)
VarCorr(lmeFit)
sigma11<-as.numeric(VarCorr(lmeFit)[1])
sigma22<-as.numeric(VarCorr(lmeFit)[2])
sigma12<-sqrt(sigma11)*sqrt(sigma22)*as.numeric(VarCorr(lmeFit)[8])
sigma<-as.numeric(VarCorr(lmeFit)[6])

# number of measurement occasions for the whole sample
N<-dim(Data)[1]
N
# number of measurement accasions for each subject
ni<-as.vector(tapply(Data$id,Data$id,length))
ni
summary(ni)
data2<- Data[!duplicated(Data$id),]
# number of subjects in the sample
n<-length(data2$id)
n
#head(data2,20)
##############
b<-ranef(lmeFit)
loglik.yb<-logLik(lmeFit)
#data2

delta1<-c()
for(i in 1:n){
if(data2$with.status.isc[i]==1){delta1[i]=1} else {delta1[i]=0}
}

delta2<-c()
for(i in 1:n){
  if(data2$with.status.uae[i]==1){delta2[i]=1} else {delta2[i]=0}
}

########################### Competing risks process
###############
#################### Function of update lambda
twostage.JMLC.lambda<-function(beta0,beta1,beta2,beta3,lambda,gamma,alpha,b,loglik.yb,delta,iters,Data,data2){
N<-dim(Data)[1]
N
# number of measurement accasions for each subject
ni<-as.vector(tapply(Data$id,Data$id,length))
ni
summary(ni)
data2<- Data[!duplicated(Data$id),]
# number of subjects in the sample
n<-length(data2$id)
n

loglik.yb<-logLik(lmeFit)

lgLik<-rep(numeric(),iters)

X1<-rep(1,n)
X2<-data2$with.time

Z.time<-matrix(c(X1,X2),ncol=2)

Ztime.b<- rep(0,n)
for (i in 1:n) {Ztime.b[i]<- Z.time[i,] %*% t(b[i,])}

cur.lgLik<- -10000
EM.lambda<-lambda
EM.gamma<-gamma
EM.alpha<-alpha

for (t in 1:iters)  {
  
  log.hazard<-delta*log(lambda) + delta*( gamma*data2$treat + alpha*(beta0 + beta1*data2$with.time+beta2*data2$treat
  +beta3*(data2$treat*data2$with.time)+Ztime.b))
# print(log.hazard)
  log.survival<-rep(0,n)
  
for (i in 1:n) {
  a1<-(lambda/(alpha*(beta1+beta3*data2$treat[i]+b[i,2])))*exp(gamma*data2$treat[i]+alpha*(beta0+beta2*data2$treat[i]+b[i,1]))
  a2<- 1- exp(alpha*(beta1+beta3*data2$treat[i]+b[i,2])*data2$with.time[i])
  log.survival[i]<-a1*a2
}
  
  
  log.p.tb<-log.hazard+log.survival
#print(log.p.yt)
   lgLik[t]<-loglik.yb+sum(log.p.tb)
# print(lgLik[t])

   if (lgLik[t] > cur.lgLik) {  { EM.lambda<-lambda} &{ cur.lgLik<-lgLik[t]  } }
else  {  {cur.lgLik <-cur.lgLik} &{break} }
# if (t==2) {break}
# update new parameter values

   S1.lambda<-rep(0,n)
   S2.lambda<-rep(0,n)

   for (i in 1:n)  {
a1<-exp(gamma*data2$treat[i]+alpha*(beta0+beta2*data2$treat[i]+b[i,1]))
b1<-alpha*(beta1+beta3*data2$treat[i]+b[i,2])
c1<- 1- exp(b1*data2$with.time[i])

{S1.lambda[i]<-(delta[i]/lambda)+(a1/b1)*c1} &
{S2.lambda[i]<- (-delta[i])/(lambda^2)} }
   
   
Sum.S1.lambda<-sum(S1.lambda,na.rm=TRUE)
Sum.S2.lambda<-sum(S2.lambda,na.rm=TRUE)

nlambda<- lambda-Sum.S1.lambda/Sum.S2.lambda
#print(nlambda)
if (nlambda< 0) {break}
lambda<-nlambda
}

EM.lambda<-lambda

a<-c(EM.lambda,EM.gamma,EM.alpha,cur.lgLik)

a
}

#twostage.JMLC.lambda(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],lmeFit$coefficients$fixed[3],lmeFit$coefficients$fixed[4],0.1,0.1,0.02,b,loglik.yb,delta1,10,Data,data2)


##################### Function of update gamma

twostage.JMLC.gamma<-function(beta0,beta1,beta2,beta3,lambda,gamma,alpha,b,loglik.yb,delta,iters,Data,data2){
N<-dim(Data)[1]
N
# number of measurement accasions for each subject
ni<-as.vector(tapply(Data$id,Data$id,length))
ni
summary(ni)
data2<- Data[!duplicated(Data$id),]

# number of subjects in the sample
n<-length(data2$id)
n
loglik.yb<-logLik(lmeFit)
lgLik<-rep(numeric(),iters)

X1<-rep(1,n)
X2<-data2$with.time


Z.time<-matrix(c(X1,X2),ncol=2)
Ztime.b<- rep(0,n)

 for (i in 1:n) {Ztime.b[i]<- Z.time[i,] %*% t(b[i,])}


cur.lgLik<- -10000

EM.lambda<-lambda
EM.gamma<-gamma
EM.alpha<-alpha

for (t in 1:iters)  {
log.hazard<-delta*log(lambda) + delta*( gamma*data2$treat + alpha*(beta0 + beta1*data2$with.time+beta2*data2$treat
+beta3*(data2$treat*data2$with.time)+Ztime.b))
# print(log.hazard)
log.survival<-rep(0,n)

   for (i in 1:n) {
a1<-(lambda/(alpha*(beta1+beta3*data2$treat[i]+b[i,2])))*exp(gamma*data2$treat[i]+alpha*(beta0+beta2*data2$treat[i] +b[i,1]))
a2<- 1- exp(alpha*(beta1+beta3*data2$treat[i]+b[i,2])*data2$with.time[i])
log.survival[i]<-a1*a2
 }
  
log.p.tb<-log.hazard+log.survival
#print(log.p.yt)
 
lgLik[t]<-loglik.yb+sum(log.p.tb)
# print(lgLik[t])

if (lgLik[t] > cur.lgLik) {  { EM.gamma<-gamma} &{ cur.lgLik<-lgLik[t]  } }
else  {  {cur.lgLik <-cur.lgLik} &{break} }
#    if (t==2) {break}
# update new parameter values
  
S1.gamma<-rep(0,n)
S2.gamma<-rep(0,n)


for (i in 1:n)  {
 a1<-lambda*exp(gamma*data2$treat[i]+alpha*(beta0+beta2*data2$treat[i]+b[i,1]))
 b1<-alpha*(beta1+beta3*data2$treat[i]+b[i,2])
 c1<- 1- exp(b1*data2$with.time[i])
  
 
 {S1.gamma[i]<-delta[i]*data2$treat[i]+(a1/b1)*c1*data2$treat[i]} &
 {S2.gamma[i]<-(a1/b1)*c1*(data2$treat[i]^2)}}


Sum.S1.gamma<-sum(S1.gamma,na.rm=TRUE)
Sum.S2.gamma<-sum(S2.gamma,na.rm=TRUE)
  


ngamma<- gamma-Sum.S1.gamma/Sum.S2.gamma
#print(ngamma)

if (ngamma< 0) {break}
gamma<-ngamma
}

EM.gamma<-gamma

#print(cur.lgLik)
#print(lgLik)
a<-c(EM.lambda,EM.gamma,EM.alpha,cur.lgLik)
a
}

#twostage.JMLC.gamma(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],lmeFit$coefficients$fixed[3],lmeFit$coefficients$fixed[4],0.1,-1,0.02,b,loglik.yb,delta1,10,Data,data2)

##################### Function of update alpha
twostage.JMLC.alpha<-function(beta0,beta1,beta2,beta3,lambda,gamma,alpha,b,loglik.yb,delta,iters,Data,data2){
N<-dim(Data)[1]
N
# number of measurement accasions for each subject
ni<-as.vector(tapply(Data$id,Data$id,length))
ni
summary(ni)
data2<- Data[!duplicated(Data$id),]
# number of subjects in the sample
n<-length(data2$id)
n
loglik.yb<-logLik(lmeFit)
lgLik<-rep(numeric(),iters)

X1<-rep(1,n)
X2<-data2$with.time

Z.time<-matrix(c(X1,X2),ncol=2)
Ztime.b<- rep(0,n)

for (i in 1:n) {Ztime.b[i]<- Z.time[i,] %*% t(b[i,])}

cur.lgLik<- -10000
EM.lambda<-lambda
EM.gamma<-gamma
EM.alpha<-alpha

for (t in 1:iters)  {
log.hazard<-delta*log(lambda) + delta*( gamma*data2$treat + alpha*(beta0 + beta1*data2$with.time+beta2*data2$treat
+beta3*(data2$treat*data2$with.time)+Ztime.b))
# print(log.hazard)

log.survival<-rep(0,n)

for (i in 1:n) {
a1<-(lambda/(alpha*(beta1+beta3*data2$treat[i]+b[i,2])))*exp(gamma*data2$treat[i]+alpha*(beta0+beta2*data2$treat[i] +b[i,1]))
a2<- 1- exp(alpha*(beta1+beta3*data2$treat[i]+b[i,2])*data2$with.time[i])
log.survival[i]<-a1*a2
}

log.p.tb<-log.hazard+log.survival
#print(log.p.yt)

lgLik[t]<-loglik.yb+sum(log.p.tb)
# print(lgLik[t])

if (lgLik[t] > cur.lgLik) {  { EM.alpha<-alpha} &{ cur.lgLik<-lgLik[t]  } }
else  {  {cur.lgLik <-cur.lgLik} &{break} }
#    if (t==2) {break}

# update new parameter values

S1.alpha<-rep(0,n)
S2.alpha<-rep(0,n)

for (i in 1:n)  {
f1<-function(x){
          vara1<-lambda*exp(gamma*data2$treat[i]+x*(beta0+beta2*data2$treat[i]+b[i,1]))
          varb1<-x*(beta1+beta3*data2$treat[i]+b[i,2])
          varc1<- 1-exp(varb1*data2$with.time[i])

 fun<-(vara1/varb1)*varc1   }

 e1<- beta0 + beta1*data2$with.time[i]+beta2*data2$treat[i]+
beta3*(data2$treat[i]*data2$with.time[i])+b[i,1]+b[i,2]*data2$with.time[i]

        
        {S1.alpha[i]<-delta[i]*e1+grad(f1,alpha)} &
        {S2.alpha[i]<-hessian(f1,alpha)}                    }
      
      Sum.S1.alpha<-sum(S1.alpha,na.rm=TRUE)
      Sum.S2.alpha<-sum(S2.alpha,na.rm=TRUE)
      
      nalpha<- alpha-Sum.S1.alpha/Sum.S2.alpha

#print(nalpha)

   alpha<-nalpha
      
    }
    EM.alpha<-alpha
    
    
    a<-c(EM.lambda,EM.gamma,EM.alpha,cur.lgLik)
    a
    
  }
  
#twostage.JMLC.alpha(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],lmeFit$coefficients$fixed[3],lmeFit$coefficients$fixed[4],0.1,0.15,-1,b,loglik.yb,delta1,10,Data,data2)
 

######################################function of Updtae parameter

update.parameter<-function(beta0,beta1,beta2,beta3,lambda,gamma,alpha,b,loglik.yb,delta,Data,data2,iters){
  #update.parameter(0.5,1,1.5,1,1,0.05,0.1,0.1,0.01,1,3,3,3,3,2,10)
  
    f.lambda<-rep(0,iters+1)
    
    f.gamma<-rep(0,iters+1)
    f.alpha<-rep(0,iters+1)
    
    f.lgLik<-rep(0,iters+1)
  
   
  f.lambda[1]<-lambda
  
  f.gamma[1]<-gamma
  f.alpha[1]<-alpha
  ########
  
  N<-dim(Data)[1]
  N
  # number of measurement accasions for each subject
  ni<-as.vector(tapply(Data$id,Data$id,length))
  ni
  #summary(ni)
  data2<- Data[!duplicated(Data$id),]
  # number of subjects in the sample
  n<-length(data2$id)
  n
  
  X1<-rep(1,n)
  X2<-data2$with.time
    
  Z.time<-matrix(c(X1,X2),ncol=2)
  
  Ztime.b<- rep(0,n) 
  for (i in 1:n) {Ztime.b[i]<- Z.time[i,] %*% t(b[i,])}
  
    
  log.hazard<- delta*log(lambda)+delta*( gamma*data2$treat + alpha*(beta0 + beta1*data2$with.time+beta2*data2$treat
                                                   +beta3*(data2$treat*data2$with.time)+Ztime.b))
  # print(log.hazard)
  
  log.survival<-rep(0,n)
  for (i in 1:n) {
    a1<-(f.lambda[1]/(f.alpha[1]*(beta1+beta3*data2$treat[i]+b[i,2])))*exp(f.gamma[1]*data2$treat[i]+f.alpha[1]*(beta0+beta2*data2$treat[i]+b[i,1]))
    a2<- 1- exp(f.alpha[1]*(beta1+beta3*data2$treat[i]+b[i,2])*data2$with.time[i])
    
    
    log.survival[i]<-a1*a2 
  }
    
  log.p.tb<-log.hazard+log.survival 
  
  f.lgLik[1]<-loglik.yb+sum(log.p.tb)
  
  # print(lgLik[t])
  
  ############
  for (i in 1:iters) {
        
  re1<-twostage.JMLC.lambda(beta0,beta1,beta2,beta3,lambda,gamma,alpha,b,loglik.yb,delta,iters,Data,data2)
    
    lambda<-re1[1]
    lgLik1<-re1[4]
    
    
    
  re2<-twostage.JMLC.gamma(beta0,beta1,beta2,beta3,lambda,gamma,alpha,b,loglik.yb,delta,iters,Data,data2)
    
    gamma<-re2[2]
    lgLik2<-re2[4]
    
  re3<-twostage.JMLC.alpha(beta0,beta1,beta2,beta3,lambda,gamma,alpha,b,loglik.yb,delta,iters,Data,data2)
      
      alpha<-re3[3]
      
      lgLik3<-re3[4]
      
    f.lambda[i+1]<-lambda

    f.gamma[i+1]<-gamma
    f.alpha[i+1]<-alpha
    
    epsilon<-10^(-8)
    
     
      if ( lgLik3 -f.lgLik[i] < epsilon*(abs(f.lgLik[i])+epsilon)) {break}
      #
      f.lgLik[i+1]<-lgLik3   
      
  }
  #print(f.lambda)  
  
  #print(f.gamma)
  
  #print(f.alpha)
  
  #print(f.lgLik)
 final.parameter<-c(lambda,gamma,alpha)
    final.parameter
}
############ Likelihood function(Compute AIC,BIC)

lk<-function(beta0,beta1,beta2,beta3,lambda1,lambda2,gamma1,gamma2,alpha1,alpha2,b,loglik.yb,delta1,delta2,iters,Data,data2){
N<-dim(Data)[1]
N
# number of measurement accasions for each subject
ni<-as.vector(tapply(Data$id,Data$id,length))
ni
summary(ni)
data2<- Data[!duplicated(Data$id),]
# number of subjects in the sample
n<-length(data2$id)
n

loglik.yb<-logLik(lmeFit)

lgLik<-rep(numeric(),iters)

X1<-rep(1,n)
X2<-data2$with.time

Z.time<-matrix(c(X1,X2),ncol=2)

Ztime.b<- rep(0,n)
for (i in 1:n) {Ztime.b[i]<- Z.time[i,] %*% t(b[i,])}


  log.hazard1<-delta1*log(lambda1) + delta1*( gamma1*data2$treat + alpha1*(beta0 + beta1*data2$with.time+beta2*data2$treat
  +beta3*(data2$treat*data2$with.time)+Ztime.b))


 log.hazard2<-delta2*log(lambda2) + delta2*( gamma2*data2$treat + alpha2*(beta0 + beta1*data2$with.time+beta2*data2$treat
  +beta3*(data2$treat*data2$with.time)+Ztime.b))
# print(log.hazard2)

  log.survival1<-rep(0,n)
  
for (i in 1:n) {
  a1<-(lambda1/(alpha1*(beta1+beta3*data2$treat[i]+b[i,2])))*exp(gamma1*data2$treat[i]+alpha1*(beta0+beta2*data2$treat[i]+b[i,1]))
  a2<- 1- exp(alpha1*(beta1+beta3*data2$treat[i]+b[i,2])*data2$with.time[i])
  log.survival1[i]<-a1*a2
}
  
  
  log.survival2<-rep(0,n)
  
for (i in 1:n) {
  a1<-(lambda2/(alpha2*(beta1+beta3*data2$treat[i]+b[i,2])))*exp(gamma2*data2$treat[i]+alpha2*(beta0+beta2*data2$treat[i]+b[i,1]))
  a2<- 1- exp(alpha2*(beta1+beta3*data2$treat[i]+b[i,2])*data2$with.time[i])
  log.survival2[i]<-a1*a2
}
  
  log.p.tb<-log.hazard1+log.hazard2+log.survival1+log.survival1
#print(log.p.yt)
   lgLik<-loglik.yb+sum(log.p.tb)
 lgLik
}

###############  Cause ISC
time00 = Sys.time()
surv.par1<-update.parameter(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],lmeFit$coefficients$fixed[3],lmeFit$coefficients$fixed[4],0.1,0.015,0.02,b,loglik.yb,delta1,Data,data2,100)
surv.par1


lambda1<-round(surv.par1[1],3)
gamma1<-round(surv.par1[2],3)
alpha1<-round(surv.par1[3],3)

################ Cause UAE
surv.par2<-update.parameter(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],lmeFit$coefficients$fixed[3],lmeFit$coefficients$fixed[4],0.1,-0.46,0.02,b,loglik.yb,delta2,Data,data2,100)
surv.par2
lambda2<-round(surv.par2[1],3)
gamma2<-round(surv.par2[2],3)
alpha2<-round(surv.par2[3],3)
Finaltime=Sys.time() - time00
Finaltime
####
est<-c(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],
lmeFit$coefficients$fixed[3],lmeFit$coefficients$fixed[4],
sigma,sigma11,sigma22,sigma12,lambda1,lambda2,gamma1,gamma2,alpha1,alpha2)

LG.lIK<-lk(est[1],est[2],est[3],est[4],est[9],est[10],est[11],est[12],est[13],est[14],b,loglik.yb,delta1,delta2,Data,data2)
AIC<--2*LG.lIK+2*14
BIC<- -2*LG.lIK+14*log(14)
AIC
BIC

########################
confidence.interval<-function(beta0,beta1,beta2,beta3,lambda,gamma,alpha,b,delta,Data,data2){

data2<- Data[!duplicated(Data$id),]
# number of subjects in the sample
n<-length(data2$id)
n

S2.lambda<-rep(0,n)
S2.gamma<-rep(0,n)
S2.alpha<-rep(0,n)

   for (i in 1:n)  {

f1<-function(x){
          vara1<-lambda*exp(gamma*data2$treat[i]+x*(beta0+beta2*data2$treat[i]+b[i,1]))
          varb1<-x*(beta1+beta3*data2$treat[i]+b[i,2])
          varc1<- 1-exp(varb1*data2$with.time[i])

 fun<-(vara1/varb1)*varc1   }
 
S2.gamma[i]<-f1(alpha)*(data2$treat[i]^2)

S2.lambda[i]<- (-delta[i])/(lambda^2)
S2.alpha[i]<-hessian(f1,alpha) 

 }
Sum.S2.lambda<-sum(S2.lambda,na.rm=TRUE)
Sum.S2.gamma<-sum(S2.gamma,na.rm=TRUE)
Sum.S2.alpha<-sum(S2.alpha,na.rm=TRUE)

v1<- -1/Sum.S2.lambda
v2<- -1/Sum.S2.gamma
v3<- -1/Sum.S2.alpha
int.lambda<-c(lambda-1.96*sqrt(v1),lambda+1.96*sqrt(v1))
int.gamma<-c(gamma-1.96*sqrt(v2),gamma+1.96*sqrt(v2))
int.alpha<-c(alpha-1.96*sqrt(v3),alpha+1.96*sqrt(v3))

interval<-cbind(int.lambda,int.gamma,int.alpha)
interval
}
confidence.interval(est[1],est[2],est[3],est[4],lambda1,gamma1,alpha1,b,delta1,Data,data2)
confidence.interval(est[1],est[2],est[3],est[4],lambda2,gamma2,alpha2,b,delta2,Data,data2)


#####################################################
################# SEPARATE MODEL

lmeFit <- lme(dose ~ time+treat * time,
              random = ~ time | id,
              data = epileptic)

##### Confidence interval
fixef(lmeFit) - qnorm(0.975) * sqrt(diag(lmeFit$varFix))
fixef(lmeFit) + qnorm(0.975) * sqrt(diag(lmeFit$varFix))
intervals(lmeFit)

####################

survdat <- epileptic[!duplicated(epileptic$id), ]
survdat$status <- rep("alive", nrow(survdat))
survdat$status[survdat$with.status.uae == 1] <- "UAE"
survdat$status[survdat$with.status.isc == 1] <- "ISC"
survdat.CR <- crLong(survdat, statusVar = "status",
                     censLevel = "alive", nameStrata = "CR")

coxFit <- coxph(Surv(with.time, status2) ~ treat*CR + strata(CR),
                data = survdat.CR, 
                x = TRUE)

summary(coxFit)
confint(coxFit)
aic.sep<-AIC(lmeFit)+AIC(coxFit)
bic.sep<-BIC(lmeFit)+BIC(coxFit)

#########Rizopoulos Joint model ###########################

time00<-Sys.time()
jointFit <- jointModel(lmeFit, coxFit,
                        timeVar = "time",
                        method = "spline-PH-GH",
                        CompRisk = TRUE,
                        interFact = list(value = ~ CR, data = survdat.CR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        verbose = FALSE)

summary(jointFit)
Finaltime=Sys.time() - time00

jointFit$coef
confint(jointFit)

confint(jointFit,parm="Longitudinal")
#########################
 
coxFit.isc <- coxph(Surv(with.time, with.status.isc) ~ treat,
                    data = survdata,x=TRUE)
coxFit.uae<- coxph(Surv(with.time, with.status.uae) ~ treat,
                    data = survdata,x=TRUE)
summary(coxFit.isc)

summary(coxFit.uae)

confint(coxFit.isc)
confint(coxFit.uae) 
   

