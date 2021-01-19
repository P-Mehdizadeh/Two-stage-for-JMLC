library(parallel)
library(mnormt)
library(MASS)
library(mvtnorm)
library(survival)
library(nlme)
library(joineR)
library(JM)
library(splines)
library(grid)

rep<-200

RESULTS=matrix(0,rep,12)
SEPARATE=matrix(0,rep,2)
JOINT=matrix(0,rep,10)
#time00 = Sys.time()
colnames(RESULTS)=c("beta0","beta1","sigma","D11","D12","D22","lambda1","lambda2","gamma1","gamma2","alpha1","alpha2")
#RESULTS2=matrix(0,10,3)
#generate data
  simulation.JMLC.model<-function(n,beta0,beta1,lambda1,lambda2,gamma1,gamma2,alpha1,alpha2,sigma,D11,D12,D22) {
    
    # n: number of subject
    # sigma =var(epsilon)
    
    # crreate survival time
    
    D<-matrix(0,2,2)
    #D=diag(.5,4)+matrix(.5,4,4)
    D[1,1]<-D11
    D[1,2]<-D12
    D[2,1]<-D12
    D[2,2]<-D22
    
    #print(D)
    r<-matrix(0,n,2)
    
    # print(r)
    
    group.n<-c(0,n)
    u1<-c(0,n)
    u2<-c(0,n)
    
    Time.p1=Time.p2<-rep(0,n)
    Time.p<-rep(0,n)
    for (i in 1:n) {
  repeat{
      r[i,]<- rmvnorm(1,rep(0,2),D)
      group.n[i]<-rbinom(1,1,0.5)
      u1[i]<-runif(1,0,1)
      u2[i]<-runif(1,0,1)
      
      a11<-lambda1*exp(gamma1*group.n[i]+alpha1*(beta0+r[i,1]))
      
      b11<-alpha1*(beta1+r[i,2])
      
      Time.p1[i]<-(1/b11)*log((b11/a11)*(-log(u1[i]))+1)
      
      a12<-lambda2*exp(gamma2*group.n[i]+alpha2*(beta0+r[i,1]))
      
      b12<-alpha2*(beta1+r[i,2])
      
      
      Time.p2[i]<-(1/b12)*log((b12/a12)*(-log(u2[i]))+1)
      
      Time.p[i]<-min(Time.p1[i],Time.p2[i])
      
            if(!is.na(Time.p[i])) break
}

    }                                                                                                                      
    
   
    #############
    
    
    # print(Time.p)
    # the censoring mechanism is induced by an exponential process with parameter 0.25
    censor.Time<-rexp(n,0.25)
    #  print(censor.Time)
    
    
    Time.n<-rep(0,n)
    death.n<-rep(0,n)
    D<-rep(0,n)
    D1<-rep(0,n)
    surt<-matrix(0,nrow=n,ncol=2)
    surt.cen<-matrix(0,nrow=n,ncol=2)
    for (i in 1:n) {
      Time.n[i]<-min(Time.p[i],censor.Time[i])
      if (Time.p[i]<=censor.Time[i]) {death.n[i]<-1}
      if( Time.p1[i]< Time.p2[i]) {D[i]=1} else {D[i]=2}
      if (Time.p[i]<=censor.Time[i]) {D1[i]=D[i]} else {D1[i]=0}  
      if (D1[i]== 1)  {surt[i,1]=Time.n[i];surt.cen[i,1]=0 } else {surt[i,1]=NA; surt.cen[i,1]=Time.n[i]}
      if (D1[i]== 2)  {surt[i,2]=Time.n[i];surt.cen[i,2]=0 } else {surt[i,2]=NA; surt.cen[i,2]=Time.n[i]}
    }
    table(D1)
    #  print(Time.n) 
    
    #print(sum(death.n))
    
    # create obstime at 0,0.5,1,1.5,2,...
    b<-0
    for (i in 1:n) { 
      if (abs(Time.n[i]-floor(Time.n[i]))<0.5) {a<-2*floor(Time.n[i])}
      else {a <- 2*floor(Time.n[i])+1}
      b  <- b+a }
    
    
    #print(b)
    
    # create variable obstime
    
    obstime <- rep(0,b+n)
    Time <- rep(0,b+n)
    group <- rep(0,b+n)
    death <- rep(0,b+n)
    patient <- rep(0,b+n)
   
    a<-rep(0,(n+b)*2)
    random<-matrix(a,ncol=2)
    s<-rep(0,(n+b)*2)
    surt1<-matrix(s,ncol=2)
    cause<-rep(0,b+n)
    x <- 0
    y <- 0
    for (i in 1:n) {
      if (abs(Time.n[i]-floor(Time.n[i]))<0.5) {ni<-floor(Time.n[i])}
      else {ni<-floor(Time.n[i])+0.5}
      
      for (j in 1: (2*ni+1)) {
        obstime[y+j] <- (j/2)-0.5  
        Time[y+j] <- Time.n[i]
        group[y+j] <- group.n[i]
        death[y+j] <- death.n[i]
        patient[y+j] <- i   
        random[y+j,1] <-r[i,1]
        random[y+j,2] <-r[i,2]
        surt1[y+j,1]<-surt[i,1]
        surt1[y+j,2]<-surt[i,2]
        cause[y+j]<-D1[i]
      }
      x <- 2*ni
      y <- y+x+1
      
    }
    #print(head(random,20))
    
    
    #create epsilon variable
    epsilon<-rnorm(n+b,0,sigma)
    
    #create longitudinal data
    longi<-rep(0,n+b)
    for (i in 1:(n+b)) {
      longi[i]<-beta0+beta1*obstime[i]+random[i,1]+random[i,2]*obstime[i]+epsilon[i]
    }
    
    
    data<-data.frame(patient,obstime,Time,group,longi,death,cause,surt1[,1],surt1[,2])
    #print(head(data,100))
    response<-cbind(data,1-mean(death.n))
    response
  }
#
data1<-simulation.JMLC.model(500,5,2,0.2,0.4,1,0.5,0.1,0.08,2,1,0.5,1)
########### function of update lambda

 twostage.JMLC.lambda<-function(beta0,beta1,lambda,gamma,alpha,b,loglik.yb,iters,delta,data1,data2){
    
    # number of measurement occasions for the whole sample  
    N<-length(data1$patient)
    
    # number of measurement accasions for each subject
    ni<-as.vector(tapply(data1$patient,data1$patient,length))
    
    # data22<-subset(data3,obstime==0)
    # number of subjects in the sample
    n<-length(data2$patient)
    
    loglik.yb<-logLik(lmeFit)
    
    lgLik<-rep(numeric(),iters) 
    
    
    X1<-rep(1,n)
    X2<-data2$Time
    
    
    Z.time<-matrix(c(X1,X2),ncol=2)
    
    Ztime.b<- rep(0,n) 
    for (i in 1:n) {Ztime.b[i]<- Z.time[i,] %*% t(b[i,])}
    
    
    group<-data2$group  
    
    cur.lgLik<- -10000  
    
    EM.lambda<-lambda 
    EM.gamma<-gamma
    EM.alpha<-alpha
    
    
    for (t in 1:iters)  {
      
      
      log.hazard<-delta*log(lambda) + delta*( gamma*group+ alpha*(beta0 + beta1*data2$Time+Ztime.b))
      # print(log.hazard)
      
      log.survival<-rep(0,n)
      for (i in 1:n) {
        a1<-(lambda/(alpha*(beta1+b[i,2])))*exp(gamma*group[i]+alpha*(beta0+b[i,1]))
        a2<- 1- exp(alpha*(beta1+b[i,2])*data2$Time[i])
        
        
        log.survival[i]<-a1*a2    }
      
      
      log.p.tb<-log.hazard+log.survival 
      
      
      #print(log.p.yt)
      lgLik[t]<-loglik.yb+sum(log.p.tb)
      
      # print(lgLik[t])
      
      
      
      
      if (lgLik[t] > cur.lgLik) {  { EM.lambda<-lambda} &{ cur.lgLik<-lgLik[t]  } }
      else  {  {cur.lgLik <-cur.lgLik} &{break} } 
      
      
      S1.lambda<-rep(0,n)
      S2.lambda<-rep(0,n)
      
      for (i in 1:n)  {
        
        a1<-exp(gamma*group[i]+alpha*(beta0+b[i,1]))    
        b1<-alpha*(beta1+b[i,2])    
        c1<- 1- exp(b1*data2$Time[i])   
        
        
        
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
    #a<-beta0.his
    
    a
  }
########## Function of update gamma 

  twostage.JMLC.gamma<-function(beta0,beta1,lambda,gamma,alpha,b,loglik.yb,iters,delta,data1,data2){
    
    # number of measurement occasions for the whole sample  
    N<-length(data1$patient)
    
    # number of measurement accasions for each subject
    ni<-as.vector(tapply(data1$patient,data1$patient,length))
    
    #  data2<-subset(data1,obstime==0)
    # number of subjects in the sample
    n<-length(data2$patient)
    
    lgLik<-rep(numeric(),iters) 
    
    X1<-rep(1,n)
    X2<-data2$Time
    
    
    Z.time<-matrix(c(X1,X2),ncol=2)
    
    Ztime.b<- rep(0,n) 
    for (i in 1:n) {Ztime.b[i]<- Z.time[i,] %*% t(b[i,])}
    
    
    group<-data2$group  
        
    cur.lgLik<- -100000  
        
    EM.lambda<-lambda  
    EM.gamma<-gamma
    EM.alpha<-alpha
    
        for (t in 1:iters)  {
      
      log.hazard<-delta*log(lambda) + delta*( gamma*group+ alpha*(beta0 + beta1*data2$Time+Ztime.b))
      # print(log.hazard)
      
      log.survival<-rep(0,n)
      for (i in 1:n) {
        a1<-(lambda/(alpha*(beta1+b[i,2])))*exp(gamma*group[i]+alpha*(beta0+b[i,1]))
        a2<- 1- exp(alpha*(beta1+b[i,2])*data2$Time[i])
        
        
        log.survival[i]<-a1*a2 }
      
      
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
        a1<-lambda*exp(gamma*group[i]+alpha*(beta0+b[i,1]))      
        b1<-alpha*(beta1+b[i,2])        
        c1<- 1- exp(b1*data2$Time[i])
        
        
        
        {S1.gamma[i]<-delta[i]*group[i]+(a1/b1)*c1*group[i]} &
        {S2.gamma[i]<-(a1/b1)*c1*(group[i]^2)}}
      
      
      
      Sum.S1.gamma<-sum(S1.gamma,na.rm=TRUE)
      Sum.S2.gamma<-sum(S2.gamma,na.rm=TRUE)
      
      ngamma<- gamma-Sum.S1.gamma/Sum.S2.gamma
      
      #print(ngamma)
      
      if (ngamma< 0) {break}
      gamma<-ngamma
      
    }
    EM.gamma<-gamma
    
    a<-c(EM.lambda,EM.gamma,EM.alpha,cur.lgLik)
    
    a
  }
  
  
###twostage.JMLC.gamma(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],0.1,0.15,0.02,b,loglik.yb,10,data1,data2)
  
  
############## Function of Update alpha 
  
 twostage.JMLC.alpha<-function(beta0,beta1,lambda,gamma,alpha,b,loglik.yb,iters,delta,data1,data2){
    
    # number of measurement occasions for the whole sample  
    N<-length(data1$patient)
    N
    # number of measurement accasions for each subject
    ni<-as.vector(tapply(data1$patient,data1$patient,length))
    ni
    
    #  data2<-subset(data1,obstime==0)
    # number of subjects in the sample
    n<-length(data2$patient)
    
    #iters<-10
    lgLik<-rep(numeric(),iters) 
        
    X1<-rep(1,n)
    X2<-data2$Time
        
    Z.time<-matrix(c(X1,X2),ncol=2)
    
    Ztime.b<- rep(0,n) 
    for (i in 1:n) {Ztime.b[i]<- Z.time[i,] %*% t(b[i,])}
    
    group<-data2$group  
    
    cur.lgLik<- -100000  
    
    EM.lambda<-lambda  
    EM.gamma<-gamma
    EM.alpha<-alpha
    
        
    for (t in 1:iters)  {
           
      log.hazard<-delta*log(lambda) + delta*( gamma*group+ alpha*(beta0 + beta1*data2$Time+Ztime.b))
      # print(log.hazard)
      
      log.survival<-rep(0,n)
      for (i in 1:n) {
        a1<-(lambda/(alpha*(beta1+b[i,2])))*exp(gamma*group[i]+alpha*(beta0+b[i,1]))
        a2<- 1- exp(alpha*(beta1+b[i,2])*data2$Time[i])
        
        
        log.survival[i]<-a1*a2
      }
            
      log.p.tb<-log.hazard+log.survival 
            
      
      #print(log.p.yt)
      lgLik[t]<-loglik.yb+sum(log.p.tb)
      
      # print(lgLik[t])
      
      
      if (lgLik[t] > cur.lgLik) {  { EM.alpha<-alpha} &{ cur.lgLik<-lgLik[t]  } }
      else  {  {cur.lgLik <-cur.lgLik} &{break}  } 
      #  if (t==2) {break}
      # update new parameter values
      
      
      S1.alpha<-rep(0,n)
      S2.alpha<-rep(0,n)
      
      for (i in 1:n)  {      
        
        
        f1<-function(x){
          vara1<-lambda*exp(gamma*group[i]+x*(beta0+b[i,1]))
          varb1<-x*(beta1+b[i,2])
          varc1<- 1-exp(varb1*data2$Time[i])
          
          fun<-(vara1/varb1)*varc1   }
        
        
        
        e1<- beta0+beta1*data2$Time[i]+b[i,1]+b[i,2]*data2$Time[i]
        
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
  
  #twostage.JMLC.alpha(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],0.1,0.15,0.02,b,loglik.yb,10,data1,data2)
  

########################################### function of update parameters  

  update.parameter<-function(beta0,beta1,lambda,gamma,alpha,b,loglik.yb,delta,data1,data2,iters){
    
    f.lambda<-rep(0,iters+1)
    
    f.gamma<-rep(0,iters+1)
    f.alpha<-rep(0,iters+1)
    
    f.lgLik<-rep(0,iters+1)
    
    
    # lambda<-ini.lambda
    
    # gamma<-ini.gamma
    # alpha<-ini.alpha
    
    
    f.lambda[1]<-lambda
    
    f.gamma[1]<-gamma
    f.alpha[1]<-alpha
    ########
    N<-length(data1$patient)
    N
    # number of measurement accasions for each subject
    ni<-as.vector(tapply(data1$patient,data1$patient,length))
    ni
    
    # number of subjects in the sample
    n<-length(data2$patient)
    
    
    X1<-rep(1,n)
    X2<-data2$Time
    
    
    Z.time<-matrix(c(X1,X2),ncol=2)
    
    Ztime.b<- rep(0,n) 
    for (i in 1:n) {Ztime.b[i]<- Z.time[i,] %*% t(b[i,])}
    
    
    group<-data2$group  
    
    log.hazard<-delta*log(lambda) + delta*( gamma*group+ alpha*(beta0 + beta1*data2$Time+Ztime.b))
    # print(log.hazard)
    
    log.survival<-rep(0,n)
    for (i in 1:n) {
      a1<-(f.lambda[1]/(f.alpha[1]*(beta1+b[i,2])))*exp(f.gamma[1]*group[i]+f.alpha[1]*(beta0+b[i,1]))
      a2<- 1- exp(f.alpha[1]*(beta1+b[i,2])*data2$Time[i])
      
      
      log.survival[i]<-a1*a2
    }
    
    
    log.p.tb<-log.hazard+log.survival 
    
    
    f.lgLik[1]<-loglik.yb+sum(log.p.tb)
    
    # print(lgLik[t])
    
    ############
    for (i in 1:iters) {
      
      
      re1<-twostage.JMLC.lambda(beta0,beta1,lambda,gamma,alpha,b,loglik.yb,iters,delta,data1,data2)
      
      
      lambda<-re1[1]
      lgLik1<-re1[4]
      
      
      
      re2<-twostage.JMLC.gamma(beta0,beta1,lambda,gamma,alpha,b,loglik.yb,iters,delta,data1,data2)
      
      gamma<-re2[2]
      lgLik2<-re2[4]
      
      re3<-twostage.JMLC.alpha(beta0,beta1,lambda,gamma,alpha,b,loglik.yb,iters,delta,data1,data2)
      
      alpha<-re3[3]
      
      lgLik3<-re3[4]
      
     
      f.lambda[i+1]<-lambda
      
      f.gamma[i+1]<-gamma
      f.alpha[i+1]<-alpha
      
      epsilon<-10^(-8)
      
      if ( abs(lgLik3 -f.lgLik[i]) < epsilon*(abs(f.lgLik[i])+epsilon)) {break}
     
      f.lgLik[i+1]<-lgLik3   

    }
    
    #print(f.lambda)  
    
    #print(f.gamma)
    
    #print(f.alpha)
    
    #print(f.lgLik)
    final.parameter<-c(lambda,gamma,alpha)
    final.parameter
    
  }
  
  ######################Generate N=200 sample for simulation
  out1<-list()
  
  for(kkk in 1:rep){
    #set.seed(kkk)
 out1[[kkk]]<-simulation.JMLC.model(500,5,2,0.2,0.4,1,0.5,0.1,0.08,2,1,0.5,1)
 print(kkk) 
}
result=out1
#save(result,file="D:\\phd\\Outputs\\result.RData")
#load(file.choose())
  

######calculate rate of censoring
rate<-c()
for(i in 1:200){
b<-out1[[i]]
rate[i]<-b[i,10]
}
 mean(rate)


################Estimate Parameters based on Two-Stage approach for JMLC and Separtae model
#############longitudinal process
sample<-result
for(kkk in 1:rep){
   data1=sample[[kkk]]
  ctrl<-lmeControl(opt='optim')
  lmeFit<-lme(longi~obstime,random=~obstime|patient,control=ctrl,data=data1)
  lmeFit
  b<-ranef(lmeFit)
  
  
  # number of measurement occasions for the whole sample  
  N<-length(data1$patient)
  N
  # number of measurement accasions for each subject
  ni<-as.vector(tapply(data1$patient,data1$patient,length))
  ni
  summary(ni)
  data2<-subset(data1,obstime==0)
  # number of subjects in the sample
  n<-length(data2$patient)
  n
  
survdat<-data2
survdat$status <- rep("alive", nrow(survdat))
survdat$status[survdat$cause == 1] <- "A"
survdat$status[survdat$cause == 2] <- "B"
survdat.CR <- crLong(survdat, statusVar = "status",
                     censLevel = "alive", nameStrata = "CR")

coxFit <- coxph(Surv(Time, status2) ~ group*CR + strata(CR),
                data = survdat.CR, 
                x = TRUE)

  loglik.yb<-logLik(lmeFit)
  #data2

  delta1<-c()
  for(i in 1:n){
    if(data2$cause[i]==1){delta1[i]=1} else {delta1[i]=0}
  } 
  
  delta2<-c()
  for(i in 1:n){
    if(data2$cause[i]==2){delta2[i]=1} else {delta2[i]=0}

  } 
 #################### competing risks process
  ##############cause 1
  
  
  surv.par1<-update.parameter(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],0.1,0.15,0.02,b,loglik.yb,delta1,data1,data2,100)
  surv.par1
  
  
  lambda1<-surv.par1[1]
  gamma1<-surv.par1[2]
  alpha1<-surv.par1[3]
  
  
  
  ########################cause 2
  
  surv.par2<-update.parameter(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],0.2,0.1,0.03,b,loglik.yb,delta2,data1,data2,100)
  surv.par2
  
  
  lambda2<-surv.par2[1]
  gamma2<-surv.par2[2]
  alpha2<-surv.par2[3]
  
  
  ######################
  
  RESULTS[kkk,]<-c(lmeFit$coefficients$fixed[1],lmeFit$coefficients$fixed[2],as.numeric(VarCorr(lmeFit)[6]),as.numeric(VarCorr(lmeFit)[1]),
                   as.numeric(VarCorr(lmeFit)[4])*as.numeric(VarCorr(lmeFit)[5])*as.numeric(VarCorr(lmeFit)[8]),as.numeric(VarCorr(lmeFit)[2]),
                   lambda1,lambda2,gamma1,gamma2,alpha1,alpha2 )
 
 SEPARATE[kkk,]<-c(coxFit$coef[1],coxFit$coef[1]+coxFit$coef[3])


  print(kkk)
  
  }



 #################Rizopolouse(2012) Joint model

for(kkk in 1:rep){
   data1=sample[[kkk]]
  ctrl<-lmeControl(opt='optim')
  lmeFit<-lme(longi~obstime,random=~obstime|patient,control=ctrl,data=data1)
  lmeFit
  #b<-ranef(lmeFit)
  N<-length(data1$patient)
  N
  # number of measurement accasions for each subject
  ni<-as.vector(tapply(data1$patient,data1$patient,length))
  ni
  summary(ni)
  data2<-subset(data1,obstime==0)
  # number of subjects in the sample
  n<-length(data2$patient)
  n
survdat<-data2
survdat$status <- rep("alive", nrow(survdat))
survdat$status[survdat$cause == 1] <- "A"
survdat$status[survdat$cause == 2] <- "B"
survdat.CR <- crLong(survdat, statusVar = "status",
                     censLevel = "alive", nameStrata = "CR")

coxFit <- coxph(Surv(Time, status2) ~ group*CR + strata(CR),
                data = survdat.CR, 
                x = TRUE)


jointFit <- jointModel(lmeFit, coxFit,
                        timeVar = "obstime",
                        method = "spline-PH-aGH",
                        CompRisk = TRUE,
                        interFact = list(value = ~ CR, data = survdat.CR),
                        #GHk = 9,
                        #GKk = 15,
                        #iter.qN = 500,
                        #numeriDeriv = "cd",
                        verbose = FALSE)
#summary(jointFit)
JOINT[kkk,]<-c(jointFit$coef$betas[1],jointFit$coef$betas[2],jointFit$coef$sigma,
jointFit$coef$D[1],jointFit$coef$D[2],jointFit$coef$D[4],jointFit$coef$gammas[1],
jointFit$coef$gammas[1]+jointFit$coef$gammas[2],jointFit$coef$alpha[1],
jointFit$coef$alpha[1]+jointFit$coef$alpha[2])
print(kkk)
}


#JOINT2 <- JOINT [which(rowSums(JOINT) != 0), ]

################
 Calculate Estimate, Bias, R.bias, RMSE and CP for Two-stage, Separate and Rizopoulos joint model


res<-matrix(0,rep,12)
res=RESULTS
############
Est<-c(mean(res[,1]),mean(res[,2]),mean(res[,3]),mean(res[,4]),mean(res[,5]),mean(res[,6]),mean(res[,7]),mean(res[,8]),mean(res[,9]),mean(res[,10]),mean(res[,11]),mean(res[,12]))
Est<-round(Est,3)
########theta:the vector of parameters
theta<-c(beta0,beta1,sigma,D11,D12,D22,lambda1,lambda2,gamma1,gamma2,alpha1,alpha2)
theta<-c(5,2,2,1,0.5,1,0.2,0.4,1,0.5,0.1,0.08)
###################
Bias<-c((mean(res[,1])-theta[1]),(mean(res[,2])-theta[2]),(mean(res[,3])-theta[3]),
        (mean(res[,4])-theta[4]),(mean(res[,5])-theta[5]),(mean(res[,6])-theta[6]),
        (mean(res[,7])-theta[7]),(mean(res[,8])-theta[8]),(mean(res[,9])-theta[9]),
        (mean(res[,10])-theta[10]),(mean(res[,11])-theta[11]),(mean(res[,12])-theta[12]))
Bias<-round(Bias,3)

R.Bias<-c((mean(res[,1])-theta[1])/theta[1],(mean(res[,2])-theta[2])/theta[2],(mean(res[,3])-theta[3])/theta[3],
          (mean(res[,4])-theta[4])/theta[4],(mean(res[,5])-theta[5])/theta[5],(mean(res[,6])-theta[6])/theta[6],
          (mean(res[,7])-theta[7])/theta[7],(mean(res[,8])-theta[8])/theta[8],(mean(res[,9])-theta[9])/theta[9],
          (mean(res[,10])-theta[10])/theta[10],(mean(res[,11])-theta[11])/theta[11],(mean(res[,12])-theta[12])/theta[12])
R.Bias<-round(R.Bias,3)
#RB<-R.Bias

##############
rmse<-c()
mse<-matrix(0,nrow=rep,ncol=12)

for(j in 1:12){
  for(i in 1:rep){
    mse[i,j]<-(res[i,j]-theta[j])^2
  }
  rmse[j]<-round(sqrt((1/rep)*sum( mse[,j])),3)
}
rmse
#############
SE<-c()
V<-c()
V<-c(var(res[,1]),var(res[,2]),var(res[,3]),var(res[,4]),
     var(res[,5]),var(res[,6]),var(res[,7]),var(res[,8]),
     var(res[,9]),var(res[,10]),var(res[,11]),var(res[,12]))
SE<-round(sqrt(V),3)
SE
############
al<-0.05
cp<-c()
m<-c()
tmat<-t(res)
low=up<-matrix(0,nrow=12,ncol=rep)
for(j in 1:12){
  m[j]<-0
  for(i in 1:rep){
    low[j,i]<-tmat[j,i]-SE[j]*qnorm(1-al/2)
    up[j,i]<-tmat[j,i]+SE[j]*qnorm(1-al/2)
    if ((low[j,i]<theta[j])&(theta[j]< up[j,i])) {m[j]=m[j]+1} else {m[j]=m[j]}
  }
  cp[j]<-m[j]/rep
}
cp


par<-c("beta0","beta1","sigma","D11","D12","D22","lambda1","lambda2","gamma1","gamma2","alpha1","alpha2")


table<-data.frame(par,theta,Est,Bias,R.Bias,SE,rmse,cp)
table
edit(table)

#################### separate parameter

sep<-matrix(0,rep,2)
sep=SEPERATE[1:200,]
est2<-c(mean(sep[,1]),mean(sep[,2]))
gamma<-c(gamma1,gamma2)
gamma<-c(1,0.5)

rmse2<-c()
mse2<-matrix(0,nrow=rep,ncol=2)

for(j in 1:2){
  for(i in 1:rep){
    mse2[i,j]<-(sep[i,j]-gamma[j])^2
  }
  rmse2[j]<-round(sqrt((1/rep)*sum( mse2[,j])),3)
}
rmse2<-round(rmse2,3)
##########

SE2<-c()
V2<-c()
V2<-c(var(sep[,1]),var(sep[,2]))
SE2<-round(sqrt(V2),3)
SE2
###########
bias2<-c(est2[1]-gamma[1],est2[2]-gamma[2])
rbias2<-c((est2[1]-gamma[1])/gamma[1],(est2[2]-gamma[2])/gamma[2])
bias2<-round(bias2,3)
############
al<-0.05
cp2<-c()
m2<-c()
tmat2<-t(sep)
low2=up2<-matrix(0,nrow=2,ncol=rep)
for(j in 1:2){
  m2[j]<-0
  for(i in 1:rep){
    low2[j,i]<-tmat2[j,i]-SE2[j]*qnorm(1-al/2)
    up2[j,i]<-tmat2[j,i]+SE2[j]*qnorm(1-al/2)
    if ((low2[j,i]<gamma[j])&(gamma[j]< up2[j,i])) {m2[j]=m2[j]+1} else {m2[j]=m2[j]}
  }
  cp2[j]<-m2[j]/rep
}
cp2<-round(cp2,3)
est2<-round(est2,3)
table2<-data.frame(gamma,est2,SE2,bias2,rbias2,rmse2,cp2)
table2
edit(table2)

####################joint model Rizopolous
JOINT2 <- JOINT [which(rowSums(JOINT) != 0), ]

joint<-matrix(0,rep,10)
joint=JOINT2
est3<-c(mean(joint[,1]),mean(joint[,2]),mean(joint[,3]),mean(joint[,4]),mean(joint[,5]),
mean(joint[,6]),mean(joint[,7]),mean(joint[,8]),mean(joint[,9]),mean(joint[,10]))
est3<-round(est3,3)

par<-c(beta0,beta1,sigma,D11,D12,D22,gamma1,gamma2,alpha1,alpha2)
par<-c(5,2,2,1,0.5,1,1,0.5,0.1,0.08)


rmse3<-c()
mse3<-matrix(0,nrow=rep,ncol=10)
y<-dim(JOINT2)[1]
for(j in 1:10){
  for(i in 1:y){
    mse3[i,j]<-(joint[i,j]-par[j])^2
  }
  rmse3[j]<-round(sqrt((1/y)*sum( mse3[,j])),3)
}
rmse3
##########

SE3<-c()
V3<-c()
V3<-c(var(joint[,1]),var(joint[,2]),var(joint[,3]),var(joint[,4]),var(joint[,5]),
var(joint[,6]),var(joint[,7]),var(joint[,8]),var(joint[,9]),var(joint[,10]))
SE3<-round(sqrt(V3),3)
SE3
###########
bias3<-c(est3[1]-par[1],est3[2]-par[2],est3[3]-par[3],est3[4]-par[4],
est3[5]-par[5],est3[6]-par[6],est3[7]-par[7],est3[8]-par[8],est3[9]-par[9],
est3[10]-par[10]
)
bias3<-round(bias3,3)

rbias3<-c((est3[1]-par[1])/par[1],(est3[2]-par[2])/par[2],(est3[3]-par[3])/par[3],(est3[4]-par[4])/par[4],
(est3[5]-par[5])/par[5],(est3[6]-par[6])/par[6],(est3[7]-par[7])/par[7],(est3[8]-par[8])/par[8],
(est3[9]-par[9])/par[9],(est3[10]-par[10])/par[10])
rbias3<-round(rbias3,3)
############
al<-0.05
cp3<-c()
m3<-c()
tmat3<-t(joint)
low3=up3<-matrix(0,nrow=10,ncol=y)
for(j in 1:10){
  m3[j]<-0
  for(i in 1:y){
    low3[j,i]<-tmat3[j,i]-SE3[j]*qnorm(1-al/2)
    up3[j,i]<-tmat3[j,i]+SE3[j]*qnorm(1-al/2)
    if ((low3[j,i]<par[j])&(par[j]< up3[j,i])) {m3[j]=m3[j]+1} else {m3[j]=m3[j]}
  }
  cp3[j]<-m3[j]/y
}
cp3<-round(cp3,3)

table3<-data.frame(par,est3,SE3,bias3,rbias3,rmse3,cp3)
table3
edit(table3)

#####################

