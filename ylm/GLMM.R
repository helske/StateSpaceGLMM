#install.packages("mvabund") #spider-aineisto
library(mvabund)
#install.packages("lme4")    #glmer-funktio yleistetyille lineaarisille sekamalleille
library(lme4)

data(spider)

# Muokataan aineisto glmer-funktiolle sopivaksi
Y<-as.matrix(spider$abund)
X<-spider$x
X<-X[,c(1,4,5,6)]
X<-rbind(X,X,X,X,X,X,X,X,X,X,X,X)
site<-rep(seq(1,28),12)
dataspider<-data.frame(c(Y),X,site)
names(dataspider)<-c("Y","soil.dry", "moss", "herb.layer", "reflection", "site")

# sovitetaan malli
fit.lme4 <- glmer(Y ~ soil.dry + moss + herb.layer + reflection + (1|site), family=poisson(link = log), 
                  data=dataspider,control=glmerControl(optimizer="bobyqa"))

summary(fit.lme4)
fit.lme4@beta
ranef(fit.lme4)$site           

## 
library(KFAS)

out<-glmm(response.var="Y",grouping.var="site",fixed=~soil.dry+moss+herb.layer+reflection,
          random=~1,data=dataspider,distribution="poisson",control=list(iprint=0))


coef(out,1,1)
out$model$
fit.lme4

kfas<-function() out<-glmm(response.var="Y",grouping.var="site",fixed=~soil.dry+moss+herb.layer+reflection,
                           random=~1,data=dataspider,distribution="poisson")
lme4 <-function() fit.lme4 <- glmer(Y ~ soil.dry + moss + herb.layer + reflection + (1|site), family=poisson(link = log), 
                                 data=dataspider,control=glmerControl(optimizer="bobyqa"))

library(microbenchmark)

microbenchmark(kfas(),lme4())

# muokataan aineisto KFAS-paketille sopivaan muotoon
Ymv<-as.matrix(as.data.frame(split(dataspider[,1],dataspider$site)))
Xmv<-split(dataspider[,-c(1,6)],dataspider$site)

# rakennetaan ensin malli käyttäen edellä saatuja estimaatteja:

model<-SSModel(Ymv ~ -1 
               +  SSMregression(rep(list(~soil.dry+moss+herb.layer+reflection),28), 
                               type="common", data=Xmv, a1=fit.lme4@beta,P1inf=0) 
               + SSMregression(~1, P1=diag(fit.lme4@theta^2,28)), distribution="poisson")

out<-KFS(model)
coef(out,1,1)
c(ranef(fit.lme4)$site)
fit.lme4@beta
# Nyt sama malli itse estimoiden

out<-glmm(Y=Ymv,fixed=~soil.dry+moss+herb.layer+reflection,random=~1,data=Xmv,distribution="poisson")

model<-SSModel(Ymv ~ -1
               + SSMregression(rep(list(~soil.dry+moss+herb.layer+reflection),28),type="common",data=Xmv) 
               + SSMregression(~1,P1=diag(NA,28)),distribution="poisson")


likfn<-function(pars,model,estimate=TRUE,nsim=0){
  diag(model$P1)[6:33]<-exp(pars)
  if(estimate){
    -logLik(model,nsim=nsim)  
  }else model
}

fit.kfas<-nlm(f=likfn,p=-1,print.level=1,model=model,nsim=0)

model<-likfn(fit.kfas$e,model,FALSE)
out<-KFS(model,nsim=0)

summary(fit.glmer)
coef(out,1,1)

cov2cor(out$V[1:5,1:5,1])


sim<-importanceSSM(model,type="signal",nsim=1000,antithetics=FALSE)

#install.packages("plotrix")
library(plotrix)
weighted.hist(exp(sim$s[10,28,]),sim$w)

f1<-function(init){
  model<-SSModel(Ymv ~ -1
                 + SSMregression(rep(list(~soil.dry+moss+herb.layer+reflection),28),type="common",data=Xmv) 
                 + SSMregression(~1,P1=diag(NA,28)),distribution="poisson")
  
  
  likfn<-function(pars,model,estimate=TRUE,nsim=0){
    diag(model$P1)[6:33]<-exp(pars)
    if(estimate){
      -logLik(model,nsim=nsim)  
    }else model
  }
  
  fit.kfas<-nlm(f=likfn,p=init,print.level=0,model=model,nsim=0)
  
  model<-likfn(fit.kfas$e,model,FALSE)
  out<-KFS(model,nsim=0)
}

f2<-function() fit.lme4 <- glmer(Y ~ soil.dry + moss + herb.layer + reflection + (1|site), family=poisson(link = log), 
                  data=dataspider,control=glmerControl(optimizer="bobyqa"))

library(microbenchmark)

microbenchmark(f1(-1),f2())
#Unit: seconds
#expr      min       lq   median       uq      max neval
#f1() 1.296971 1.305978 1.313591 1.351389 1.557538   100
#f2() 2.112143 2.139694 2.148727 2.159271 2.438743   100

microbenchmark(f1(2),f2())
#Unit: seconds
#expr      min       lq   median       uq      max neval
#f1(2) 2.519638 2.527309 2.535835 2.579823 3.037389   100
#f2() 2.104172 2.115046 2.124332 2.137799 2.409674   100
