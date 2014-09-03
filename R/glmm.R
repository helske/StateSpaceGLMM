#' Fit Generalized Linear Mixed-Effects Models Using State Space Framework
#'
#' Function \code{glmm} estimates GLMM using methods basec on state space modelling. 
#'
#' @export
#' @import nloptr
#' @import KFAS
#' @param response.var Name of the response variable in data.
#' @param grouping.var Name of the grouping variable in data. Only one grouping variable is allowed.
#' @param data Optional data frame environment containing the variables in the model.
#' @param distribution Character vector defining the distributions used for each group. 
#' Either length of one (same distribution for all groups) or length of p, number of groups. Default is "gaussian".
#' @param method "REML", "REML-ML" or "ML". "REML-ML" uses estimates obtained from REML for ML estimation.
#' Note that ML is much slower.
#' @param init.random Initial values for random effect covariances.
#' @param init.dispersion Initial values for dispersion paremeters.
#' @param init.fixed Initial values for fixed effects.
#' @param correlating.effects Logical. Default is TRUE.
#' @param estimate.dispersion Logical. Is the dispersion parameter estimated or fixed.
#' @param return.model Logical. Is the estimated model returned in output.
#' @param nsim Integer. Number of independent samples used in importance sampling. 
#' Default is 0, which corresponds to Laplace approximation.
#' @param maxiter Integer. Number of iterations for in iterative weighted least squares.
#' 
#' @examples
#' data(butterfly)
#' 
#' require(lme4)
#' system.time(fit.SSglmm.reml<-glmm(response.var="Colias", grouping.var="site", fixed=~habitat + building + urbanveg,random=~1,
#'                                  REML=TRUE,data=butterfly,distribution="poisson",return.model=TRUE))
#'
#' system.time(fit.SSglmm.ml<-glmm(response.var="Colias", grouping.var="site", fixed=~habitat + building + urbanveg,random=~1,
#'                                REML=FALSE,data=butterfly,distribution="poisson",return.model=TRUE))
#'
#' fit.glmer<-glmer(Colias ~ habitat + building + urbanveg + (1|site), 
#'                 family=poisson, data=butterfly,control=glmerControl(optimizer = "bobyqa"))
#'
#' sqrt(fit.SSglmm.reml$random$P)
#' sqrt(fit.SSglmm.ml$random$P)
#' VarCorr(fit.glmer)
#'
#' fixef(fit.glmer)
#' fit.SSglmm.ml$fixed$coef
#' fit.SSglmm.reml$fixed$coef
#'
#' require(mvabund)
#' data(spider)
#' x<-apply(spider$x,2,rep,times=12)
#' dataspider<-data.frame(species=rep(names(spider$abund),each=28),abund=unlist(spider$abund),x)
#' out<-glmm(response.var="abund", grouping.var="species", fixed=~soil.dry+bare.sand, random=~1, 
#' data=dataspider,distribution="poisson")
#' 
#' out2<-glmm(response.var="abund", grouping.var="species", fixed=~soil.dry+bare.sand,
#' data=dataspider,distribution="poisson")
#' fit<-manyglm(mvabund(spider$abund)~soil.dry+bare.sand,data=data.frame(spider$x),family="poisson")
#' 
#' 

glmm<-function(response.var, grouping.var, fixed.formula, random.formula, data, distribution,
               init.random, init.dispersion, init.fixed, correlating.effects=TRUE, REML=TRUE,
               estimate.dispersion, return.model=TRUE,nsim=0, maxiter=50,maxeval=1000,xtol_rel=1e-6,...){
  
  
  #modify data
  Y<-split(data[response.var],data[grouping.var])
  Y<-matrix(unlist(Y,use.names=FALSE),ncol=length(Y),dimnames=list(NULL,paste0(names(Y[[1]]),".",names(Y))))
  data<-split(data,data[grouping.var])
  
  
  p<-ncol(Y)
  n<-nrow(Y)
  
  
  if(missing(init.dispersion)){
    init.dispersion<-rep(1,p)
  } else {
    if(length(init.dispersion)==1){
      init.dispersion<-rep(init.dispersion,p)
    } else if(length(init.dispersion)!=p) stop("Number of initial values for dispersion parameters is not equal to number of groups.")
  }
  
  if(missing(distribution))
    distribution<-rep("gaussian",p)
  if(all(distribution=="gaussian")){
    gaussian<-TRUE
  }else gaussian<-FALSE
  
  if(missing(estimate.dispersion))
    estimate.dispersion<-any(distribution%in%c("negative binomial","gamma"))
  
  
  # Find correct number of parameters
  
  if(missing(random.formula)){
    k.fix<-p*ncol(model.matrix(fixed.formula,data=data[[1]]))
    k.rand<-0
    correlating.effects<-FALSE  
    if(gaussian){
      model<-SSModel(Y~-1
                     + SSMregression(rep(list(fixed.formula),p),data=data), H=diag(p))
    } else{
      model<-SSModel(Y~-1
                     + SSMregression(rep(list(fixed.formula),p),data=data), distribution=distribution,u=init.dispersion)
    } 
    
  } else {
    
    k.fix<-ncol(model.matrix(fixed.formula,data=data[[1]]))
    k.rand<-ncol(model.matrix(random.formula,data=data[[1]]))
    
    if(k.rand<2) correlating.effects<-FALSE
    
    if(missing(init.random)){
      if(correlating.effects){ 
        init.random <- diag(k.rand)[upper.tri(diag(k.rand),TRUE)] 
      }else{ 
        init.random <- rep(1, k.rand)
      }
    }else {
      if(length(init.random)==1){
        init.random<-rep(init.random,if(correlating.effects) k.rand*(k.rand+1)/2 else k.rand)
      } else if(length(init.random)!=(if(correlating.effects) k.rand*(k.rand+1)/2 else k.rand))
        stop("Incorrect number of initial values for random effect covariances.")
      
    }
    
    if(gaussian){
      model<-SSModel(Y~-1
                     + SSMregression(rep(list(fixed.formula),p),type="common",data=data)
                     + SSMregression(rep(list(random.formula),p),P1=diag(NA,p*k.rand),data=data),
                     H=diag(p))
    } else{
      model<-SSModel(Y~-1
                     + SSMregression(rep(list(fixed.formula),p),type="common",data=data)
                     + SSMregression(rep(list(random.formula),p),P1=diag(NA,p*k.rand),data=data),
                     distribution=distribution,u=init.dispersion)
    }
    
  }
    
  
  if(k.rand>0 || estimate.dispersion || gaussian){ #something to estimate other than fixed effects
    
    if(REML){
      if(correlating.effects){
        likfn<-function(pars,model,estimate=TRUE,nsim=0,maxiter=50){
          if(estimate.dispersion)
            model$u[]<-matrix(pars[1:p],n,p,byrow=TRUE)          
          if(k.rand>0){
            P1<-matrix(0,k.rand,k.rand)
            P1[upper.tri(P1,TRUE)]<-pars[estimate.dispersion*p+1:length(init.random)]
            P1<-crossprod(P1)
            model$P1[(k.fix+1):(k.fix+p*k.rand),(k.fix+1):(k.fix+p*k.rand)]<-
              as.matrix(.bdiag(replicate(p,P1,simplify=FALSE)))
          }
          if(gaussian)
            diag(model$H[,,1])<-pars[length(pars)]
          if(estimate){
            -logLik(model,nsim=nsim,maxiter=maxiter)  
          }else model
        }
        inits<-c(if(estimate.dispersion) init.dispersion else NULL, init.random,if(gaussian) 1 else NULL)
        tmp<-diag(k.rand)
        lower<-c(rep(0,(estimate.dispersion)*p),
                 if(correlating.effects) 
                   ifelse(tmp[upper.tri(tmp,TRUE)]==0,-Inf,0) else rep(0,k.rand),if(gaussian) 0 else NULL)
        
        
        fit<-nloptr(eval_f=likfn,x0=inits,lb=lower,model=model,estimate=TRUE,
                                 nsim=nsim,maxiter=maxiter,opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=maxeval,xtol_rel=xtol_rel,...))      
        
      } else {
        likfn<-function(pars,model,estimate=TRUE,nsim=0,maxiter=50){
          if(estimate.dispersion)
            model$u[]<-matrix(pars[1:p],n,p,byrow=TRUE)
          if(k.rand>0)
            diag(model$P1)[(k.fix+1):(length(model$a1))]<-pars[estimate.dispersion*p+1:k.rand]  
          if(gaussian)
            diag(model$H[,,1])<-pars[length(pars)]
          if(estimate){
            -logLik(model,nsim=nsim,maxiter=maxiter)  
          }else model
        }
        inits<-c(if(estimate.dispersion) init.dispersion else NULL, init.random,if(gaussian) 1 else NULL)
        
        
        fit<-nloptr(eval_f=likfn,x0=inits,lb=rep(0,length(inits)),model=model,estimate=TRUE,
                                 nsim=nsim,maxiter=maxiter,opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=maxeval,xtol_rel=xtol_rel,...))       
        
      }
    } else {
      if(missing(init.fixed)) init.fixed<-0
      if(length(init.fixed)==1){
        init.fixed<-rep(init.fixed,k.fix)
      } else if(length(init.fixed)!=k.fix)
        stop("Incorrect number of initial values for fixed effects.")
      
      model$P1inf[]<-0
      
      if(correlating.effects){
        
        likfn<-function(pars,model,estimate=TRUE,nsim=0,maxiter=50){
          if(estimate.dispersion)
            model$u[]<-matrix(pars[1:p],n,p,byrow=TRUE)
          
          if(k.rand>0){
            P1<-matrix(0,k.rand,k.rand)
            P1[upper.tri(P1,TRUE)]<-pars[estimate.dispersion*p+1:length(init.random)]
            P1<-crossprod(P1)
            model$P1[(k.fix+1):(k.fix+p*k.rand),(k.fix+1):(k.fix+p*k.rand)]<-
              as.matrix(.bdiag(replicate(p,P1,simplify=FALSE)))
          }
          model$a1[1:k.fix]<-pars[estimate.dispersion*p+length(init.random)+1:k.fix]
          if(gaussian)
            diag(model$H[,,1])<-pars[length(pars)]
          if(estimate){
            -logLik(model,nsim=nsim,maxiter=maxiter)  
          }else model
        }
        inits<-c(if(estimate.dispersion) init.dispersion else NULL, init.random,
                 init.fixed,if(gaussian) 1 else NULL)
        tmp<-diag(k.rand)
        lower<-c(rep(0,(estimate.dispersion)*p),
                 if(correlating.effects) 
                   ifelse(tmp[upper.tri(tmp,TRUE)]==0,-Inf,0) else rep(0,k.rand),rep(-Inf,k.fix),
                 if(gaussian) 0 else NULL)
        
        
        fit<-nloptr(eval_f=likfn,x0=inits,lb=lower,model=model,estimate=TRUE,
                                 nsim=nsim,maxiter=maxiter,opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=maxeval,xtol_rel=xtol_rel,...))       
        
        
      } else {
        likfn<-function(pars,model,estimate=TRUE,nsim=0,maxiter=50){
          if(estimate.dispersion)
            model$u[]<-matrix(pars[1:p],n,p,byrow=TRUE)
          if(k.rand>0)
            diag(model$P1)[(k.fix+1):(length(model$a1))]<-pars[estimate.dispersion*p+1:k.rand]  
          model$a1[1:k.fix]<-pars[estimate.dispersion*p+k.rand+1:k.fix]
          if(gaussian)
            diag(model$H[,,1])<-pars[length(pars)]
          if(estimate){
            -logLik(model,nsim=nsim,maxiter=maxiter)  
          }else model
        }
        inits<-c(if(estimate.dispersion) init.dispersion else NULL, init.random,init.fixed, if(gaussian) 1 else NULL)
        
        
        fit<-nloptr(eval_f=likfn,x0=inits,lb=c(rep(0,length(inits)-k.fix-gaussian),rep(-Inf,k.fix),if(gaussian) 0 else NULL),
                                 model=model,estimate=TRUE,
                                 nsim=nsim,maxiter=maxiter,opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=maxeval,xtol_rel=xtol_rel,...))
        
        
        
      }
    } 
    
    
  } else if(!REML){
    
    if(missing(init.fixed)) init.fixed<-0
    if(length(init.fixed)==1){
      init.fixed<-rep(init.fixed,k.fix)
    } else if(length(init.fixed)!=k.fix)
      stop("Incorrect number of initial values for fixed effects.")
    
    model$P1inf[]<-0      
    
    likfn<-function(pars,model,estimate=TRUE,nsim=0,maxiter=50){         
      model$a1[1:k.fix]<-pars[1:k.fix]
      if(gaussian)
        diag(model$H[,,1])<-pars[length(pars)]
      if(estimate){
        -logLik(model,nsim=nsim,maxiter=maxiter)  
      }else model
    }
    inits<-c(init.fixed,if(gaussian) 1 else NULL)
    
    
    fit<-nloptr(eval_f=likfn,x0=inits,lb=c(rep(-Inf,k.fix),if(gaussian) 0 else NULL),
                             model=model,estimate=TRUE,
                             nsim=nsim,maxiter=maxiter,opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=maxeval,xtol_rel=xtol_rel,...))
    
    
    
  }
  
  model<-likfn(fit$solution,model,FALSE)   
  
  if(fit$status!=1)
    warning(fit$message)  
  
  
  kfs.out<-KFS(model,smoothing=c("state","mean"),nsim=nsim,maxiter=maxiter)  
  
 
  ll<-logLik(model,nsim=nsim,maxiter=maxiter)
  
  results<-list(fixed=list(coef=coef(kfs.out,1,1)[1:k.fix],V_fixed=kfs.out$V[1:k.fix,1:k.fix,1]),
                random=if(k.rand>0) list(effects=coef(kfs.out,1,1)[-(1:k.fix)],V_random=kfs.out$V[-(1:k.fix),-(1:k.fix),1],
                                         P=model$P1[(k.fix+1):(k.fix+k.rand),(k.fix+1):(k.fix+k.rand)]) else NULL,
                fitted=list(fitted=kfs.out$mu,V=kfs.out$V_muhat),
                dispersions=if(!gaussian) model$u[1,] else model$H[1,1,1],logLik=ll,
                call= match.call(expand.dots = FALSE), model=if(return.model) model else NULL)
  class(results)<-"glmm.results"
  results
}