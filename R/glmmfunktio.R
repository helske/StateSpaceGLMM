glmm<-function(response.var, grouping.var, fixed.formula, random.formula, data, distribution, REML=TRUE, 
               init.random=NULL, init.dispersion=NULL, init.fixed=NULL, correlating.effects=TRUE,
               estimate.dispersion, return.model=FALSE,nsim=0, maxiter=50,...){
  
  
  Y<-split(data[response.var],data[grouping.var])
  Y<-matrix(unlist(Y,use.names=FALSE),ncol=length(Y),dimnames=list(NULL,paste0(names(Y[[1]]),".",names(Y))))
  data<-split(data,data[grouping.var])
  
  if(missing(estimate.dispersion))
    estimate.dispersion<-any(distribution%in%c("negative binomial","gamma"))
  
  p<-ncol(Y)
  n<-nrow(Y)
  if(is.null(init.dispersion))
    init.dispersion<-rep(1,p)
  
  
  if(missing(random.formula)){
    k.fix<-p*ncol(model.matrix(fixed.formula,data=data[[1]]))
    k.rand<-0    
    model<-SSModel(Y~-1
                   + SSMregression(rep(list(fixed.formula),p),data=data),
                   distribution=distribution,u=init.dispersion)
    
  } else {
    
    k.fix<-ncol(model.matrix(fixed.formula,data=data[[1]]))
    k.rand<-ncol(model.matrix(random.formula,data=data[[1]]))
    if(is.null(init.random))
      init.random<-rep(1,if(correlating.effects) k.rand*(k.rand+1)/2 else k.rand)
    
    model<-SSModel(Y~-1
                   + SSMregression(rep(list(fixed.formula),p),type="common",data=data)
                   + SSMregression(rep(list(random.formula),p),P1=diag(NA,p*k.rand),data=data),
                   distribution=distribution,u=init.dispersion)
  }
  
  
  if(k.rand>0 || estimate.dispersion || !REML){
    
    if(is.null(init.fixed) && !REML)
      init.fixed<-rep(1,k.fix)
    
    if(k.rand<2) correlating.effects<-FALSE
    
    if(REML){
      if(correlating.effects){
        likfn<-function(pars,model,estimate=TRUE,nsim=0,maxiter=50){
          if(estimate.dispersion)
            model$u[]<-matrix(pars[1:p],n,p,byrow=TRUE)
          if(k.rand>0)
            P1<-matrix(0,k.rand,k.rand)
          P1[upper.tri(P1,TRUE)]<-pars[estimate.dispersion*p+1:length(init.random)]
          P1<-crossprod(P1)
          model$P1[]<-as.matrix(.bdiag(replicate(p,P1,simplify=FALSE)))
          if(estimate){
            -logLik(model,nsim=nsim,maxiter=maxiter)  
          }else model
        }
        inits<-c(if(estimate.dispersion) init.dispersion else NULL, init.random)
        tmp<-diag(k.rand)
        lower<-c(rep(0,(estimate.dispersion)*p),
                 if(correlating.effects) 
                   ifelse(tmp[upper.tri(tmp,TRUE)]==0,-Inf,0) else rep(0,k.rand))
        
        fit<-bobyqa(f=likfn,p=inits,lower=lower,model=model,
                    nsim=nsim,maxiter=maxiter,...)
      } else {
        likfn<-function(pars,model,estimate=TRUE,nsim=0,maxiter=50){
          if(estimate.dispersion)
            model$u[]<-matrix(pars[1:p],n,p,byrow=TRUE)
          if(k.rand>0)
            diag(model$P1)[(k.fix+1):(length(model$a1))]<-pars[estimate.dispersion*p+1:k.rand]  
          if(estimate){
            -logLik(model,nsim=nsim,maxiter=maxiter)  
          }else model
        }
        inits<-c(if(estimate.dispersion) init.dispersion else NULL, init.random)
        
        fit<-bobyqa(f=likfn,p=inits,lower=rep(0,inits),model=model,nsim=nsim,maxiter=maxiter,...)
      }
    } else {     
      model$P1inf[]<-0
      if(correlating.effects){
        
        likfn<-function(pars,model,estimate=TRUE,nsim=0,maxiter=50){
          if(estimate.dispersion)
            model$u[]<-matrix(pars[1:p],n,p,byrow=TRUE)
          if(k.rand>0)
            P1<-matrix(0,k.rand,k.rand)
          P1[upper.tri(P1,TRUE)]<-pars[estimate.dispersion*p+1:length(init.random)]
          P1<-crossprod(P1)
          model$P1[]<-as.matrix(.bdiag(replicate(p,P1,simplify=FALSE)))
          model$a1[1:k.fix]<-pars[estimate.dispersion*p+length(init.random)+1:k.fix]
          if(estimate){
            -logLik(model,nsim=nsim,maxiter=maxiter)  
          }else model
        }
        inits<-c(if(estimate.dispersion) init.dispersion else NULL, init.random,init.fixed)
        tmp<-diag(k.rand)
        lower<-c(rep(0,(estimate.dispersion)*p),
                 if(correlating.effects) 
                   ifelse(tmp[upper.tri(tmp,TRUE)]==0,-Inf,0) else rep(0,k.rand),rep(-Inf,k.fix))
        
        fit<-bobyqa(f=likfn,p=inits,lower=lower,model=model,
                    nsim=nsim,maxiter=maxiter,...)
      } else {
        likfn<-function(pars,model,estimate=TRUE,nsim=0,maxiter=50){
          if(estimate.dispersion)
            model$u[]<-matrix(pars[1:p],n,p,byrow=TRUE)
          if(k.rand>0)
            diag(model$P1)[(k.fix+1):(length(model$a1))]<-pars[estimate.dispersion*p+1:k.rand]  
          model$a1[1:k.fix]<-pars[estimate.dispersion*p+k.rand+1:k.fix]
          if(estimate){
            -logLik(model,nsim=nsim,maxiter=maxiter)  
          }else model
        }
        inits<-c(if(estimate.dispersion) init.dispersion else NULL, init.random,init.fixed)
        
        fit<-bobyqa(f=likfn,p=inits,lower=c(rep(0,length(inits)-k.fix),rep(-Inf,k.fix)),model=model,
                    nsim=nsim,maxiter=maxiter,...)
      }
    }
    model<-likfn(fit$p,model,FALSE)
    
    
    
  }
  out<-KFS(model,smoothing=c("state","mean"),nsim=nsim,maxiter=maxiter)
  
  list(fixed=list(coef=coef(out,1,1)[1:k.fix],V_fixed=out$V[1:k.fix,1:k.fix,1]),
       random=if(k.rand>0) list(effects=coef(out,1,1)[-(1:k.fix)],V_random=out$V[-(1:k.fix),-(1:k.fix),1],
                                P=model$P1[(k.fix+1):(k.fix+k.rand),(k.fix+1):(k.fix+k.rand)]) else NULL,
       fitted=list(fitted=out$mu,V=out$V_muhat),
       dispersions=model$u[1,],logLik=-fit$fval,
       call= match.call(expand.dots = FALSE), model=if(return.model) model else NULL)
  
  
}