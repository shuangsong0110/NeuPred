# Infer the z value's distribution using EM with finite mixture model
# Input: z value-z[vector], number of non-null components-K[scalar,4]
#        Maximum iteration number-maxIter[scalar,1000]
#        Tolerance-tol[scalar, 1e-3], Show iteration info-info[bool, T]
# Output: The proportion of null components-pi0[scalar]
#         The estimated variance of null component-sigma02[scalar, 0 if empNull==F]
#         The proportion of non-null components-Pi1[vector]
#         The estimated variance of non-null components-sigma2[vector]
#         The probability of being null component-lfdr[vector]
#		  The probability of being each of the non-null components-ltdr[list]
#         Iteration number-iter[scalar]
#         Expectation of composite (Marginal) log likelihood-Qval[scalar]
#         Log likelihood-logLik[scalar]
snpEM<-function(z, K=3, maxIter=1000, tol=1e-3, beta0=length(z)/20, empNull=FALSE, info=TRUE){
  #For real data usage, a certain penality to pi0, such as beta0=length(z)/5 is needed.
  #For simulated data, this penalization will increase bias.
  m<-length(z)
  #initialization
  pi0_0<-0.95
  Pi1_0<-rep((1-pi0_0)/K,K)
  sigma2_0<-rgamma(K,1)
  sigma02_0<-0

  pi0_t<-pi0_0
  Pi1_t<-Pi1_0
  sigma2_t<-sigma2_0
  sigma02_t<-sigma02_0
  h<-list()
  h0<-0
  tmpH0<-0
  tmpH<-list()
  nanVal<-function(x) ifelse(is.nan(x),0,x)
  for(iter in 1:maxIter){
    #E step
    for(i in 0:K){
      if(i==0) tmpH0<-pi0_t*dnorm(z, sd=sqrt(1+sigma02_t))
      else tmpH[[i]]<-Pi1_t[i]*dnorm(z,sd=sqrt(1+sigma2_t[i]))
    }
    norH<-tmpH0+Reduce('+',tmpH)
    h0<-nanVal(tmpH0/norH)
    h<-lapply(tmpH,FUN=function(x) return(nanVal(x/norH)))

    #M step
    pi0_t1<-(sum(h0)+beta0)/(m+beta0)
    Pi1_t1<-sapply(h,FUN=sum)/(m+beta0)
    sigma2_t1<-sapply(h, FUN=function(x){
      if(sum(x)==0) return(0)
      else return(max(sum(x*z^2)/sum(x)-1,0))
    } )
    if(empNull==TRUE){
      sigma02_t1<-max(sum(h0*z^2)/sum(h0)-1,0)
    }else{
      sigma02_t1<-sigma02_t
    }

    if( (abs(nanVal((pi0_t1-pi0_t)/pi0_t))<tol) && sqrt(nanVal(sum((Pi1_t1-Pi1_t)^2)/sum(Pi1_t^2)))<tol &&
        Reduce('+',Map(function(x,y){return(ifelse(sum(y^2)==0, 0, sqrt(sum((x-y)^2)/sum(y^2))))},sigma2_t1,sigma2_t))<tol && nanVal((sigma02_t1-sigma02_t)/sigma02_t)<tol) break
    else{
      pi0_t<-pi0_t1
      Pi1_t<-Pi1_t1
      sigma2_t<-sigma2_t1
      sigma02_t<-sigma02_t1
    }
  }

  #	Qval<-ifelse(tmpH0==0,0,sum(h0*log(tmpH0)))+sum(Reduce('+',Map('*',h,lapply(tmpH,function(x){return(ifelse(x==0,0,log(x)))}))))
  logF<-function(x) ifelse(x==0,0,log(x))
  Qval<-sum(h0*logF(tmpH0))+sum(Reduce('+',Map('*',h,lapply(tmpH,logF))))
  #	library(MCMCpack)
  logLik <- sum(log(dnormMix(z, pi0_t, Pi1_t, sigma2_t)))#+log(ddirichlet(c(pi0_t, Pi1_t), c(beta0+1, rep(1,K))))
  if(info){
    cat('pi0:',pi0_t,'\n')
    cat('sigma02:',sigma02_t,'\n')
    cat('Pi1:\n')
    print(Pi1_t)
    cat('sigma^2:\n')
    print(sigma2_t)
    cat('Iteration:',iter,'\nLog-likelihood:',logLik,'\n')
  }
  return(list(pi0=pi0_t,sigma02=sigma02_t,Pi1=Pi1_t,sigma2=sigma2_t,lfdr=h0,ltdr=h,iter=iter,Qval=Qval, logLik=logLik))
}
normMix<-function(q, sigma02, Pi1, sigma2){
  return( (1-sum(Pi1))*pnorm(q, 0, sqrt(1+sigma02))+Reduce('+', Map(function (x, y) x*pnorm(q, 0, sqrt(1+y)), Pi1, sigma2)) )
}

rnormMix<-function(n, sigma02, Pi1, sigma2){
  K <- length(Pi1)
  indicator <- sample(0:K, n, replace=T, prob=c(1-sum(Pi1), Pi1))
  return(list(ind=indicator, z=sapply(indicator, function(x) if(x==0) rnorm(1,0,sqrt(1+sigma02)) else rnorm(1, 0, sqrt(1+sigma2[x]))) ) )
}

dnormMix <- function(x, sigma02, Pi1, sigma2){
  return( (1-sum(Pi1))*dnorm(x, 0, sqrt(1+sigma02))+Reduce('+', Map(function (xx, yy) xx*dnorm(x, 0, sqrt(1+yy)), Pi1, sigma2)) )
}

