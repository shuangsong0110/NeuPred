#' @import stats
#' @import NPrior
#' @import nlshrink
get.bl <- function(path, summs,chr0, l, MCMC,BURN,n,shrink,prior,scale=T,eta0=1){
 # library(EBPRS)
  LD <- read_plink(paste0(path,'/ldblocks/chr',chr0,'/LD_region_',l))
  LD$bim$order.ld <- 1:nrow(LD$bim)
  colnames(LD$bim)[c(2,5,6)] <- c('SNP','A1','A2')

  dat <- merge(summs,LD$bim,by='SNP')
  #if(!'beta'%in%colnames(dat)){dat$beta <- log(dat$or)}
  dat <- dat[!is.na(dat$BETA),]
  sig <- agtc(dat$A1.x,dat$A2.x,dat$A1.y,dat$A2.y)
  dat$BETA <- dat$BETA*sig
  dat <- dat[sig!=0,]
  dat <- dat[order(dat$order.ld),]
  bb <- -qnorm(dat$P/2)*sign(dat$BETA)/sqrt(n)
  bb[bb>35/sqrt(n)]=35/sqrt(n)
  bb[bb< -35/sqrt(n)]= -35/sqrt(n)
  bb.sd <- dat$BETA/bb
  #bb=bb[dat$order.ld]
  X2.temp <- bedNA(LD$bed[,dat$order.ld])
 # library(nlshrink)
  if(shrink){
    R <- try(cov2cor(nlshrink_cov(X2.temp)))
    if(is.null(dim(R))){
      print('TRY ERROR')
      R <- cor(X2.temp)
      R[which(abs(R)<2/sqrt(nrow(LD$fam)))] <- 0
    }
    R[which(is.na(R))] <- 0
    diag(R) <- 1
    temp <- eigen(R,symmetric = T)
    temp$values[temp$values<0.1] <- 0
  }else{
    R <- cor(X2.temp)
    R[which(is.na(R))] <- 0
    diag(R) <- 1
    temp <- eigen(R,symmetric = T)
    temp$values[temp$values<1e-6] <- 0
  }
  eff.num <- min(length(which(temp$values>0)),nrow(LD$fam))
  uu <- (t(temp$vectors))[1:eff.num,]
  vv <- diag(sqrt(1/(temp$values[1:eff.num])))
  if(n>nrow(summs)){
    K <- sqrt(nrow(summs))*vv%*%uu
  }else{
    K <- sqrt(n)*vv%*%uu
  }
  xx <- as.matrix(K%*%R)
  yy <- as.matrix(K%*%bb)
  # print('model start')
  p <- ncol(xx)
  if(prior==4){
    fit11 <-  NPrior_run(xx,yy, prior = "SpSL-C", a0 = 1, b0 = 1, prior_sig_type = 1,N=MCMC,BURN=BURN, verbose=F,eta=eta0)
  }
  if(prior==2){
    fit11 <-  NPrior_run(xx,yy,prior = "HS", eta = 1,N=MCMC,BURN=BURN, verbose=F)
  }
  if(prior==9){
    fit11 <-  NPrior_run(xx,yy, prior = "SpSL-L",  prior_sig_type = 1,N=MCMC,BURN=BURN, verbose=F,eta=eta0)
  }
  if(scale){
  dat$THETA <- fit11$ThetaHat
  }else{
    dat$THETA <- fit11$ThetaHat * bb.sd
  }
  return(dat)
}
