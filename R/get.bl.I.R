#' @import stats
#' @import NPrior
get.bl.I <- function(path, chr0, l, MCMC,BURN,prior,scale=T){
  train <- read_plink(paste0(path,'/ldblocks/chr',chr0,'/LD_region_',l))
  X <- train$bed
  dat <- train$bim
 # dat1 <- dat[,c('SNP','CHR','A1.y','A2.y','P','THETA')]
  colnames(dat) <- c('CHR','SNP','morgans','pos','A1.y','A2.y')
  X1 <- bedNA(X)
  xx.sd <- apply(X1,2,sd,na.rm=T)
  xx.sd[xx.sd==0] <- sqrt(2*0.5*0.5)
  xx.sd[is.na(xx.sd)] <- sqrt(2*0.5*0.5)

  xx <- scale(bedNA(X))
  yy <- scale(train$fam$V6)
  yy[which(is.na(yy))] <- mean(yy,na.rm=T)
  if(prior==4){
    fit11 <-  NPrior_run(xx,yy, prior = "SpSL-C", a0 = 1, b0 = 1, prior_sig_type = 1,N=MCMC,BURN=BURN, verbose=F)
  }
  if(prior==2){
    fit11 <-  NPrior_run(xx,yy,prior = "HS", eta = 1,N=MCMC,BURN=BURN, verbose=F)
  }
  if(prior==9){
    fit11 <-  NPrior_run(xx,yy, prior = "SpSL-L",  prior_sig_type = 1,N=MCMC,BURN=BURN, verbose=F)
  }
  if(scale){
  dat$THETA <- fit11$ThetaHat
  }else{
    dat$THETA <- fit11$ThetaHat / xx.sd
  }
  return(dat)
}
