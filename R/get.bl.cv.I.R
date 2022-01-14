#' @import stats
#' @import NPrior
#' @import nlshrink
get.bl.cv.I <- function(l,path, chr0, MCMC,BURN,ncv1,shrink=F){
  # library(EBPRS)
  if(file.exists(paste0(path,'/ldblocks/chr',chr0,'/LD_region_',l,'.fam'))){
    LD <- read_plink(paste0(path,'/ldblocks/chr',chr0,'/LD_region_',l))
    # LD$bim$order.ld <- 1:nrow(LD$bim)
    # colnames(LD$bim)[c(2,5,6)] <- c('SNP','A1','A2')
    #
    # dat <- merge(summs,LD$bim,by='SNP')
    xx <- scale(bedNA(LD$bed))
    yy <- scale(LD$fam$V6)
    #  print('model start')
    p <- ncol(xx)
    nn <- nrow(xx)
    set.seed(0)
    cvind <- sample(ncv1,nn,replace = T)
    rr44 <- rr22 <- rr99 <- ntest <- list()
    sse44 <- sse22 <- sse99 <- mse44 <- mse22 <- mse99 <- list()
    yy.v <- yy.s2 <- yy.s4 <- yy.s9 <- list()
    for(cv1 in 1:ncv1){
      fit22 <-  NPrior_run(xx[cvind!=cv1,],yy[cvind!=cv1],prior = "SpSL-C", a0 = 1, b0 = 1, prior_sig_type = 1,N=MCMC,BURN=BURN, verbose=F)
      fit33 <-  NPrior_run(xx[cvind!=cv1,],yy[cvind!=cv1], prior = "HS", eta = 1,N=MCMC,BURN=BURN,verbose=F)
      fit44 <-  NPrior_run(xx[cvind!=cv1,],yy[cvind!=cv1], prior = "SpSL-L",  prior_sig_type = 1,N=MCMC,BURN=BURN,verbose=F)
      ntest[[cv1]] <- length(which(cvind==cv1))
      rr44[[cv1]] <- cor(yy[cvind==cv1],xx[cvind==cv1,]%*%apply(fit22$THETA,1,mean))
      rr22[[cv1]] <- cor(yy[cvind==cv1],xx[cvind==cv1,]%*%apply(fit33$THETA,1,mean))
      rr99[[cv1]] <- cor(yy[cvind==cv1],xx[cvind==cv1,]%*%apply(fit44$THETA,1,mean))
      mse44[[cv1]] <- mean((scale(yy[cvind==cv1])-scale(xx[cvind==cv1,]%*%apply(fit22$THETA,1,mean)))^2)
      mse22[[cv1]] <- mean((scale(yy[cvind==cv1])-scale(xx[cvind==cv1,]%*%apply(fit33$THETA,1,mean)))^2)
      mse99[[cv1]] <- mean((scale(yy[cvind==cv1])-scale(xx[cvind==cv1,]%*%apply(fit44$THETA,1,mean)))^2)
      sse44[[cv1]] <- sum((scale(yy[cvind==cv1])-scale(xx[cvind==cv1,]%*%apply(fit22$THETA,1,mean)))^2)
      sse22[[cv1]] <- sum((scale(yy[cvind==cv1])-scale(xx[cvind==cv1,]%*%apply(fit33$THETA,1,mean)))^2)
      sse99[[cv1]] <- sum((scale(yy[cvind==cv1])-scale(xx[cvind==cv1,]%*%apply(fit44$THETA,1,mean)))^2)
      yy.v[[cv1]] <- yy[cvind==cv1]
      yy.s2[[cv1]] <- xx[cvind==cv1,]%*%fit33$ThetaHat
      yy.s4[[cv1]] <- xx[cvind==cv1,]%*%fit22$ThetaHat
      yy.s9[[cv1]] <- xx[cvind==cv1,]%*%fit44$ThetaHat
    }
    return(list(ntest=ntest,rr44=rr44,rr22=rr22,rr99=rr99,
                mse22=mse22,mse44=mse44,mse99=mse99,sse22=sse22,sse44=sse44,sse99=sse99,
                yy.v=yy.v,yy.s4=yy.s4, yy.s2=yy.s2, yy.s9=yy.s9))
  }}
