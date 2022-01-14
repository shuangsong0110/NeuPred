#' @import stats
#' @import data.table
#' @import nlshrink

get.prs <- function(path,chr,select.prior,cores,auto=T){
  path1 <- paste0(path,'/testfiles')
  if(auto){
    S.all <- list()
    rr <- 0
    for(chr0 in chr){
      test <- read_plink(paste0(path1,'/chr',chr0))
      for(prior1 in c(2,4,9)){
        pp <- c('HS','C','L')[which(c(2,4,9)==prior1)]
        path1 <- paste0(path,'/testfiles')
        res <- fread(paste0(path,'/posterior/chr',chr0,'_',pp,'.txt'))
        colnames(test$bim)[2] <- 'SNP'
        test$bim$order <- 1:nrow(test$bim)
        ww <- merge(test$bim,res,by='SNP')
        sig <- agtc(ww$V5,ww$V6,ww$A1.y,ww$A2.y)
        ww$THETA <- ww$THETA*sig
        X <- bedNA(test$bed[,ww$order])
        X1 <- scale(X)
        X1[which(is.na(X1))] <- 0
        if(rr==0){
          S.all[[which(c(2,4,9)==prior1)]] <- X1%*%ww$THETA
        }else{
          S.all[[which(c(2,4,9)==prior1)]] <- S.all[[which(c(2,4,9)==prior1)]]+X1%*%ww$THETA
        }
      }
      rr <- rr+1
    }
    weight <- unlist(lapply(S.all,sd))
    S <- NULL
    tt <- 1
    corr <- c()
    #select.prior <- c()
    for(chr0 in chr){
      #prior1 <- find.prior(chr0=chr0,path=path,summs=summs,n=n,cores=cores)
      prior1 <- select.prior[tt]
      pp <- c('HS','C','L')[which(c(2,4,9)==prior1)]
      #select.prior[tt] <- pp
      path1 <- paste0(path,'/testfiles')
      #for(chr0 in chr){
      test <- read_plink(paste0(path1,'/chr',chr0))
      res <- fread(paste0(path,'/posterior/chr',chr0,'_',pp,'.txt'))
      colnames(test$bim)[2] <- 'SNP'
      test$bim$order <- 1:nrow(test$bim)
      ww <- merge(test$bim,res,by='SNP')
      sig <- agtc(ww$V5,ww$V6,ww$A1.y,ww$A2.y)
      ww$THETA <- ww$THETA*sig
      X <- bedNA(test$bed[,ww$order])
      X1 <- scale(X)
      X1[which(is.na(X1))] <- 0
      corr[tt] <- cor(test$fam$V6,X1%*%ww$THETA)
      tt <- tt+1
      if(is.null(S)){S <- X1%*%ww$THETA/weight[which(c(2,4,9)==prior1)]
      }else{S <- S+X1%*%ww$THETA/weight[which(c(2,4,9)==prior1)]}
    }
    return(list(S=S,select.prior=select.prior,weight=weight,corr=corr))

  }else{
    pp <- c('HS','C','L')[which(c(2,4,9)==select.prior[1])]
    S <- NULL
    for(chr0 in chr){
      test <- read_plink(paste0(path1,'/chr',chr0))
      res <- fread(paste0(path,'/posterior/chr',chr0,'_',pp,'.txt'))
      colnames(test$bim)[2] <- 'SNP'
      test$bim$order <- 1:nrow(test$bim)
      ww <- merge(test$bim,res,by='SNP')
      sig <- agtc(ww$V5,ww$V6,ww$A1.y,ww$A2.y)
      ww$THETA <- ww$THETA*sig
      X <- bedNA(test$bed[,ww$order])
      X1 <- scale(X)
      X1[which(is.na(X1))] <- 0
      if(is.null(S)){S <- X1%*%ww$THETA}else{S <- S+X1%*%ww$THETA}
    }
    return(list(S=S))

  }

}
