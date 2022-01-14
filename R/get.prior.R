#' @import stats
#' @import data.table
get.prior <- function(chr,path,summs,n,cores,ncv1=5){
  tt <- 1
  select.prior <- c()
  for(chr0 in chr){
    nfile <- as.numeric(read.table(paste0(path,'/ldblocks/chr',chr0,'/nfile.txt')))
    rr <- mclapply(1:nfile,get.bl.cv,summs=summs,ncv1=ncv1,chr0=chr0,path=path,mc.cores =cores,MCMC=1000,BURN=500,n=n,mc.preschedule = F)
    setwd(path)
    dir.create('./cvres')
    save(rr,file=paste0('./cvres/chr',chr0,'.res.RData'))
    #path <- paste0('/home/songs/test/testforneupred/',trait,'/cv',cv)
    # for(chr0 in chr){
    load(paste0('./cvres/chr',chr0,'.res.RData'))
    y.v <- y2 <- y4 <- y9 <- c()
    for(i in 1:length(rr)){
      y.v <- c(y.v,unlist(rr[[i]]$yy.v))
      y2 <- c(y2,unlist(rr[[i]]$yy.s2))
      y4 <- c(y4,unlist(rr[[i]]$yy.s4))
      y9 <- c(y9,unlist(rr[[i]]$yy.s9))
    }
    res2 <- c( mean((scale(y.v)-scale(y2))^2),
               mean((scale(y.v)-scale(y4))^2),
               mean((scale(y.v)-scale(y9))^2))
    print('start')
    # print(res1)
    sc <- list()
    sc[[1]] <- y2
    sc[[2]] <- y4
    sc[[3]] <- y9
    res1 <- c(ks.test(y.v,y2)$statistic,ks.test(y.v,y4)$statistic,ks.test(y.v,y9)$statistic)

    print(which.min(res2))
    ww <- c(1:3)[-which.min(res2)]
    if((cor.diff.test(y.v,sc[[which.min(res2)]],y.v,sc[[ww[1]]])$p.value.onesided<0.05)&
       (cor.diff.test(y.v,sc[[which.min(res2)]],y.v,sc[[ww[2]]])$p.value.onesided<0.05)){
      res <- res2
      print('true')
    }else if(cor.diff.test(y.v,sc[[which.min(res2)]],y.v,sc[[ww[1]]])$p.value.onesided<0.05){
      res <- res1
      res[ww[1]] <- 10
    }else if(cor.diff.test(y.v,sc[[which.min(res2)]],y.v,sc[[ww[2]]])$p.value.onesided<0.05){
      res <- res1
      res[ww[2]] <- 10
    }else{
      # print(res2)
      # print(which.min(res2))
      res <- res1
    }
    print(res)
    print(which.min(res))
    prior0 <- c('Neu-HS','Neu-SpSL-C','Neu-SpSL-L')

    print(paste0('chr',chr0,' select prior (sig): ', prior0[which.min(res)]))
    select.prior[tt] <- c(2,4,9)[which.min(res)]
    tt <- tt+1
  }
  return(select.prior)
}



