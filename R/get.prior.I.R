#' @import stats
#' @import data.table
get.prior.I <- function(chr,path,cores,ncv1=5){
  tt <- 1
  select.prior <- c()
  for(chr0 in chr){
    nfile <- as.numeric(read.table(paste0(path,'/ldblocks/chr',chr0,'/nfile.txt')))
    rr <- mclapply(1:nfile,get.bl.cv.I,ncv1=ncv1,chr0=chr0,path=path,mc.cores =cores,MCMC=1000,BURN=500,mc.preschedule = F)
    setwd(path)
    dir.create('./cvres')
    save(rr,file=paste0('./cvres/chr',chr0,'.res.RData'))
    #path <- paste0('/home/songs/test/testforneupred/',trait,'/cv',cv)
    # for(chr0 in chr){
    load(paste0('./cvres/chr',chr0,'.res.RData'))
    y.v <- y2 <- y4 <- y9 <- c()
    rr2 <- rr4 <- rr9 <- c()
    for(ncv in 1:ncv1){
    for(i in 1:length(rr)){
      if(i ==1){
        y.v <- rr[[i]]$yy.v[[ncv]]
        y2 <- rr[[i]]$yy.s2[[ncv]]
        y4 <- rr[[i]]$yy.s4[[ncv]]
        y9 <- rr[[i]]$yy.s9[[ncv]]
      }else{
       # y.v <- y.v+rr[[i]]$yy.v[[ncv]]
        y2 <- y2+rr[[i]]$yy.s2[[ncv]]
        y4 <- y4+ rr[[i]]$yy.s4[[ncv]]
        y9 <- y9+ rr[[i]]$yy.s9[[ncv]]
      }
    }
      rr2[ncv] <- cor(y.v,y2)^2
      rr4[ncv] <- cor(y.v,y4)^2
      rr9[ncv] <- cor(y.v,y9)^2
    }
    res <- c(mean(rr2),mean(rr4),mean(rr9))
    print(res)
    print(which.max(res))
    prior0 <- c('Neu-HS','Neu-SpSL-C','Neu-SpSL-L')

    print(paste0('chr',chr0,' select prior (sig): ', prior0[which.max(res)]))
    select.prior[tt] <- c(2,4,9)[which.max(res)]
    tt <- tt+1
  }
  return(select.prior)
}



