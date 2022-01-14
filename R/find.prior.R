find.prior <- function(chr0,path,summs,n,cores,ncv1=5){
  nfile <- as.numeric(read.table(paste0(path,'/ldblocks/chr',chr0,'/nfile.txt')))
  rr <- mclapply(1:nfile,get.bl.cv,summs=summs,ncv1=ncv1,chr0=chr0,path=path,mc.cores =cores,MCMC=1000,BURN=500,n=n)
  setwd(path)
  dir.create('./cvres')
  save(rr,file=paste0('./cvres/chr',chr0,'.res.RData'))
  select.prior <- c()
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
    do.sc <- function(x){return(2*(x-max(x))/(max(x)-min(x))+2)}
    res <- c( mean((do.sc(y.v)-do.sc(y2))^2)+ mean((scale(y.v)-scale(y2))^2),
              mean((do.sc(y.v)-do.sc(y4))^2)+ mean((scale(y.v)-scale(y4))^2),
              mean((do.sc(y.v)-do.sc(y9))^2)+ mean((scale(y.v)-scale(y9))^2))
    prior0 <- c('Neu-HS','Neu-SpSL-C','Neu-SpSL-L')
    print(paste0('chr',chr0,' select prior: ', prior0[which.min(res)]))
    select.prior <- c(2,4,9)[which.min(res)]
  #  tt <- tt+1
 # }
  return(select.prior)
}

#ww <- find.prior(4,path,summs,n,cores=25)
