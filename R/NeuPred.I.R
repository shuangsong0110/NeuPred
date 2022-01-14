#' @title Main function
#' @description A unified Bayesian framework to calculate PRS based on individual-level genotype data
#'
#' @param trainpath training data path (plink bfile)
#' @param ethnic 'EUR'/'AFR'/'ASN', The ethnic of GWAS summary statistics.
#' @param plinkpath The full path to plink software.
#' @param path The full path to save the result files.
#' @param prior 'auto'/'L'/'C'/'HS'. If not specified, when external.ld=T, the prior would be set to 'L'; when external.ld=F, meaning that the LD is accurately estimated, the prior would be set to 'auto', and the algorithm will use a CV strategy to automatically select the best-performing prior for each chromosome using only training data (details are provided in the paper).
#' @param testpath The full path to tht test data. If the test data is provided, the algorithm will calculate the overlapped SNPs between test data and training data, and then present the predictive r^2 and AUC (for binary traits), and give a ROC plot. If the test data is not provided, the algorithm will derive the posterior effect sizes for all SNPs in training dataset.
#' @param MCMC MCMC iteration times. Default=10000.
#' @param BURN MCMC burning times. Default=2000.
#' @param parallel T/F.  To decide whether to run the algorithm in parallel to speed up, default=T. Note that parallel=T will improve the computational effiency, but pay attention to the limits of memory.
#' @param cores Number of cores for running MCMC when parallel=T, the optional is min(5, # available cores).
#' @param chr Default=1:22.
#' @param plot Default=F. Provide an AUC plot for binary traits.
#' @param scale Default=T. If the effect sizes are used for standardized genotypes, please set scale=T; if the effect sizes are used for raw genotypes (0,1,2), such as calculating scores with PLINK, please set scale=F
#' @param tmp Default=T.  If tmp=T, the temp files will be kept, including the LD blocks, test files, etc. If tmp=F, the temp files will be deleted after the results have been generated.
#' @return a list containing:
#' $res: The predictive r2 and AUC with selected prior.
#' $res.all: The predictive r2 and AUC with all three priors including Neuronized SpSL-Laplace prior, Neuronized SpSL-Cauchy prior, and Neuronized Horseshoe prior.
#' $select.prior: The prior automatically selected.
#' $S: estimated polygenic risk score.
#' @import EBPRS
#' @import parallel
#' @import data.table
#' @import utils
#' @import stats
#' @importFrom pROC roc
#'
#' @export
#'
#' @references
#' Song, S.,  Hou, L. and Jun, S.L. A flexible Bayesian regression approach  for accurate polygenic risk prediction. \emph{Submitted}.
#'

NeuPred.I.run <- function(trainpath,
                        ethnic='EUR',
                        plinkpath,
                        path,
                        prior='auto',
                        testpath=NULL,
                        MCMC=1e4,
                        BURN=2e3,
                        parallel=T,
                        cores=5,
                        chr=1:22,
                        plot=F,
                        scale=T,
                        tmp=T){
  #summs: snpid,chr, a1, a2, beta, p
  #LD: bfile path
  #pos: chr start stop
  trainfam <- fread(paste0(trainpath,'.fam'))
  n <- nrow(trainfam)
  rm(trainpath)
  print('If the sample size of the training data is very large, we also suggest generate summary statistics using PLINK and use NeuPred_run based on summary statistics and in-sample LD, which will save a lot of time and have little information loss.')
  res.all <- list()
  #shrink <- external.ld
  if(is.null(path)){path <- getwd()}else{dir.create(path,recursive = T)}
  if(is.null(path)){print('please input a valid path')}
  #if((is.null(LDpath))&(is.null(LDpath.chr))){stop('Please specify a valid path for LD reference.')}
  if((prior=='L')||(prior=='Laplace')){
    prior <- 9
  }else if((prior=='C')||(prior=='Cauchy')){
    prior <- 4
  }else if((prior=='HS')||(prior=='Horseshoe')){
    prior <- 2
  }
  if(is.null(testpath)){
    print('No test data are provided, thus the posterior effect sizes for all SNPs in training data will be estimated. The posterior effect sizes are stored in the folder ',path,'/posterior/ ')
    plot <- F
  }
  # vars <- colnames(summs)
  # if(!'BETA'%in%vars){
  #   if('OR'%in%vars){
  #     summs$BETA <- log(summs$OR)}else{
  #       stop("Please recheck the form of GWAS summary statistics!")
  #     }
  # }

  max.cores <- detectCores(all.tests = FALSE, logical = TRUE)
  if(cores>=max.cores){
    print(paste0(cores,' cores specified, but only ', max.cores, ' cores detected, thus ',max.cores, ' run in parallel'))
    cores <- max.cores
  }else{
    print(paste0(cores,' cores run in parallel'))
  }
  # asn <- eur <- afr <- NULL
  if(ethnic=='ASN'){
    #   data('asn')
    pos <- asn
  }else if(ethnic=='AFR'){
    #   data('afr')
    pos <- afr
  }else{
    #    data('eur')
    pos <- eur
  }
  #print(head(pos))
  # temp <-  processTRAIN.I(trainpath=trainpath,testpath=testpath,LDpath=trainpath, LDpath.chr=NULL, pos=pos,path=path,plinkpath=plinkpath,chr=chr)
  processLD.I(trainpath=trainpath,testpath=testpath,LDpath=trainpath, LDpath.chr=NULL, pos=pos,path=path,plinkpath=plinkpath,chr=chr)

  if(!is.null(testpath)){
    temp <- processTEST(plinkpath=plinkpath,path=path,testpath=testpath,chr=chr)
  }
  setwd(path)
  path2 <- paste0(path,'/posterior/')
  dir.create(path2,recursive = T)
  for(chr0 in chr){
    if(prior=='auto'){
      for(prior1 in c(2,4,9)){
        pp <- c('HS','C','L')[which(c(2,4,9)==prior1)]
        #if(!file.exists(paste0(path2,'/chr',chr0,'_',pp,'.txt'))){
        nfile <- as.numeric(read.table(paste0(path,'/ldblocks/chr',chr0,'/nfile.txt')))
        theta0 <- mclapply(1:nfile,get.bl.I,scale=scale,mc.cores=cores,path=path, chr0=chr0,MCMC=MCMC,BURN=BURN,prior=prior1,mc.preschedule=F)
        for(j in 1:nfile){
          if(j==1){
            dat <- theta0[[1]]
          }else{
            dat <- rbind(dat,theta0[[j]])
          }
        }
        dat1 <- dat[,c('SNP','CHR','A1.y','A2.y','THETA')]
        write.table(dat1,paste0(path2,'/chr',chr0,'_',pp,'.txt'),quote=F,row.names = F,col.names = T)
        #}
      }
    }else{
      pp <- c('HS','C','L')[which(c(2,4,9)==prior)]
      # if(!file.exists(paste0(path2,'/chr',chr0,'_',pp,'.txt'))){
      nfile <- as.numeric(read.table(paste0(path,'/ldblocks/chr',chr0,'/nfile.txt')))
      theta0 <- mclapply(1:nfile,get.bl.I,mc.cores=cores,path=path, chr0=chr0,MCMC=MCMC,BURN=BURN,prior=prior,mc.preschedule=F)
      for(j in 1:nfile){
        if(j==1){
          dat <- theta0[[1]]
        }else{
          dat <- rbind(dat,theta0[[j]])
        }
      }
      # print(head(dat))
      dat1 <- dat[,c('SNP','CHR','A1.y','A2.y','THETA')]

      write.table(dat1,paste0(path2,'/chr',chr0,'_',pp,'.txt'),quote=F,row.names = F,col.names = T)
      #   }
    }
  }
  if(prior=='auto'){
    select.prior <- get.prior.I(chr,path,cores,ncv1=5)
  }else{
    select.prior <- rep(prior,length(chr))
  }
  print(select.prior)
  if(prior=='auto'){auto=T}else{auto <- F}
  if(!is.null(testpath)){
    temp <- get.prs(path=path,chr=chr,select.prior=select.prior,cores=cores,auto=auto)
    S <- temp$S
    #select.prior <- temp$select.prior
    fam <- fread(paste0(path,'/testfiles/chr',chr0,'.fam'))
    type <- length(unique(fam$V6))
    if(type!=2){plot=F}
    if(plot){
      tab <- unique(fam$V6)
      if(length(tab)!=2){plot <- F}
    }
    res <- validate(S,fam$V6)
  }

  if(plot){
    rocobj1<- roc(fam$V6,S[,1])
    # print(class(S))
    # print(head(S))
    # print(class(fam$V6))
    pdf(paste0(path,'/AUC_plot.pdf'))
    plot(rocobj1,print.auc=TRUE,auc.polygon=TRUE,grid=c(0.1,0.2),grid.col=c('green','red'),
         max.auc.polygon=T,auc.polygon.col='skyblue',print.thres=F)
    dev.off()
  }
  if((prior=='auto')&(!is.null(testpath))){
    kk <- 1
    for(prior1 in c(2,4,9)){
      temp <- get.prs(path=path,chr=chr,select.prior=select.prior,cores=cores,auto=F)
      S <- temp$S
      res.all[[kk]] <- validate0(S,fam$V6)
      kk <- kk+1
    }
    if(!tmp){
      system(paste0('rm ',path,'/merge.*'))
      system(paste0('rm -r ',path,'/testfiles'))
      system(paste0('rm -r ',path,'/ldblocks'))
      system(paste0('rm -r ',path,'/cvres'))
    }
    names(res.all) <- c('HS','SpSL-C','SpSL-L')
    print(paste0('NeuPred with automatically selected prior achieves performance of predictive r2: ',round(res$r2,3)))
    print(paste0('The comprehensive results are stored in ',path,'/result.RData'))
    result=list(S=S,res=res,res.all=res.all,select.prior=select.prior)
    save(result,file=paste0(path,'/result.RData'))
    return(list(S=S,res=res,res.all=res.all,select.prior=select.prior))

  }else if((prior!='auto')&(!is.null(testpath))){
    if(!tmp){
      system(paste0('rm ',path,'/merge.*'))
      system(paste0('rm -r ',path,'/testfiles'))
      system(paste0('rm -r ',path,'/ldblocks'))
      system(paste0('rm -r ',path,'/cvres'))
    }
    prior00 <- c('Neu-HS','Neu-SpSL-C','Neu-SpSL-L')[which(c(2,4,9)==prior)]
    print(paste0('NeuPred with prior ',prior00,' achieves performance of predictive r2: ',round(res$r2,3)))
    print(paste0('The comprehensive results are stored in ',path,'/result.RData'))
    result <- list(S=S,res=res,select.prior=select.prior)
    save(result,file=paste0(path,'/result.RData'))

    return(list(S=S,res=res,select.prior=prior))
  }else if((prior=='auto')&(is.null(testpath))){
    return(list(select.prior=select.prior))
    print(paste0('No test data provided. The posterior effect sizes are saved in the folder ',path,'/posterior/, and NeuPred automatically selects ', select.prior, ' for the chromosomes provided.'))
  }else{
    print(paste0('No test data provided. The posterior effect sizes with ',prior,' prior are saved in the folder ',path,'/posterior/'))
  }


}
