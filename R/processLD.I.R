#' @import data.table
#' @import utils
processLD.I <- function(trainpath,testpath=NULL,LDpath,LDpath.chr=NULL,pos,path=NULL,plinkpath,chr){
  summs <- fread(paste0(trainpath,'.bim'))
  colnames(summs)[2] <- 'SNP'
  if(!is.null(path)){setwd(path)}else{path <- getwd()}
  if(is.null(testpath)){snp <- summs$SNP}else{
    testbim <- fread(paste0(testpath,'.bim'))
    snp <- intersect(summs$SNP,testbim$V2)
  }

  write.table(snp,'snp.txt',quote=F,row.names=F,col.names=F)
  if(!is.null(LDpath)){
    system(paste0(plinkpath, '  --bfile ',LDpath,
                  ' --extract snp.txt --allow-no-sex --make-bed --out merge'),show.output.on.console=F)
    for(chr0 in chr){
      dir.create(paste0(path,'/ldblocks/chr',chr0),recursive = T)
      setwd(paste0(path,'/ldblocks/chr',chr0))
      temp <- length(which(pos$chr==chr0))
      if(temp>0){
        bed <- pos[which(pos$chr==chr0),]
      }else{
        temp <- length(which(pos$chr==paste0('chr',chr0)))
        if(temp>0){
          bed <- pos[which(pos$chr==paste0('chr',chr0)),]
          bed$chr <- rep(chr0,nrow(bed))
        }else{
          print(paste0('no pos file for chr ',chr0))
          next
        }
      }

      ldbim <- fread(paste0(path,'/merge.bim'))
      snp.pos <- sort(ldbim$V4[ldbim$V1==chr0])
      print(head(snp.pos))
      print(head(bed))
      index <- rep(0,length(snp.pos))
      chr1 <- stop1 <- start1 <- c()
      w <- 1
      len <- 0
      indw <- 1
      start1[1] <- bed$start[1]
      for(k in 1:length(snp.pos)){
        while((index[k]==0)&(w<=nrow(bed))){
          if((snp.pos[k]>=bed$start[w])&(snp.pos[k]<bed$stop[w])){
            index[k] <- indw
            len <- len+1
          }else{
            w <- w+1
            if(len>100){
              stop1[indw] <- bed$stop[w-1]
              indw <- indw+1
              start1[indw] <- bed$start[w]
              len <- 0
            }}}}
      stop1[length(start1)] <- bed$stop[nrow(bed)]
      bed <- data.frame(chr=rep(chr0,length(start1)),start=start1,stop=stop1,ind=rep('R1',length(start1)))

      for(k in 1:nrow(bed)){
        write.table(bed[k,],paste0('region_',k,'.txt'),quote=F,row.names = F,col.names = F)
        # system(paste0('/home/songs/plink1.9/plink '))
        system(paste0(plinkpath, '  --bfile ',path,'/merge --extract range ',path,'/ldblocks/chr',chr0,'/region_',k,'.txt ',
                      ' --allow-no-sex --make-bed --out LD_region_',k),show.output.on.console=F)
      }
      write.table(nrow(bed),'nfile.txt',quote=F,row.names=F,col.names=F)
    }
  }else{
    for(chr0 in chr){
      LDpath <- paste0(LDpath.chr,chr0)
      system(paste0(plinkpath, '  --bfile ',LDpath,
                    ' --extract snp.txt --allow-no-sex --make-bed --out merge'),show.output.on.console=F)

      dir.create(paste0(path,'/ldblocks/chr',chr0),recursive = T)
      setwd(paste0(path,'/ldblocks/chr',chr0))
      temp <- length(which(pos$chr==chr0))
      if(temp>0){
        bed <- pos[which(pos$chr==chr0),]
      }else{
        temp <- length(which(pos$chr==paste0('chr',chr0)))
        if(temp>0){
          bed <- pos[which(pos$chr==paste0('chr',chr0)),]
          bed$chr <- rep(chr0,nrow(bed))
        }else{
          print(paste0('no pos file for chr ',chr0))
          next
        }
      }

      ldbim <- fread(paste0(path,'/merge.bim'))
      snp.pos <- sort(ldbim$V4[ldbim$V1==chr0])
      index <- rep(0,length(snp.pos))
      chr1 <- stop1 <- start1 <- c()
      w <- 1
      len <- 0
      indw <- 1
      start1[1] <- bed$start[1]
      for(k in 1:length(snp.pos)){
        while((index[k]==0)&(w<=nrow(bed))){
          if((snp.pos[k]>=bed$start[w])&(snp.pos[k]<bed$stop[w])){
            index[k] <- indw
            len <- len+1
          }else{
            w <- w+1
            if(len>100){
              stop1[indw] <- bed$stop[w-1]
              indw <- indw+1
              start1[indw] <- bed$start[w]
              len <- 0
            }}}}
      stop1[length(start1)] <- bed$stop[nrow(bed)]
      bed <- data.frame(chr=rep(chr0,length(start1)),start=start1,stop=stop1,ind=rep('R1',length(start1)))

      for(k in 1:nrow(bed)){
        write.table(bed[k,],paste0('region_',k,'.txt'),quote=F,row.names = F,col.names = F)
        # system(paste0('/home/songs/plink1.9/plink '))
        system(paste0(plinkpath, '  --bfile ',path,'/merge --extract range ',path,'/ldblocks/chr',chr0,'/region_',k,'.txt ',
                      ' --allow-no-sex --make-bed --out LD_region_',k),show.output.on.console=F)
      }
      write.table(nrow(bed),'nfile.txt',quote=F,row.names=F,col.names=F)
    }

  }
}
