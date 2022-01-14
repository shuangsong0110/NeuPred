processTEST <- function(plinkpath,path,testpath,chr){
  path1 <- paste0(path,'/testfiles')
  dir.create(path1,recursive = T)
  setwd(path1)
  for(chr0 in chr){
    system(paste0(plinkpath,' --bfile ',testpath,' --chr ',chr0,' --allow-no-sex --make-bed --out chr',chr0))
  }
}

