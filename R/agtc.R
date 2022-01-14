agtc <- function(a1,a2,b1,b2){
  sig <- rep(1,length(a1))
  for(i in 1:length(a1)){
    if((is.na(a1[i]))||(is.na(a2[i]))||(is.na(b1[i]))||(is.na(b2[i]))){
      sig[i] <- 0
    }
    else if((a1[i]==b1[i])&(a2[i]==b2[i])){
      sig[i] <- 1
    }else if((a1[i]==b2[i])&(a2[i]==b1[i])){
      sig[i] <- -1
    }
    else{
      if(b1[i]=="A"){temp1 <- "T"
      }else if(b1[i]=="T"){temp1 <- "A"
      }else if(b1[i]=="G"){temp1 <- "C"
      }else if(b1[i]=="C"){temp1 <- "G"
      }else{temp1 <- 'NA'}

      if(b2[i]=="A"){temp2 <- "T"
      }else if(b2[i]=="T"){temp2 <- "A"
      }else if(b2[i]=="G"){temp2 <- "C"
      }else if(b2[i]=="C"){temp2 <- "G"
      }else{temp2 <- 'NA'}

        if((a1[i]==temp1)&(a2[i]==temp2)){
        sig[i] <- 1
      }else if((a1[i]==temp2)&(a2[i]==temp1)){
        sig[i] <- -1
      }else{ sig[i] <- 0}
    }
    # if(i %% 10000==0){
    #   print(i)
    # }
  }
  #cat(length(sig)," variants in total. \n")
  #cat("Skip ", length(which(sig==0)), " variants due to allele code mismatch. \n")
  return(sig)
}
