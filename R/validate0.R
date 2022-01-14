validate0 <- function (score, truey, verbose=F){
  if (length(table(truey)) == 2) {
    auc <- printAUC(score, truey)
    r2 <- cor(score, truey)^2
  #  cat("The R2 equals to ", r2, "\n")
  #  cat("The AUC equals to", auc, "\n")
    return(list(r2 = r2, AUC = auc))
  }
  else if (length(table(score)) == 2) {
    auc <- printAUC(truey, score)
    r2 <- cor(score, truey)^2
    cat("Please check the order of truey and score. \n")
   # cat("The R2 equals to ", r2, "\n")
   # cat("The AUC equals to", auc, "\n")
    return(list(r2 = r2, AUC = auc))
  }
  else if (length(table(truey)) == 3) {
    loss <- which(truey == -9)
    score <- score[-loss]
    truey <- truey[-loss]
  }
  else {
    r2 <- cor(score, truey)^2
  #  cat("The R2 equals to ", r2, "\n")
    return(list(r2 = r2))
  }
}
