#' @import methods
#' @import ROCR

printAUC <- function(S,testy){
  #library(ROCR)
  pred <- prediction(S, testy)
  auc <- performance(pred,'auc')
  auc=unlist(slot(auc,"y.values"))
  if(auc<0.5){
    pred <- prediction(S, -testy)
    auc <- performance(pred,'auc')
    auc=unlist(slot(auc,"y.values"))
  }
  return(auc)
}
