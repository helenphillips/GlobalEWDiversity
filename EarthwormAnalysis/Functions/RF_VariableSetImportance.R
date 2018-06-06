## Using this page
## https://stats.stackexchange.com/questions/92419/relative-importance-of-a-set-of-predictors-in-a-random-forests-classification-in
## Although in comments it says it is not quite right, this page confirms the correct approch
## https://stats.stackexchange.com/questions/311488/summing-feature-importance-in-scikit-learn-for-a-set-of-features

rf.obj <- spR_rf
var.share <- function(rf.obj, members) {
  count <- table(rf.obj$forest$bestvar)[-1]
  names(count) <- names(rf.obj$forest$ncat)
  share <- count[members] / sum(count[members]) ## weighting
  return(share)
}

group.importance <- function(rf.obj, groups) {
  var.imp <- as.matrix(sapply(groups, function(g) {
    sum(importance(rf.obj, 2)[g, ]*var.share(rf.obj, g))
  }))
  colnames(var.imp) <- "MeanDecreaseGini"
  return(var.imp)
}

######## My own code for creating a plot??

# add order of importances
OrderImportance <- function(groupedImportances){
  order <- order(groupedImportances)
  return(order)
}
