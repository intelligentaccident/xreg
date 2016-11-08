#' Predict function for xreg
#' @param object xreg-object
#' @param newdata Data.frame with new data, or list of data.frames with new data
#' @param type determines what type of predictions are wanted. Not yet implemented
#' @param return_lik_df logical argument, if true, the full data.frame with likelihood-estimates will be returned
#' @param return_vector logical argument, if true, a single vector of predicted values will be returned
#' @description Predict function for xreg-objects.
#' @author Kim Rand-Hendriksen
predict.xreg <- function(object, newdata = NULL, type = c("link", "response", "terms"), return_lik_df = FALSE, return_vector = FALSE) {
  preds <- do.call(xlik, c(list(controlList = object$controlList, dataList = newdata, return_first_df = TRUE, return_lik_df = return_lik_df), as.list(object$coef)))
  if(length(newdata) == 1 && return_vector) return(preds[[1]]$Xb)
  return(preds)
}
