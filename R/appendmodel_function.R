#' Change left-hand side of formula
#' @param object object containing formula of interest
#' @param prepend character that will precede the model expression
#' @param append character that will succede the model expression
#' @description Prepends or appends to the right-hand-side of a formula
#' @aliases appendmodel.formula
#' @author Kim Rand
#' @examples
#' appendmodel(~ a + b, "exp(", ")")
#' appendmodel(z ~ a + b, " y + ")
#' @export
appendmodel <- function(object, prepend = "", append = "") UseMethod("appendmodel")
#' @export
appendmodel.formula <- function(object, prepend = "", append = ""){
  LHS <- ""
  RHS <- as.character(object[length(object)])
  if(length(object) == 3) LHS <- as.character(object[2])
  RHS <- paste(prepend, RHS, append)
  return(formula(paste(LHS, " ~ ", RHS), env = globalenv()))
}
