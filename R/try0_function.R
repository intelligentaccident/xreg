#' Wrapper for try which returns second argument if not succesful
#' @param expr expression to be tested
#' @param ifError return value if test fails
#' @description Returns second argument if expression in first argument fails.
#' @author Kim Rand-Hendriksen
#' @examples
#' try0(as.formula(list()))
#' try0(as.formula("a ~ b))
try0 <- function(expr, ifError = NULL) {
  tr <- try(expr, TRUE)
  succ <- (class(tr) != "try-error")
  if(!succ) return(ifError)
  return(tr)
}
