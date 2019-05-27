#' Change left-hand side of formula
#' @param object object containing formula of interest
#' @param newtarget character for new left-hand side
#' @param oldtarget character to restrict to specific existing target. If NULL, any target will be changed.
#' @description Changes the left-hand side of a formula
#' @aliases change_target.formula
#' @author Kim Rand-Hendriksen
#' @examples
#' formula1 <- oldtarget ~ x + y
#' change_target(formula1, "newtarget")
#' change_target(formula1, "newtarget", "not_this_target")
#' @export
change_target <- function(object, newtarget, oldtarget = NULL) UseMethod("change_target")
change_target.formula <- function(object, newtarget, oldtarget = NULL) {
  fenv <- attr(object, ".Environment")
  if(length(as.list(object)) == 2) retval <- as.formula(paste0(newtarget, "~", as.character(object[2])))
  else {
    target <- object[[2]]
    if(is.null(oldtarget)) oldtarget <- target
    if(oldtarget != target) return(object)
    retval <- as.formula(paste0(newtarget, "~", as.character(object[3])))
  }
  attr(retval, ".Environment") <- attr(object, ".Environment")
  return(retval)
}
