#' Internal function used by xreg
#' @param ... optional arguments to be forwarded to x_mle, optim, xregControl, etc.
#' @param controlList one or more xregControl-objects (concatenated using c() ) to control maximum-likelihood estimation for different types or sets of data provided in dataList. Can alternatively provide a single formula or list of formulas, in which case an xregControl object will be generated on the fly. An xreg-object can be provided, in which case the xregControlList will be extracted. See examples
#' @param dataList Named list of data.frames for which maximum sum-likelihood should be estimated. Names must correspond to names of xregControl objects. If no names are provided, order will be used to decide which xregControl-object to use.
#' @param print_sum Optional argument to force continuous printing of the progress of fitting. Can be useful for debugging.
#' @param return_lik_df if true, will return the data.frame(s) on which likelihood are estimated
#' @param return_first_df if true, eill return the first fit to initial values
#' @param latent_classes Optional integer parameter defining the number of latent classes to be fitted
#' @param latent_class_parameter Optional character vector with names corresponding to the parameters that are to be varied in the latent classes (if latent_classes > 1)
#' @param ... optional arguments to be forwarded to x_mle, optim, xregControl, etc.
#' @description Allows a list of formulas and fitting functions to be used for each type or set of data.
#' @aliases summary.xreg, print.xreg
#' @author Kim Rand
#' @export
xlik <- function(...,
                 controlList,
                 dataList,
                 print_sum = FALSE,
                 return_lik_df = FALSE,
                 return_first_df = FALSE,
                 latent_classes = 0,
                 latent_class_parameters = character()
) {
  args <- list(...)
  
  
  
  inner_loop <- function(...,
                         formulas,
                         valueVar,
                         data_df,
                         p_fun,
                         p_aggregation_fun) {
    
    formula_loop <- function(cur_suffix = "") {
      
      for(formulanum in 1:length(formulas)) {
        formula <- formulas[[formulanum]]
        new_var <- as.character(formula[[2]])
        if(new_var == valueVar) new_var = "Xb"
        d_df[, new_var] <- with(args, with(d_df, eval(parse(text = formula))))
      }
      if(!"predict" %in% names(args)) {
        d_df$p <- p_fun(d_df)
        d_df$wp <- p_aggregation_fun(d_df)
        colnames(d_df)[match(c("Xb", "p", "wp"), colnames(d_df))] <- paste0(c("Xb", "p", "wp"), cur_suffix)  
      }
      
      return(d_df)
    }
    
    
    args <- list(...)
    #print(names(args))
    d_df <- data_df
    
    d_df <- formula_loop()
    
    return(d_df)
    
  }
  
  
  resList <- list()
  if(print_sum) {
    
    this_call <- match.call()
    #print(names(this_call))
    
    #print(names(args))
    this_call[!names(this_call) %in% names(args)] <- NULL
    
    print(unlist(as.list(this_call)))
    
  }
  
  for(dataName in names(dataList)) {
    data_df <- dataList[[dataName]]
    #print(colnames(data_df))
    control <- controlList[[dataName]]
    resList[[dataName]] <- with(control, {
      return(inner_loop(...,
                        valueVar = valueVar,
                        formulas = formulas,
                        data_df = data_df,
                        p_fun = p_fun,
                        p_aggregation_fun = p_aggregation_fun))
    })
  }
  
  sumList <- list()
  ret_sum <- 0
  
  for(i in 1:length(resList)) ret_sum <- ret_sum + (sumList[[names(resList)[i]]] <- sum(resList[[i]]$wp))
  
  if(return_first_df) return(resList)
  
  if(print_sum) {

    print(paste("Sum: ",ret_sum))
    
  }
  
  if(is.na(ret_sum)) ret_sum <- 9999999999999999999
  
  return(ret_sum)
}