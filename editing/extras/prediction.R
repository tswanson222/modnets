##########
predict.ggm <- function(object, newdata = NULL, type = 'link', ...){
  dots <- tryCatch({list(...)}, error = function(e){list()})
  if('fitobj' %in% names(object)){
    fitobj <- object$fitobj
    nn <- names(fitobj)
    if(is(fitobj[[1]], 'glmnet')){
      if(is.null(newdata)){
        data <- as.matrix(object$data)
        colnames(data) <- paste0('V', 1:ncol(data), '.')
        stop('Under development')
      } else {
        stop('Under development')
      }
    }
  } else {
    message('No models detected\nRe-fit model with saveMods = TRUE to obtain predicted values')
    out <- tryCatch({predictNet(object)}, error = function(e){NULL})
  }
}
