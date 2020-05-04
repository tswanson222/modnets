##### netInts: get interactions from fit objects
netInts <- function(fit, n = 'temporal', threshold = FALSE, avg = FALSE, 
                    rule = 'none', r = NULL, empty = TRUE){
  rules <- c('none', 'or', 'and')
  eout <- function(fit, empty = TRUE){
    n <- tryCatch({ncol(net(fit))}, error = function(e){TRUE})
    cond <- isTRUE(empty | isTRUE(n) | is.null(n) | is.na(n))
    return(switch(2 - cond, list(), diag(0, n)))
  }
  if(is(fit, 'splitNets')){return(fit$ints)}
  if(is(fit, 'mgm')){
    m <- ifelse(!is.null(r), r, fit$call$moderators[length(fit$call$moderators)])
    mgmInt <- function(x, m){
      ni <- out <- x$interactions$indicator
      if(length(ni) == 2){
        ni <- t(apply(ni[[2]], 1, function(z) setdiff(z, m)))
        vals <- unlist(x$interactions$weightsAgg[[2]]) * x$interactions$signs[[2]]
        out <- matrix(0, ncol(net(x)) - 1, ncol(net(x)) - 1)
        for(i in seq_along(vals)){out[ni[i, 1], ni[i, 2]] <- vals[i]}
        out <- out + t(out)
      }
      return(out)
    }
    out <- mgmInt(fit, m)
    if(!is.null(r) & !is(out, 'list')){out <- out[-r, -r]}
    if(is(out, 'list')){return(eout(fit, empty))} else {return(out)}
  }
  if(is(fit, 'mgmSim')){if('b2' %in% names(fit)){return(fit$b2)} else {return(eout(fit, empty))}}
  if(is(fit, 'bootNet') | (isTRUE(attr(fit, 'resample')) & 'fit0' %in% names(fit))){fit <- fit$fit0}
  if(isTRUE(n) & isTRUE(threshold)){avg <- TRUE}
  if(isTRUE(n)){
    if(tolower(as.character(threshold)) %in% rules){rule <- threshold}
    n <- threshold <- "temporal"
  }
  n <- match.arg(tolower(n), c("temporal", "between"))
  atts <- names(attributes(fit))
  if("mlGVAR" %in% atts){fit <- switch(2 - isTRUE(n == "temporal"), fit$fixedNets, fit$betweenNet)}
  if("SURnet" %in% names(fit)){fit <- fit$SURnet}
  if("lmerVAR" %in% atts){
    stopifnot("ints" %in% names(fit$coefs))
    out <- fit$coefs$ints$coefs
    pvals <- fit$coefs$ints$Pvals
  } else if(any(c('mlVARsim', 'simMLgvar') %in% atts)){
    stopifnot("mm" %in% names(fit))
    out <- fit$mm$mb2
  } else {
    if(!'interactions' %in% names(fit)){return(eout(fit, empty))}
    out <- fit$interactions[[1]]
    if(isTRUE(attr(fit, "ggm"))){
      pvals <- out[[2]][[2]]
      out <- out[[2]][[1]]
    } else {
      pvals <- fit$interactions$pvals
    }
  }
  if(threshold != FALSE & !"mlVARsim" %in% atts){
    rule <- match.arg(tolower(rule), rules)
    if(!is.numeric(threshold)){threshold <- .05}
    if(rule == 'none'){
      out <- out * ifelse(pvals <= threshold, 1, 0)
    } else if(rule == 'or'){
      out <- out * ifelse(pvals <= threshold | t(pvals) <= threshold, 1, 0)
    } else if(rule == 'and'){
      out <- out * ifelse(pvals <= threshold & t(pvals) <= threshold, 1, 0)
    }
  }
  if(avg){out <- (t(out) + out)/2}
  if(!is.null(r)){out <- out[-r, -r]}
  if(is(out, 'list')){return(eout(fit, empty))}
  return(out)
}
