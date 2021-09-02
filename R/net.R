#' Get adjacency matrices from fit objects
#'
#' Used all the time
#'
#' @param fit model
#' @param n character
#' @param threshold numeric or logical
#' @param rule chatacter
#' @param binary logical
#' @param nodewise logical
#' @param d numeric
#' @param r something
#' @param avg logical
#' @param empty logical
#' @param mselect logical
#'
#' @return A network matrix
#' @export
#'
#' @examples
#' 1 + 1
net <- function(fit, n = "beta", threshold = FALSE, rule = "OR",
                binary = FALSE, nodewise = FALSE, d = 14, r = NULL){
  if(inherits(fit, c('splitNets', 'try-error'))){return(NULL)}
  if(is(fit, 'bootNet') | (isTRUE(attr(fit, 'resample')) & 'fit0' %in% names(fit))){fit <- fit$fit0}
  if(is(fit, 'mgmSim')){return(fit[[grep('^trueNet$|^b1$|^adjMat$', names(fit))]])}
  if(is(fit, "matrix")){return(fit)}
  if(is(fit, "mgm")){
    out <- fit$pairwise$wadj * replace(fit$pairwise$signs, is.na(fit$pairwise$signs), 0)
    if(!is.null(r)){out <- out[-r, -r]}
    return(out)
  }
  if(isTRUE(n)){
    if(tolower(threshold) %in% c('and', 'or')){rule <- threshold}
    n <- threshold <- "beta"
  }
  n <- match.arg(tolower(n), c(
    "beta", "contemporaneous", "between", "pdc", "pcc", "kappa", "temporal"))
  n1 <- switch(n, between = "", beta = , pdc = , temporal = "$temporal", "$contemporaneous")
  n2 <- switch(n, kappa = "kappa", pdc = "PDC$adjMat", "adjMat")
  n3 <- paste0("fit", n1, "$", n2)
  if(isTRUE(attr(fit, "mlGVAR"))){fit <- if(n == "between"){fit$betweenNet} else {fit$fixedNets}}
  if("SURnet" %in% names(fit)){fit <- fit$SURnet}
  if(isTRUE(attr(fit, "ggm"))){n3 <- ifelse(nodewise, "fit$nodewise$adjNW", "fit$adjMat")}
  if("fixedPCC" %in% names(fit)){
    threshold <- FALSE
    n <- switch(n, contemporaneous = "pcc", temporal = "beta", n)
    if(n != "between" & "fixedResults" %in% names(fit)){fit <- fit$fixedResults}
    out <- fit[[grep(n, tolower(names(fit)))[1]]]
    if(n == "beta" & ncol(out) != nrow(out)){out <- out[, -1]}
    if(n == "between"){diag(out) <- 0}
    if(n == "pdc"){out <- t(out)}
    if(!is.null(d)){out <- round(out, d)}
  } else {
    if(isTRUE(attr(fit, "lmerVAR")) & n == "between"){n3 <- paste0("fit$between$adjMat")}
    out <- eval(parse(text = n3))
  }
  if(threshold != FALSE){
    if(!is.numeric(threshold)){threshold <- .05}
    rule <- match.arg(tolower(rule), c("or", "and"))
    atts <- names(attributes(fit))
    if("SURnet" %in% atts){
      stopifnot(!n == "between")
      n4 <- ifelse(n %in% c("beta", "temporal", "pdc"), "temporal", "contemporaneous")
      n5 <- ifelse(n4 == "temporal", "fit$temporal$coefs", "fit$contemporaneous")
      pvals <- eval(parse(text = n5))$pvals
      if(ncol(pvals) != nrow(pvals)){
        if(any(grepl("lag", colnames(pvals)))){
          colnames(pvals) <- gsub("[.]lag1[.]$", "", colnames(pvals))
        }
        pvals <- pvals[, intersect(colnames(pvals), gsub("[.]y$", "", rownames(pvals)))]
      }
    } else {
      n4 <- ifelse(n %in% c("beta", "temporal", "pdc"), "Beta",
                   ifelse(n == "between", "gammaOmega", "gammaTheta"))
      pvals <- if("ggm" %in% atts){fit$nodewise$pvalsNW} else {fit$coefs[[n4]]$Pvals}
      if("ggm" %in% atts){n4 <- ifelse(nodewise & rule == "or", "Beta", "")}
    }
    if(any(is.na(pvals))){pvals[is.na(pvals)] <- 1}
    if(n4 %in% c("temporal", "Beta")){
      out <- out * ifelse(pvals <= threshold, 1, 0)
    } else if(rule == "or"){
      out <- out * ifelse(pvals <= threshold | t(pvals) <= threshold, 1, 0)
    } else if(rule == "and"){
      out <- out * ifelse(pvals <= threshold & t(pvals) <= threshold, 1, 0)
    }
  }
  if(binary){out <- abs(sign(out))}
  if(!is.null(r)){out <- out[-r, -r]}
  return(out)
}

#' @rdname net
#' @export
netInts <- function(fit, n = 'temporal', threshold = FALSE, avg = FALSE,
                    rule = 'none', r = NULL, empty = TRUE, mselect = NULL){
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
  } else if(!'interactions' %in% names(fit)){
    out <- tryCatch({mmat(fit = fit, m = mselect)}, error = function(e){TRUE})
    if(isTRUE(out)){return(eout(fit, empty))}
    pvals <- out$pvals
    out <- out$betas
  } else {
    out <- fit$interactions[[1]]
    if(isTRUE(attr(fit, "ggm"))){
      pvals <- out[[2]][[2]]
      out <- out[[2]][[1]]
    } else {
      rule <- 'none'
      pvals <- fit$interactions$pvals
      if(ncol(out) != nrow(out) & !is.null(mselect)){
        if(isTRUE(mselect)){mselect <- fit$call$moderators[1]}
        m1 <- paste0(':', mselect); m2 <- paste0(mselect, ':')
        mm <- paste0(c(m1, m2), collapse = '|')
        my <- paste0(gsub('[.]y$', '', rownames(out)), collapse = '|')
        mb1 <- out[, grep(mm, colnames(out))]
        mb1 <- mb1[, grep(my, colnames(mb1))]
        mp1 <- pvals[, grep(mm, colnames(pvals))]
        mp1 <- mp1[, grep(my, colnames(mp1))]
        out <- mb1
        pvals <- mp1
      }
    }
  }
  if(threshold != FALSE & !"mlVARsim" %in% atts){
    rule <- match.arg(tolower(rule), rules)
    if(isTRUE(attr(fit, 'ggm')) & rule == 'none'){rule <- 'or'}
    if(!is.numeric(threshold)){threshold <- .05}
    if(rule == 'none' | isTRUE(attr(fit, 'SURnet'))){
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
