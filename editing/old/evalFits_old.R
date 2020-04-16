##### evalFits: summarize iterated fitting from cluster
evalFits <- function(out1, compare = FALSE, combine = TRUE, threshold = TRUE, avg = TRUE, 
                     na = 2, getOriginal = FALSE, rmFailed = TRUE, multi = FALSE){
  if(is.logical(na)){na <- as.numeric(!na)}
  mconds <- (!identical(compare, FALSE) | na %in% 0:1) & (combine %in% c(TRUE, 'all'))
  if(!identical(multi, FALSE) & isTRUE(mconds)){
    N <- switch(2 - isTRUE(multi), as.numeric(gsub('[^0-9]', '', names(out1))), multi)
    stopifnot(is.numeric(N) & length(N) == length(out1))
    call <- replace(as.list(match.call())[-1], 'multi', list(multi = FALSE))
    out <- lapply(seq_along(out1), function(i){
      out0 <- do.call(evalFits, replace(call, 'out1', list(out1 = out1[[i]])))
      if(identical(compare, FALSE)){out0 <- out0$out2}
      return(out0)
    })
    out <- cbind.data.frame(do.call(rbind, out), N = rep(N, sapply(out, nrow)))
    out <- out[order(out$N), ]
    rownames(out) <- 1:nrow(out)
    return(out)
  }
  if(is.character(combine)){
    combine <- match.arg(tolower(combine), c('none', 'rules', 'type', 'all'))
    combine <- switch(combine, none = FALSE, all = TRUE, 
                      type = ifelse(compare, 'type', FALSE),
                      rules = ifelse(compare, 'rules', TRUE))
  }
  if(is.null(names(out1))){names(out1) <- paste0('iter', 1:length(out1))}
  if(length(unique(sapply(out1, length))) > 1){
    if(rmFailed){
      niter <- length(out1)
      atts <- attributes(out1)
      keep <- which(sapply(out1, length) == max(sapply(out1, length)))
      out1 <- out1[keep]
      mostattributes(out1) <- atts
      names(out1) <- paste0('iter', keep)
      message(paste0(niter - length(out1), '/', niter, ' iterations removed'))
    } else {
      message('evalFits failed')
      return(out1)
    }
  }
  cna <- function(out){lapply(out, function(z) replace(z, is.na(z), 0))}
  anyna <- function(x, lev = 1){switch(lev, any(sapply(x, function(z) any(is.na(z)))), any(sapply(x, anyna, lev = 1)))}
  rules <- c('OR', 'AND')
  if(!identical(compare, FALSE)){
    outs1 <- setNames(lapply(c('netInts', 'net'), function(FUN){
      setNames(lapply(rules, function(rule){
        args <- list(threshold = threshold, rule = rule)
        if(FUN == 'netInts'){args <- append(args, list(avg = avg, empty = FALSE))}
        lapply(out1, function(iter){
          lapply(iter[setdiff(names(iter), c('fitNCT', 'fitFGL'))], function(z){
            do.call(match.fun(FUN), append(args, list(fit = z)))
          })
        })
      }), rules)
    }), c('interactions', 'pairwise'))
    if(isTRUE(compare) | !is.character(compare)){
      compare <- c('Correlation', 'MAE', 'MSE')
    } else {
      inds <- c('cosine', 'correlation', 'mse', 'rmse', 'ssd', 'mae', 'msd')
      if('all' %in% compare){compare <- inds}
      compare <- Hmisc::capitalize(match.arg(tolower(compare), inds, several.ok = TRUE))
      if(any(!startsWith(compare, 'C'))){
        compare[!startsWith(compare, 'C')] <- toupper(compare[!startsWith(compare, 'C')])
      }
    }
    outs <- lapply(outs1, function(z) lapply(z, function(x){
      nn <- switch(2 - !is.null(names(x)), names(x), seq_along(x))
      do.call(rbind.data.frame, lapply(nn, function(xn){
        fits <- x[[xn]]
        names(fits) <- c('x', 'Prune', attr(out1, 'criterion'))
        if(is.character(xn)){xn <- as.numeric(gsub('[^0-9]', '', xn))}
        setNames(data.frame(sapply(compare, function(ind){
          sapply(2:length(fits), function(i){
            matrixDist(fits[[1]], fits[[i]], ind = ind, directed = !avg)
          })
        }), factor(names(fits)[-1], levels = names(fits)[-1]), xn), c(compare, 'fit', 'iter'))
      }))
    }))
    if(!identical(combine, FALSE)){
      outs <- lapply(outs, function(dat){
        d1 <- cbind.data.frame(do.call(rbind, dat), rule = factor(
          rep(names(dat), each = nrow(dat[[1]])), levels = names(dat)))
        data.frame(fit = rep(d1$fit, length(compare)), rule = rep(d1$rule, length(compare)),
                   measure = factor(rep(compare, each = nrow(d1)), levels = compare),
                   value = unlist(c(d1[, compare]), use.names = FALSE), 
                   iter = rep(d1$iter, length(compare)))
      })
      if(isTRUE(combine) | combine == 'type'){
        outs <- structure(cbind.data.frame(do.call(rbind, outs), type = factor(
          rep(Hmisc::capitalize(names(outs)), each = nrow(outs[[1]])), 
          levels = rev(Hmisc::capitalize(names(outs))))), 
          row.names = 1:(nrow(outs[[1]]) * 2))
        if(combine == 'type'){
          outs <- setNames(lapply(rules, function(z){
            z1 <- subset(outs, rule == z)[, setdiff(colnames(outs), 'rule')]
            rownames(z1) <- 1:nrow(z1)
            return(z1)
          }), rules)
        }
      }
    } else {
      outs <- lapply(outs, function(k) lapply(k, function(d1){
        data.frame(fit = rep(d1$fit, length(compare)),
                   measure = factor(rep(compare, each = nrow(d1)), levels = compare),
                   value = unlist(c(d1[, compare]), use.names = FALSE), 
                   iter = rep(d1$iter, length(compare)))
      }))
    }
    if(getOriginal){outs <- append(list(out1 = out1), outs)}
    return(outs)
  }
  out2 <- tryCatch({lapply(out1, function(z){
    x <- z$x$b2
    xx <- expand.grid(fit = names(z)[grep('[0-9]', names(z))], 
                      rule = rules, stringsAsFactors = FALSE)
    xx <- xx[order(xx$fit), ]
    f12 <- setNames(lapply(seq_len(nrow(xx)), function(i){
      mod <- z[[xx$fit[i]]]
      w0 <- c('which.lam', 'criterion')
      w1 <- which(w0 %in% names(mod$call))
      w2 <- switch(length(w1) + 1, 'Prune', mod$call$which.lam, switch(
        length(unique(mod$call[w0[w1]])), mod$call$criterion, 'CV'))
      if(identical(w2, 'min')){w2 <- 'CV'}
      thresh <- ifelse(is.logical(threshold), threshold, isTRUE(w2 == 'Prune'))
      est1 <- netInts(mod, rule = xx$rule[i], avg = avg, threshold = thresh, empty = FALSE)
      est2 <- structure(performance(est1, x, ind = ifelse(avg, 'between', 'beta'), mcc = TRUE), criterion = w2)
      return(est2)
    }), apply(xx, 1, paste, collapse = '_'))
    xxx <- append(f12, list(
      f3 = structure(z$fitNCT$performance, criterion = 'NCT'), 
      f4 = structure(z$fitFGL$performance, criterion = 'FGL')))
    return(xxx)
  })}, error = function(e){list()})
  outs <- lapply(seq_len(anyna(out2, 2) + 1), function(i){
    out3 <- tryCatch({lapply(seq_len(unique(sapply(out2, length))), function(z){
      z1 <- data.frame(do.call(rbind, lapply(out2, '[[', z)))
      nn <- c('Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'MCC')
      z2 <- unlist(c(z1), use.names = FALSE)
      z3 <- cbind.data.frame(value = z2, measure = factor(rep(nn, each = nrow(z1)), levels = nn))
      rownames(z3) <- 1:nrow(z3)
      return(z3)
    })}, error = function(e){list()})
    if(i == 2){out3 <- cna(out3)}
    out4 <- tryCatch({setNames(lapply(seq_along(rules), function(z){
      n <- length(out3)
      ind <- c(seq(z, (n - 2), length(rules)), (n - 1):n)
      nn <- unique(unique(lapply(out2, function(zz) unname(sapply(zz, attr, 'criterion'))))[[1]])
      z1 <- do.call(rbind.data.frame, out3[ind])
      z1 <- cbind.data.frame(z1, type = factor(
        rep(nn, each = nrow(out3[[1]])), levels = nn))
      attr(z1, 'box') <- TRUE
      return(z1)
    }), rules)}, error = function(e){list()})
    out5 <- tryCatch({lapply(out4, function(z){
      z1 <- aggregate(z$value, list(type = z$type, measure = z$measure), FUN = function(x){
        c(mean = mean(x, na.rm = TRUE), 
          se = sd(x, na.rm = TRUE)/sqrt(length(x)), 
          sd = sd(x, na.rm = TRUE),
          se2 = (sd(x, na.rm = TRUE)/sqrt(length(x))) * qnorm(.975))
      })
      z1 <- cbind.data.frame(z1[, 1:2], z1[, 3])
      attr(z1, 'aggregate') <- TRUE
      return(z1)
    })}, error = function(e){list()})
    return(setNames(list(out4, out5), paste0('out', 2:3)))
  })
  if(combine){
    outs <- tryCatch({lapply(outs, function(x) lapply(x, function(dat){
      atts <- setdiff(unique(sapply(dat, function(z) names(attributes(z)))), c('names', 'class', 'row.names'))
      dat <- structure(cbind.data.frame(do.call(rbind, dat), rule = factor(
        rep(names(dat), each = nrow(dat[[1]])), levels = names(dat))), 
        row.names = 1:sum(sapply(dat, nrow)))
      attributes(dat)[atts] <- TRUE
      return(dat)
    }))}, error = function(e){outs})
  }
  x01 <- paste0('x', 0:1)
  outs <- switch(length(outs), outs[[1]], setNames(outs, x01))
  if(identical(names(outs), x01) & na %in% 0:1){outs <- outs[[na + 1]]}
  if(getOriginal){outs <- append(outs, list(out1 = out1))}
  return(outs)
}
