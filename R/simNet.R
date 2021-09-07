#' Simulate network structure and data
#'
#' Has some problems
#'
#' @param N Sample size
#' @param p Number of nodes
#' @param m moderator
#' @param m2 numeric
#' @param b1 vector or matrix?
#' @param b2 vector or matrix?
#' @param sparsity numeric
#' @param intercepts vector
#' @param nIter numeric
#' @param msym logical
#' @param onlyDat logical
#' @param pbar logical
#' @param div numeric
#' @param gibbs logical
#' @param ordinal logical
#' @param nLevels numeric
#' @param mord logical
#' @param time logical
#' @param mbinary logical
#' @param minOrd numeric
#' @param m1 numeric
#' @param m1_range numeric
#' @param m2_range numeric
#' @param modType character
#' @param lags numeric or logical
#' @param V numeric
#' @param skewErr logical or numeric
#' @param onlyNets logical
#' @param netArgs list
#' @param nCores numeric
#' @param cluster character
#' @param getChains logical
#' @param const numeric
#' @param fixedPar something
#' @param V2 numeric
#' @param ... Additional arguments.
#'
#' @return A buncha stuff
#' @export
#'
#' @examples
#' 1 + 1
simNet <- function(N = 100, p = 5, m = FALSE, m2 = .1, b1 = NULL, b2 = NULL,
                   sparsity = .5, intercepts = NULL, nIter = 250, msym = FALSE,
                   onlyDat = FALSE, pbar = TRUE, div = 10, gibbs = TRUE,
                   ordinal = FALSE, nLevels = 5, mord = FALSE, time = TRUE,
                   mbinary = FALSE, minOrd = 3, m1 = NULL, m1_range = NULL,
                   m2_range = c(.1, .3), modType = 'none', lags = NULL, V = 2,
                   skewErr = FALSE, onlyNets = FALSE, netArgs = NULL,
                   nCores = 1, cluster = 'SOCK', getChains = FALSE,
                   const = 1.5, fixedPar = NULL, V2 = 1, ...){
  t1 <- Sys.time()
  args <- tryCatch({list(...)}, error = function(e){list()})
  if(!is.null(netArgs)){list2env(netArgs[setdiff(names(netArgs), 'data')], envir = environment())}
  if(!is.null(lags) & !identical(as.numeric(lags), 0)){
    args1 <- append(list(nTime = N, nPerson = 1, nNode = p, lag = 1, GGMsparsity = sparsity), as.list(match.call())[-1])
    if(!identical(m, FALSE)){args1 <- replace(args1, 'm', ifelse(mbinary, 'binary', ifelse(mord, 'ordinal', m)))}
    if('nPerson' %in% names(args)){args1 <- replace(args1, 'nPerson', args$nPerson)}
    args1 <- append(args1, args[setdiff(names(args), names(args1))])
    out <- do.call(mlGVARsim, args1[intersect(names(args1), formalArgs('mlGVARsim'))])
    return(out)
  }
  skewErr <- ifelse(is.numeric(skewErr), skewErr, ifelse(!identical(skewErr, FALSE), 3, FALSE))
  gibbsFun <- switch(2 - identical(skewErr, FALSE), function(x){rnorm(1, x, 1)},
                     function(x){sn::rsn(1, x, alpha = skewErr)[1]})
  if(!is.null(m1_range)){if(!is.numeric(m1_range) | length(m1_range) != 2){m1_range <- NULL}}
  if(!is.numeric(m2_range) | length(m2_range) != 2){m2_range <- c(.1, .3)}
  modType <- match.arg(tolower(modType), c('none', 'full', 'partial', 'full2', 'partial2', 'zero'))
  if(grepl('2', modType) & is.null(m1)){m1 <- .5}
  if(is.null(intercepts)){
    intercepts <- rep(0, p)
  } else if('skewed' %in% intercepts){
    if(length(intercepts) == 1 | length(intercepts) > 2){intercepts <- '3'}
    intercepts <- replicate(p, sn::rsn(1, alpha = as.numeric(setdiff(intercepts, 'skewed'))))
  } else if(length(intercepts) != p){
    intercepts <- rnorm(p)
  }
  mfun <- switch(
    2 - mbinary,
    function(size = 1){sample(x = 0:1, size = size)},
    function(mean = 0){rnorm(n = 1, mean = mean)}
  )
  if(isTRUE(m)){m <- ifelse(mbinary, 1, mfun())}
  if('skewed' %in% m){
    if(length(m) == 1 | length(m) > 2){m <- '3'}
    mskew <- as.numeric(setdiff(m, 'skewed'))
    mfun <- function(mean = 0){sn::rsn(n = 1, xi = mean, alpha = mskew)[1]}
    m <- mfun()
  }
  m <- ifelse(identical(m, FALSE), 0, m)
  if(is.null(b1)){
    b1 <- simPcor(p, sparsity, constant = const, finalRange = m1_range)
    if(!is.null(fixedPar)){b1[b1 != 0] <- fixedPar[1] * sign(b1[b1 != 0])}
  }
  if(is.null(b2)){
    mat <- b2 <- diag(0, p)
    pn <- (p * (p - 1))/2
    if(V2 == 1 & m2 < 1){m2 <- ceiling(m2 * pn)}
    while(all(b2 == 0) & modType != 'zero'){
      if(m2 >= 0 & m2 < 1){
        b2 <- sample(0:1, pn, TRUE, prob = c(1 - m2, m2))
      } else if(m2 >= 1){
        b2 <- numeric(pn)
        b2[sample(1:pn, ifelse(m2 > pn, pn, round(m2)))] <- 1
      }
      b2 <- b2 * sample(c(-1, 1), pn, TRUE, prob = c(.5, .5))
      mat[upper.tri(mat)] <- b2 * runif(pn, min(m2_range), max(m2_range))
      b2 <- as.matrix(Matrix::forceSymmetric(mat))
      if(msym){b2 <- simPcor(p, sparsity = 1 - m2, constant = const, finalRange = m2_range)}
      if(!is.null(fixedPar)){
        b2[b2 != 0] <- switch(length(fixedPar), fixedPar, fixedPar[2]) * sign(b2[b2 != 0])
      }
    }
    if(all(m == 0) | modType == 'zero'){b2 <- diag(0, p)}
  }
  if(!modType %in% c('none', 'zero') & !all(b2 == 0)){
    while(ifelse(grepl('full', modType), !all(b1[b2 != 0] == 0), any(b1[b2 != 0] == 0))){
      b1 <- simPcor(p, sparsity, constant = const, finalRange = m1_range)
      if(!is.null(fixedPar)){b1[b1 != 0] <- fixedPar[1] * sign(b1[b1 != 0])}
    }
  }
  diag(b1) <- diag(b2) <- 0
  if(is.null(m1) | identical(m1, 0) | all(m == 0)){
    m1 <- rep(0, p)
  } else {
    if(isTRUE(m1)){m1 <- .5}
    if(length(m1) == 1){
      if(is.null(m1_range)){m1_range <- c(.1, .4)}
      m10 <- runif(p, min(m1_range), max(m1_range)) * sample(c(-1, 1), p, TRUE, prob = c(.5, .5))
      if(m1 >= 0 & m1 < 1){
        m10 <- m10 * sample(0:1, p, TRUE, prob = c(1 - m1, m1))
      } else if(m1 >= 1){
        m10[sample(1:p, ifelse(m1 > p, 0, round(p - m1)))] <- 0
      }
      if(grepl('2', modType) & !all(b2 == 0)){
        mm01 <- apply(b2, 1, function(z) any(z != 0))
        while(ifelse(grepl('full', modType), !all(m10[mm01] == 0), any(m10[mm01] == 0))){
          m10 <- runif(p, min(m1_range), max(m1_range))
          if(m1 >= 0 & m1 < 1){
            m10 <- m10 * sample(0:1, p, TRUE, prob = c(1 - m1, m1))
          } else if(m1 >= 1){
            m10[sample(1:p, ifelse(m1 > p, 0, round(p - m1)))] <- 0
          }
        }
      }
      m1 <- m10
      if(!is.null(fixedPar) & !all(m1 == 0)){
        m1[m1 != 0] <- fixedPar[1] * sign(m1[m1 != 0])
      }
    }
  }
  if(onlyNets & !onlyDat){
    out <- structure(list(
      b1 = b1, b2 = b2, intercepts = intercepts, m = m, m1 = m1),
      m2 = m2, modType = modType, m1_range = m1_range, m2_range = m2_range)
    return(out)
  }
  if(nCores > 1 | isTRUE(nCores)){
    if(isTRUE(nCores)){nCores <- parallel::detectCores()}
    if(grepl('Windows', sessionInfo()$running)){cluster <- 'SOCK'}
    if(tolower(cluster) != 'mclapply'){
      cluster <- match.arg(toupper(cluster), c('SOCK', 'FORK'))
      cl <- parallel::makeCluster(nCores, type = cluster)
    } else {
      cl <- nCores
    }
    sampFun <- function(case, nIter, p, m, b1, b2, intercepts, div,
                        V, m1, skewErr, getChains){
      sampling <- matrix(NA, nrow = nIter, ncol = p)
      sampling[1, ] <- rnorm(p)
      m0 <- c(ifelse(m == 0, 0, mfun(m)), numeric(nIter - 1))
      for(iter in 2:nIter){
        m0[iter] <- ifelse(m == 0, 0, mfun(m))
        vv <- switch(V, 1:p, sample(1:p, p, replace = FALSE))
        for(v in vv){
          v_mu <- 0
          v_ps <- which(b1[v, ] != 0)
          if(length(v_ps) > 0){
            for(vp in v_ps){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, sampling[v_it, vp] * b1[v, vp])
            }
          }
          if(m1[v] != 0){v_mu <- c(v_mu, m0[iter] * m1[v])}
          v_ps2 <- which(b2[v, ] != 0)
          if(length(v_ps2) > 0){
            for(vp in v_ps2){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, ((sampling[v_it, vp] * m0[iter]) * b2[v, vp]))
            }
          }
          v_mu <- intercepts[v] + sum(v_mu)
          sampling[iter, v] <- gibbsFun(v_mu)
          if(any(abs(sampling[!is.na(sampling)]) > div)){stop('Sampler diverged')}
        }
      }
      out <- list(data = sampling[nIter, ], M = m0[nIter])
      if(getChains){
        out$chains <- structure(cbind(sampling, m0), dimnames = NULL)
        if(all(m == 0)){out$chains <- out$chains[, -ncol(out$chains)]}
      }
      return(out)
    }
    if(cluster == 'SOCK'){
      parallel::clusterExport(cl, c('mfun', 'sampFun', 'gibbsFun'), envir = environment())
    }
    out <- tryCatch({
      if(pbar){
        pbapply::pboptions(type = 'timer', char = '-')
        pbapply::pblapply(1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1,
                          b2 = b2, intercept = intercepts, div = div, V = V,
                          m1 = m1, skewErr = skewErr, getChains = getChains, cl = cl)
      } else if(tolower(cluster) != 'mclapply'){
        parallel::parLapply(
          cl, 1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1, b2 = b2,
          intercepts = intercepts, div = div, V = V, m1 = m1,
          skewErr = skewErr, getChains = getChains)
      } else {
        parallel::mclapply(
          1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1, b2 = b2,
          intercepts = intercepts, div = div, V = V, m1 = m1,
          skewErr = skewErr, getChains = getChains, mc.cores = nCores)
      }},
      error = function(e){TRUE}
    )
    if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
    if(isTRUE(out)){stop('Sampler diverged')}
    data <- do.call(rbind, lapply(out, '[[', 'data'))
    M <- unlist(lapply(out, '[[', 'M'))
    if(getChains){chains <- do.call(abind::abind, c(lapply(out, '[[', 'chains'), along = 0))}
    rm(cl)
  } else if(isTRUE(gibbs) | !all(m == 0)){
    data <- matrix(NA, nrow = N, ncol = p); M <- c()
    chains <- array(NA, dim = c(N, nIter, p + as.numeric(all(m != 0))))
    if(pbar){pb <- txtProgressBar(max = N, style = 3)}
    for(case in 1:N){
      sampling <- matrix(NA, nrow = nIter, ncol = p)
      sampling[1, ] <- rnorm(p)
      m0 <- c(ifelse(m == 0, 0, mfun(m)), numeric(nIter - 1))
      for(iter in 2:nIter){
        m0[iter] <- ifelse(m == 0, 0, mfun(m))
        vv <- switch(V, 1:p, sample(1:p, p, replace = FALSE))
        for(v in vv){
          v_mu <- 0
          v_ps <- which(b1[v, ] != 0)
          if(length(v_ps) > 0){
            for(vp in v_ps){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, sampling[v_it, vp] * b1[v, vp])
            }
          }
          if(m1[v] != 0){v_mu <- c(v_mu, m0[iter] * m1[v])}
          v_ps2 <- which(b2[v, ] != 0)
          if(length(v_ps2) > 0){
            for(vp in v_ps2){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, ((sampling[v_it, vp] * m0[iter]) * b2[v, vp]))
            }
          }
          v_mu <- intercepts[v] + sum(v_mu)
          sampling[iter, v] <- gibbsFun(v_mu)
          if(any(abs(sampling[!is.na(sampling)]) > div)){stop('Sampler diverged')}
        }
      }
      if(getChains){
        chains[case, , ] <- structure(switch(
          2 - all(m == 0), sampling, cbind(sampling, m0)),
          dimnames = NULL)
      }
      data[case, ] <- sampling[nIter, ]
      M[case] <- m0[nIter]
      if(pbar){
        setTxtProgressBar(pb, case)
        if(case == N){close(pb)}
      }
    }
  } else {
    Sigma <- cov2cor(solve(diag(p) - b1))
    if(identical(skewErr, FALSE)){
      data <- mvtnorm::rmvnorm(n = N, mean = intercepts, sigma = Sigma)
    } else {
      data <- sn::rmsn(n = N, xi = intercepts, Omega = Sigma, alpha = rep(skewErr, p))
    }
  }
  out <- list(b1 = b1, b2 = b2, intercepts = intercepts)
  if(ordinal){
    for(i in 1:ncol(data)){
      ord <- c()
      while(length(unique(ord)) < minOrd){
        ord <- as.numeric(cut(data[, i], sort(c(-Inf, rnorm(nLevels - 1), Inf))))
      }
      data[, i] <- ord
    }
  }
  if(all(m == 0)){
    out$b2 <- NULL
  } else {
    if(mord & !mbinary){
      if(minOrd < 2){minOrd <- 2}
      while(length(unique(mord)) < minOrd){
        mord <- as.numeric(cut(M, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
      }
      M <- mord
    }
    data <- cbind(data, M)
  }
  out <- append(list(data = data.frame(data)), out)
  if(onlyDat){out <- data.frame(data)} else if(getChains){out$chains <- chains}
  if(!all(m == 0)){
    if(onlyDat){
      attributes(out)[c('m', 'm1')] <- list(m = m, m1 = m1)
    } else {
      out <- append(out, list(m = m, m1 = m1))
    }
    attributes(out)[c('m2', 'modType')] <- list(m2, modType)
  }
  class(out) <- c(ifelse(onlyDat, 'data.frame', 'list'), 'mgmSim')
  attr(out, 'time') <- t2 <- Sys.time() - t1
  if(time){print(Sys.time() - t1)}
  return(out)
}

##### simNet2: wrapper for simNet to prevent failures
simNet2 <- function(..., nets = NULL, maxiter = 10, pbar = TRUE, div = 10,
                    divup = TRUE, maxdiv = 10^3, errors = FALSE, nostop = FALSE){
  args <- tryCatch({list(...)}, error = function(e){list()})
  defaults <- list(pbar = pbar, div = div, time = FALSE, onlyDat = FALSE)
  args2 <- args <- replace(args, names(defaults), defaults)
  # This was the adjustment made...
  if(!'lags' %in% names(args)){
    args$onlyNets <- TRUE
    args2$netArgs <- switch(2 - is.null(nets), do.call(simNet, args), nets)
  }
  tt <- t <- out <- 0
  div0 <- div
  while(identical(out, 0)){
    out <- tryCatch({do.call(simNet, args2)}, error = function(e){0})
    t <- t + 1
    if(t == 2){
      if(!is.null(nets)){
        if(divup & div < maxdiv){
          args$div <- args2$div <- div <- div * div0
          t <- 0
        } else {
          stop('Need to generate new network matrices')
        }
      } else {
        tt <- tt + 1; t <- 0
        if(!identical(errors, FALSE)){
          if(tt == 1 & div == div0){errors <- list()}
          errors <- append(errors, setNames(
            list(args2$netArgs), paste0('tt', tt, '_div', div)))
        }
        if(tt == maxiter){
          if(divup & div < maxdiv){
            args$div <- args2$div <- div <- div * div0
            tt <- 0
          } else {
            if(nostop){
              return(list())
            } else {
              stop('Parameters may be intractable')
            }
          }
        }
        args2$netArgs <- do.call(simNet, args)
      }
    }
  }
  attributes(out)[c('div', 't', 'tt')] <- list(div, t, tt)
  if(!identical(errors, FALSE)){attributes(out)$errors <- errors}
  return(out)
}
