### ======================================================================== ###
### ========================== SIMNETS FUNCTIONS =========================== ###
### ======================================================================== ###
##### SimNets: simulate multilevel GVAR data to evaluate parameter recovery
SimNets <- function(N, nsims = 10, ntime = 50, p = 3, m = "random", FUN = "mlGVAR",
                    inds = c("cosine", "mse", "mae"), saveMods = FALSE, 
                    betweenArgs = list(moderators = NULL), d = 14,
                    verbose = TRUE, ...){
  t1 <- Sys.time()
  A0 <- list(...)
  A1 <- append(list(nTime = ntime, nPerson = N, nNode = p, m = m), A0)
  A1 <- A1[intersect(formalArgs("mlGVARsim"), names(A1))]
  A2 <- append(list(data = NULL, m = p + 1, betweenArgs = betweenArgs, verbose = FALSE), A0)
  FUN <- match.arg(FUN, c("mlGVAR", "lmerVAR"))
  FUNargs <- setdiff(formalArgs(FUN), "...")
  if("selectFUN" %in% names(A2) & FUN == "mlGVAR"){
    if(isTRUE(A2$selectFUN)){A2$selectFUN <- "varSelect"}
    FUNargs <- setdiff(union(formalArgs(switch(match.arg(
      A2$selectFUN[1], c("varSelect", "resample", "stability", "bootstrap")), 
      varSelect = "varSelect", "resample")), FUNargs), "...")
  }
  saveVars <- isTRUE("sampMethod" %in% FUNargs)
  A2 <- A2[intersect(FUNargs, names(A2))]
  if(is.null(m)){A2$m <- NULL}
  nit <- paste0("iter", 1:nsims)
  allNets <- c("temporal", "contemporaneous", "between")
  output <- setNames(lapply(seq_along(N), function(z){
    A1$nPerson <- N[z]
    if(verbose){
      cat("\n----------------------\n SAMPLE SIZE: N = ", N[z], "\n----------------------\n", sep = "")
      message("Simulating datasets")
      pb1 <- txtProgressBar(max = nsims, style = 3)  
    }
    X <- setNames(lapply(seq_len(nsims), function(i){
      x <- do.call(mlGVARsim, A1)
      nx <- setNames(lapply(allNets, function(ii){net(x, ii, d = d)}), allNets)
      if(!is.null(m)){nx$interactions <- netInts(x)}
      if(verbose){setTxtProgressBar(pb1, i)}
      return(list(x = x, nets = list(nets0 = nx)))
    }), nit)
    if(verbose){
      close(pb1)
      message("Fitting models")
      pb2 <- txtProgressBar(max = nsims, style = 3)
    }
    Y <- setNames(lapply(seq_len(nsims), function(i){
      A2$data <- switch(2 - is.null(m), X[[i]]$x$data, X[[i]]$x$mm$data)
      x <- do.call(match.fun(FUN), A2)
      nx0 <- setNames(lapply(allNets, function(ii) net(x, ii)), allNets)
      nx1 <- setNames(lapply(allNets, function(ii) net(x, ii, T)), allNets)
      nx2 <- setNames(lapply(allNets, function(ii) net(x, ii, T, "and")), allNets)
      if(!is.null(m)){
        nx0$interactions <- netInts(x)
        nx1$interactions <- nx2$interactions <- netInts(x, T)
      }
      if(verbose){setTxtProgressBar(pb2, i)}
      return(list(x = x, nets = list(nets0 = nx0, nets1 = nx1, nets2 = nx2)))
    }), nit)
    if(verbose){close(pb2)}
    out <- setNames(lapply(seq_len(nsims), function(n1){
      setNames(lapply(seq_along(Y[[n1]]$nets), function(n2){
        indout1 <- do.call(rbind, lapply(seq_along(Y[[n1]]$nets[[n2]]), function(n3){
          sapply(inds, function(n4){matrixDist(
            Y[[n1]]$nets[[n2]][[n3]], X[[n1]]$nets$nets0[[n3]], 
            n4, isTRUE(n3 %in% c(1, 4)))})
        }))
        indout2 <- do.call(rbind, lapply(seq_along(Y[[n1]]$nets[[n2]]), function(n3){
          performance(Y[[n1]]$nets[[n2]][[n3]], X[[n1]]$nets$nets0[[n3]])
        }))
        netind <- names(Y[[n1]]$nets[[n2]])
        indout <- cbind.data.frame(net = netind, iter = rep(n1, length(netind)), indout1, indout2)
        return(indout)
      }), paste0("nets", 0:2))
    }), nit)
    out0 <- do.call(rbind, lapply(out, '[[', 1))
    out1 <- do.call(rbind, lapply(out, '[[', 2))
    out2 <- do.call(rbind, lapply(out, '[[', 3))
    rownames(out0) <- rownames(out1) <- rownames(out2) <- 1:nrow(out0)
    results <- list(nets0 = out0, nets1 = out1, nets2 = out2)
    if(saveVars){Y1 <- lapply(Y, function(z) z$x$varMods)}
    if(!saveMods){
      X <- lapply(X, '[[', "nets")
      Y <- lapply(Y, '[[', "nets")
    }
    output <- list(sims = X, fits = Y, results = results)
    if(saveVars){output$varMods <- Y1}
    attr(output, "allNets") <- allNets
    return(output)
  }), paste0("N", N))
  outcall <- as.list(match.call())[-1]
  if(any(sapply(outcall, class) %in% c("symbol", "name"))){
    oc <- which(sapply(outcall, class) %in% c("symbol", "name"))
    for(j in seq_along(oc)){
      outcall[[oc[j]]] <- as.character(outcall[[oc[j]]])
      if(outcall[[oc[j]]] %in% c("T", "F")){
        outcall[[oc[j]]] <- as.logical(outcall[[oc[j]]])
      }
    }
  }
  if(any(sapply(outcall, class) == "call")){
    oc <- which(sapply(outcall, class) == "call")
    for(j in seq_along(oc)){outcall[[oc[j]]] <- eval(outcall[[oc[j]]])}
  }
  outcall <- append(outcall, list(
    allNets = switch(2 - is.null(m), allNets, c(allNets, "interactions"))))
  output <- append(list(call = outcall), output)
  attr(output, "time") <- Sys.time() - t1
  if(verbose){cat("\n"); print(Sys.time() - t1)}
  return(output)
}

##### TEMPORARY 
threeFits <- function(dat, m, scale = TRUE, saveMods = FALSE, ...){
  args <- tryCatch({list(...)}, error = function(e){list()})
  args[c('scale', 'saveMods', 'lags')] <- list(scale, saveMods, 1)
  args <- args[intersect(names(args), formalArgs('fitNetwork'))]
  out <- setNames(lapply(1:3, function(z){
    X <- switch(z, dat[, -m], dat, dat)
    M <- switch(z, NULL, NULL, m)
    covars <- switch(z, NULL, m, NULL)
    args0 <- replace(args, c('data', 'moderators', 'covariates'), list(data = X, moderators = M, covariates = covars))
    do.call(fitNetwork, args0)
  }), paste0('fit', 1:3))
  return(out)
}
getVars <- function(dat, m, inds = c('CV', 'AIC', 'BIC', 'EBIC'), 
                    vars = NULL, fit = TRUE, ex = NULL, scale = TRUE,
                    saveMods = FALSE, nCores = 1, ...){
  args <- tryCatch({list(...)}, error = function(e){list()})
  input <- list(dat, m, 1, scale, saveMods)
  newvars <- FALSE
  if(is.null(vars)){
    args1 <- args[intersect(names(args), formalArgs('varSelect'))]
    args1[c('data', 'm', 'lags', 'scale')] <- input[1:4]
    if(nCores > 1){
      pbapply::pboptions(type = 'timer', char = '-')
      vars <- pbapply::pblapply(setdiff(inds, ex), function(z){
        args0 <- replace(args1, 'criterion', list(criterion = z))
        do.call(varSelect, args0)
      }, cl = nCores)
    } else {
      vars <- lapply(setdiff(inds, ex), function(z){
        args0 <- replace(args1, 'criterion', list(criterion = z))
        do.call(varSelect, args0)
      })
    }
    newvars <- TRUE
  }
  if(fit & !is.null(vars)){
    args2 <- args[intersect(names(args), formalArgs('fitNetwork'))]
    args2[c('data', 'moderators', 'lags', 'scale', 'saveMods')] <- input
    pbapply::pboptions(type = 'timer', char = '-')
    fits <- pbapply::pblapply(vars, function(z){
      args0 <- replace(args2, 'type', list(type = z))
      do.call(fitNetwork, args0)
    }, cl = nCores)
  }
  out <- list()
  if(newvars){out$vars <- vars}
  if(fit){out$fits <- fits}
  return(out)
}
intSelect <- function(x1, x2 = NULL, len = FALSE){
  if(!is.null(x2)){
    n <- 8 - length(x2$fits)
    x1 <- append(x1, setNames(x2$fits, paste0('fit', n:7)))
    if(!len){return(x1)}
  }
  k <- lapply(lapply(x1, selected), '[[', 'ints')
  k <- lapply(lapply(k, unlist), function(z) setdiff(z, ''))
  k
}


##### simNet: simulate network structures and data
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
  args$onlyNets <- TRUE
  args2$netArgs <- switch(2 - is.null(nets), do.call(simNet, args), nets)
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


##### ordinalize: turn continuous data into ordinal data
ordinalize <- function(data, m = NULL, nLevels = 5, thresholds = NULL, 
                       mthresh = NULL, mord = TRUE, minOrd = 3){
  if(!is(data, 'matrix') & !is(data, 'data.frame') & is(data, 'numeric')){data <- matrix(data, ncol = 1)}
  if(!is.null(m)){
    if(isTRUE(m)){m <- ncol(data)}
    m0 <- m
    m <- data[, m0]
    data <- data[, -m0]
    if(mord){
      tick <- 0
      while(length(unique(mord)) < minOrd){
        thresh <- switch(2 - is.null(mthresh), rnorm(nLevels - 1), unlist(mthresh))
        mord <- as.numeric(cut(m, sort(c(-Inf, thresh, Inf))))
        tick <- tick + 1
        if(tick == 10 & !is.null(mthresh)){mthresh <- NULL}
      }
      m <- mord
    }
  }
  for(vv in 1:ncol(data)){
    tick <- ord <- 0
    while(length(unique(ord)) < minOrd){
      thresh <- switch(2 - is.null(thresholds), rnorm(nLevels - 1), thresholds[[vv]])
      ord <- as.numeric(cut(data[, vv], sort(c(-Inf, thresh, Inf))))
      tick <- tick + 1
      if(tick == 10 & !is.null(thresholds)){thresholds <- NULL}
    }
    data[, vv] <- ord
  }
  if(!is.null(m)){data$M <- m}
  return(data)
}

##### summarize: aggregate simulation output
summarize <- function(obj, ind = "mean", ...){
  N <- names(obj)[grepl("N", names(obj))]
  fun <- match.fun(match.arg(ind, c("mean", "median", "sd", "var", "quantile")))
  output <- setNames(lapply(N, function(n) setNames(lapply(obj[[n]]$results, function(z){
    z <- do.call(rbind, lapply(obj$call$allNets, function(i){
      i <- subset(z, net == i)[, -(1:2)]
      apply(i, 2, fun, na.rm = TRUE, ...)
    }))
    rownames(z) <- obj$call$allNets
    return(z)
  }), names(obj[[n]]$results))), N)
  if(length(output) == 1){output <- output[[1]]}
  return(output)
}

##### zeros: return all networks or get the number of zeros for each parameter
zeros <- function(obj, type = "sims", net = "temporal", ind = 2, n = 1, all = FALSE){
  if(is.numeric(type)){type <- c("sims", "fits")[type]}
  if(!isTRUE(tryCatch({match.arg(
    type, c("temporal", "contemporaneous", "between", "interactions"))}, 
    error = function(e){TRUE}))){
    net <- type
    type <- "sims"
  }
  type <- match.arg(type, c("sims", "fits"))
  if(type == "sims"){ind <- 1}
  if(any(grepl("N", names(obj)))){
    allNets <- obj$call$allNets
    if(is.numeric(n)){
      if(n < length(obj)){
        obj <- obj[[n + 1]]
      } else {
        obj <- obj[[grep(paste0("^", n), gsub("N", "", names(obj)))]]
      }
    } else {
      obj <- obj[[grep(n, names(obj))]]
    }
  } else {
    allNets <- attr(obj, "allNets")
  }
  if(is.numeric(net)){net <- allNets[net]}
  net <- match.arg(net, allNets)
  out0 <- do.call(abind::abind, c(lapply(obj[[type]], function(z){
    z[[ind]][[net]]}), along = 3))
  out1 <- t(sapply(seq_len(dim(out0)[1]), function(i){
    sapply(seq_len(dim(out0)[2]), function(j){
      sum(out0[i, j, ] == 0)})}))
  if(all){return(out0)} else {return(out1)}
}

##### sumPlot
sumPlot <- function(obj, ind = "cosine", nets = 1, n = 1, ylim = c(0, 1), ...){
  if(any(grepl("N", names(obj)))){
    if(length(obj) == 2){
      obj <- obj[[2]]
    } else {
      if(is.numeric(n)){
        if(n < length(obj)){
          obj <- obj[[n + 1]]
        } else {
          obj <- obj[[grep(paste0("^", n), gsub("N", "", names(obj)))]]
        }
      } else {
        obj <- obj[[grep(n, names(obj))]]
      }
    }
  }
  obj <- obj$results
  if(is.numeric(nets)){nets <- paste0("nets", 0:2)[nets]}
  nets <- match.arg(nets, paste0("nets", 0:2))
  type <- as.numeric(gsub("nets", "", nets))
  type <- switch(type + 1, "", "(thresholded, OR rule)", "(thresholded, AND rule)")
  obj <- obj[[nets]]
  ind <- match.arg(ind, colnames(obj)[-(1:2)])
  ff <- as.formula(paste0(ind, " ~ net"))
  boxplot(ff, obj, main = paste(ind, type), ylim = ylim, ...)
}

##### simFits: fit a series of mlGVAR models (OLD)
simFits <- function(data, crits, m = NULL, methods = NULL, ...){
  if(any(crits == "all")){
    crits0 <- c("aic", "bic", "ebic", "cv", "cv", "cp", "adjr2")
    crits <- c(crits0, crits[-1])
  }
  crits <- match.arg(tolower(crits), c(
    "aic", "bic", "ebic", "cv", "cp", "adjr2", "r2", "rss"), 
    several.ok = TRUE)
  lam <- rep("min", length(crits))
  seeds <- rep(1, length(crits))
  if(any(crits[duplicated(crits)] == "cv")){
    lam[which(crits == "cv")[2]] <- "1se"
    seeds[which(crits == "cv")] <- 2
  }
  if(is.null(methods)){
    methods <- rep("glmnet", length(crits))
  } else if(length(methods) == 1){
    methods <- rep(methods, length(crits))
  }
  methods[grepl("p|r", crits)] <- "subset"
  x <- length(crits) + 1
  if(!is.null(m)){
    if(!is(m, "list")){m <- rep(list(m), x)}
    if(any(!sapply(m, is.null))){
      methods[which(!sapply(m, is.null)) - 1] <- "glinternet"
    }
  } else {
    m <- rep(list(NULL), x)
  }
  fits0 <- lapply(seq_len(x), function(i){
    if(i == 1){
      pb <- txtProgressBar(min = 0, max = x, style = 3)
      assign("pb", pb, envir = parent.frame(2))
      xx <- mlGVAR(data, m = m[[i]], subjectNets = FALSE, verbose = FALSE, ...)
      setTxtProgressBar(pb, i)
    } else {
      xx <- suppressWarnings(
        mlGVAR(data, m = m[[i]], subjectNets = FALSE, selectFUN = TRUE, 
               verbose = FALSE, criterion = crits[i - 1], which.lam = lam[i - 1], 
               varSeed = switch(seeds[i - 1], NULL, 666), 
               method = methods[i - 1], ...))
      setTxtProgressBar(pb, i)
    }
    return(xx)
  })
  fnames <- paste0("fit", 0:(x - 1))
  names(fits0) <- fnames
  inds <- data.frame(net = fnames, crit = c(0, crits), 
                     lam = c(0, lam), method = c(0, methods))
  if(!is.null(m)){
    m0 <- unlist(lapply(m, function(mm) paste(colnames(data)[mm], collapse = "|")))
    if(any(m0 == "")){m0[m0 == ""] <- 0}
    inds$mods <- m0
  }
  y <- as.data.frame(SURll(fits0))
  y$net <- rownames(y)
  inds0 <- merge(y, inds, by = "net", sort = FALSE)
  fits1 <- fits0[which(!duplicated(inds0[, "LL"]))]
  inds1 <- inds0[which(!duplicated(inds0[, "LL"])), ]
  out <- list(fits = fits1, inds = inds1, fits0 = fits0, inds0 = inds0)
  list2env(out, .GlobalEnv)
  return(message("\nfits, fits0, inds, and inds0 in .GlobalEnv"))
}

##### rcheck:
rcheck <- function(x, quant = .9){
  x1 <- range(x)
  x2 <- range(abs(setdiff(x, 0)))
  x3 <- unname(quantile(x, probs = c(1 - quant[1], quant[1])))
  structure(rbind.data.frame(x1, x2, x3), names = c('lower', 'upper'), 
            row.names = c('full', 'nonZero', paste0('q', quant * 100)))
}

##### genPm2
genPm2 <- function(p, m2, niter = 100){
  pm2 <- cbind.data.frame(p = rep(p, each = length(m2)), 
                          m2 = ceiling(apply(expand.grid(m2, (p * (p - 1)/2)), 1, prod)))
  pm3 <- cbind(iter = rep(1:niter, nrow(pm2)), p = rep(pm2$p, each = niter), 
               m2 = rep(pm2$m2, each = niter))
  return(pm3)
}


processFits <- function(x, fits, threshold = TRUE, rule = 'and', avg = TRUE, 
                        inds = 'all', mprop = FALSE, multi = FALSE, saveIter = TRUE){
  if(multi){
    call <- replace(as.list(match.call())[-1], c('multi', 'mprop', 'saveIter'), 
                    list(multi = FALSE, mprop = TRUE, saveIter = FALSE))
    out <- lapply(seq_along(x), function(i){
      args <- replace(call, c('x', 'fits'), list(x = x[[i]], fits = fits[[i]]))
      do.call(processFits, args)
    })
    out1 <- do.call(rbind.data.frame, lapply(out, '[[', 'out1'))
    out1$m2 <- factor(out1$m2, levels = sort(unique(out1$m2)))
    out2 <- do.call(rbind.data.frame, lapply(out, '[[', 'out2'))
    out2$m2 <- factor(out2$m2, levels = sort(unique(out2$m2)))
    return(list(out1 = out1, out2 = out2))
  }
  n <- length(x)
  N <- nrow(x[[1]]$data)
  xn <- setNames(data.frame(apply(do.call(rbind, strsplit(names(x), '_')), 2, function(z){
    as.numeric(gsub('[^0-9]', '', z))
  })), c('iter', 'p', 'm2'))
  if(mprop){xn$m2 <- round(xn$m2/(xn$p * (xn$p - 1)/2), 1)}
  real1 <- lapply(x, net)
  real2 <- lapply(x, netInts)
  est1 <- lapply(fits, net, threshold = threshold, rule = rule)
  est2 <- lapply(fits, netInts, threshold = threshold, rule = rule, avg = avg, empty = FALSE)
  allinds <- c('cosine', 'mse', 'mae', 'correlation')
  if('all' %in% inds){inds <- allinds}
  inds <- match.arg(tolower(inds), allinds, several.ok = TRUE)
  if('mse' %in% inds){inds <- gsub('mse', 'MSE', inds)}
  if('mae' %in% inds){inds <- gsub('mae', 'MAE', inds)}
  inds <- Hmisc::capitalize(inds)
  out1 <- setNames(data.frame(do.call(cbind, lapply(inds, function(ind) sapply(seq_along(real1), function(i){
    matrixDist(real1[[i]], est1[[i]], ind = ind, directed = FALSE)
  })))), inds)
  if(any(is.na(out1))){out1[is.na(out1)] <- 0}
  out2 <- setNames(data.frame(do.call(cbind, lapply(inds, function(ind) sapply(seq_along(real1), function(i){
    matrixDist(real2[[i]], est2[[i]], ind = ind, directed = !avg)
  })))), inds)
  if(any(is.na(out2))){out2[is.na(out2)] <- 0}
  out3 <- do.call(rbind, lapply(1:n, function(z) performance(est1[[z]], real1[[z]], inds = 'between')))
  out4 <- do.call(rbind, lapply(1:n, function(z) performance(est2[[z]], real2[[z]], inds = ifelse(avg, 'between', 'beta'))))
  out3 <- structure(data.frame(out3), names = Hmisc::capitalize(colnames(out3)), row.names = 1:nrow(out3))
  out4 <- structure(data.frame(out4), names = Hmisc::capitalize(colnames(out4)), row.names = 1:nrow(out4))
  colnames(out3)[5] <- colnames(out4)[5] <- 'MCC'
  out3 <- out3[,-4]; out4 <- out4[,-4]
  inds1 <- rep(factor(rep(inds, each = n), levels = inds), 2)
  inds2 <- rep(factor(rep(colnames(out3), each = n)), 2)
  out1 <- unname(unlist(out1))
  out2 <- unname(unlist(out2))
  out3 <- unname(unlist(out3))
  out4 <- unname(unlist(out4))
  type <- factor(rep(c('Pairwise', 'Interactions'), c(length(out1), length(out2))), levels = c('Pairwise', 'Interactions'))
  out01 <- cbind.data.frame(value = c(out1, out2), measure = inds1, type = type, rbind(xn, xn), N = N)
  out02 <- cbind.data.frame(value = c(out3, out4), measure = inds2, type = type, rbind(xn, xn), N = N)
  if(!saveIter){
    out01 <- out01[, setdiff(colnames(out01), 'iter')]
    out02 <- out02[, setdiff(colnames(out02), 'iter')]
  }
  return(list(out1 = out01, out2 = out02))
}

getAllOuts <- function(path = '~/Desktop/dataSim/results'){
  setwd(path)
  dirs <- dir()[grep('[0-9]', dir())]
  dats <- appd(setNames(lapply(dirs, function(z){
    setNames(lapply(paste0(z, '/', dir(z)[grep('dat', dir(z))]), readRDS), 
             gsub('.*/|.RDS', '', paste0(z, '/', dir(z)[grep('dat', dir(z))])))
  }), paste0(dirs, '_dat')))
  fits <- appd(setNames(lapply(dirs, function(z){
    setNames(lapply(paste0(z, '/', dir(z)[grep('fit', dir(z))]), readRDS), 
             gsub('.*/|.RDS', '', paste0(z, '/', dir(z)[grep('fit', dir(z))])))
  }), paste0(dirs, '_fit')))
  list2env(list(fits = fits, dats = dats), .GlobalEnv)
}

### ======================================================================== ###
### ============================ SIM RESULTS =============================== ###
### ======================================================================== ###
##### compareMods: compare simulated networks with originals
compareMods <- function(fits, Sim, threshold = FALSE, Res = NULL, ind = "correlation", 
                        gamma = .5, mats = NULL, rule = "OR"){
  ind <- match.arg(tolower(ind), c(
    "cosine", "mse", "rmse", "ssd", "mae", "msd", "correlation"))
  atts <- names(attributes(fits))
  if("between" %in% names(Sim)){
    if(any(diag(Sim$between) == 1)){diag(Sim$between) <- 0}
    if(any(c("mlGVAR", "lmerVAR") %in% atts)){fits <- list(fits)}
    if(is.null(mats)){mats <- c("kappa", "pcc", "beta", "pdc", "between")}
    out <- sapply(fits, function(z){
      K <- lapply(mats, function(k){
        k0 <- net(z, k, threshold, rule)
        if(k == "pdc"){k0 <- t(k0)}
        return(k0)
      })
      out2 <- sapply(seq_along(K), function(i){
        matrixDist(K[[i]], Sim[[i + 1]], ind = ind, 
                   directed = isTRUE(mats[i] %in% c("beta", "pdc")))
      })
      names(out2) <- names(Sim)[2:6]
      return(out2)
    })
    if(isTRUE(Res)){
      Res <- graphicalVAR::mlGraphicalVAR(
        Sim$data, Sim$vars, idvar = Sim$idvar, 
        subjectNetworks = FALSE, gamma = gamma)
      assign("Res", Res, envir = .GlobalEnv)
    }
    if(is(Res, "mlGraphicalVAR")){
      r1 <- matrixDist(Res$fixedResults$kappa, Sim$fixedKappa, ind = ind)
      r2 <- matrixDist(Res$fixedPCC, Sim$fixedPCC, ind = ind)
      r3 <- matrixDist(Res$fixedResults$beta[, -1], Sim$fixedBeta, ind = ind)
      r4 <- matrixDist(Res$fixedPDC, Sim$fixedPDC, ind = ind)
      r5 <- matrixDist(Res$betweenNet, Sim$between, ind = ind)
      out <- cbind(out, Res = c(r1, r2, r3, r4, r5))
    }
    out <- t(out)
    if(nrow(out) <= 2){
      if(length(fits) == 2){
        rownames(out) <- if(!is.null(names(fits))){names(fits)} else {c("fit1", "fit2")}
      } else {
        rownames(out) <- switch(nrow(out), "fit", c("fit", "Res"))
      }
    }
  } else if(is(Sim, "mlVARsim")){
    if(any(c("mlGVAR", "lmerVAR") %in% atts)){fits <- list(fits)}
    x <- c("temporal", "contemporaneous", "between")
    Sim <- lapply(x, function(z){
      #z1 <- t(as.matrix(mlVAR::getNet(Sim, z)))
      z1 <- net(Sim, z)
      dimnames(z1) <- NULL
      return(z1)
    })
    out <- sapply(fits, function(z){
      nets <- lapply(x, function(n){
        X <- switch(2 - is(z, "mlVAR"), mlVAR::getNet, net)
        return(X(z, n, threshold = threshold, rule = rule))
      })
      out2 <- setNames(mapply(
        matrixDist, nets, Sim, ind = ind, 
        directed = c(TRUE, FALSE, FALSE)
      ), x)
      return(out2)
    })
    out <- t(out)
  } else {
    mats <- c("kappa", "beta", "PCC", "PDC")
    nets <- lapply(mats, function(z){
      z0 <- net(fits, z, threshold, rule)
      if(z == "PDC"){z0 <- t(z0)}
      return(z0)
    })
    names(nets) <- mats
    if(length(Sim) > length(mats)){
      Sim <- list(kappa = Sim$fixedKappa, beta = Sim$fixedBeta, 
                  PCC = Sim$fixedPCC, PDC = Sim$fixedPDC)
    }
    out <- mapply(matrixDist, nets, Sim, ind = ind, 
                  directed = rep(c(FALSE, TRUE), 2))
    if(is(Res, "graphicalVAR")){
      Res1 <- Res[mats]
      Res1$beta <- Res1$beta[, -1]
      out1 <- mapply(matrixDist, Res1, Sim, ind = ind, 
                     directed = rep(c(FALSE, TRUE), 2))
      out <- rbind(out, out1)
      rownames(out) <- c("fit", "Res")
    }
  }
  output <- list(out = out)
  if(ifelse(is(out, "matrix"), ifelse(nrow(out) > 1, TRUE, FALSE), FALSE)){
    rnks <- apply(out, 2, function(z){
      order(order(z, decreasing = ifelse(startsWith(ind, "c"), TRUE, FALSE)))
    })
    rownames(rnks) <- rownames(out)
    output$ranks <- rnks
  }
  names(output)[1] <- ifelse(ind == "correlation", "cor", ind)
  return(output)
}

##### performance: sensitivity, specificity, precision, accuracy
performance <- function(est, trueMod, threshold = FALSE, combine = FALSE, 
                        inds = "all", rule = "OR", getVals = FALSE, 
                        mcc = TRUE, rmNAs = TRUE){
  if(combine){
    nn <- switch(2 - is.null(names(est)), paste0("fit", 1:length(est)), names(est))
    xx <- lapply(lapply(est, performance, trueMod, threshold, 
                        inds = inds, rule = rule), function(z) data.frame(t(z)))
    nn2 <- colnames(xx[[1]])
    xx <- lapply(seq_along(nn2), function(z){
      data.frame(t(do.call(cbind, lapply(xx, '[', z))), row.names = nn)})
    names(xx) <- nn2
    return(xx)
  }
  if(all(inds == "all")){
    inds <- c("kappa", "beta", "PCC", "PDC")
    atts <- names(attributes(est))
    if(any(c("mlGVAR", "lmerVAR") %in% atts) | is(est, "mlGraphicalVAR")){
      inds <- c(inds, "between")
    }
    if(is(est, "matrix") & is(trueMod, "matrix")){inds <- "beta"}
  }
  inds <- match.arg(tolower(inds), c('kappa', 'beta', 'pcc', 'pdc', 'between'), several.ok = TRUE)
  if(any(grepl('^p', inds))){inds[grepl('^p', inds)] <- toupper(inds[grepl('^p', inds)])}
  tt <- FALSE
  if(length(inds) > 1){
    if(isTRUE(attr(trueMod, "simMLgvar")) | is(trueMod, "simMLgvar")){
      trueMod <- setNames(lapply(inds, function(z){
        n2 <- as.matrix(net(trueMod, z))
        if(z == "between"){diag(n2) <- 0}
        dimnames(n2) <- NULL
        return(n2)
      }), inds)
      #tt <- !is(est, "mlGraphicalVAR")
    } else if(is(trueMod, "mlVARsim")){
      inds <- c("temporal", "contemporaneous", "between")
      trueMod <- setNames(lapply(inds, function(z){
        #n2 <- t(as.matrix(mlVAR::getNet(trueMod, z)))
        n2 <- net(trueMod, z)
        dimnames(n2) <- NULL
        return(n2)
      }), inds)
    }
    nets1 <- setNames(lapply(inds, function(z){
      n1 <- as.matrix(net(est, z, threshold, rule))
      if(tt & z == "PDC"){n1 <- t(n1)}
      dimnames(n1) <- NULL
      return(n1)
    }), inds)
  } else {
    if(!is(est, "matrix")){est <- net(est, inds)}
    nets1 <- setNames(list(est), inds)
    trueMod <- setNames(list(trueMod), inds)
  }
  p <- unique(sapply(trueMod, nrow))
  trueVals <- lapply(inds, function(z){
    if(z %in% c("PCC", "between", "contemporaneous", "kappa")){
      z1 <- trueMod[[z]][lower.tri(trueMod[[z]])]
    } else {
      z1 <- as.vector(trueMod[[z]])
    }
    return(list(truePos = which(z1 != 0), trueNeg = which(z1 == 0)))
  })
  estVals <- lapply(inds, function(z){
    if(z %in% c("PCC", "between", "contemporaneous", "kappa")){
      z1 <- nets1[[z]][lower.tri(nets1[[z]])]
    } else {
      z1 <- as.vector(nets1[[z]])
    }
    return(list(estPos = which(z1 != 0), estNeg = which(z1 == 0)))
  })
  tp <- mapply(function(x1, x2){
    length(intersect(x1[[1]], x2[[1]]))}, estVals, trueVals)
  tn <- mapply(function(x1, x2){
    length(intersect(x1[[2]], x2[[2]]))}, estVals, trueVals)
  fp <- sapply(lapply(estVals, '[[', 1), length) - tp
  fn <- sapply(lapply(estVals, '[[', 2), length) - tn
  sensitivity <- tp/(tp + fn)
  specificity <- tn/(tn + fp)
  precision <- tp/(tp + fp)
  accuracy <- (tp + tn)/(tp + fp + tn + fn)
  out <- cbind.data.frame(sensitivity, specificity, precision, accuracy)
  if(mcc){
    numerator <- ((tp * tn) - (fp * fn))
    denom <- (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    mcc <- numerator/ifelse(denom == 0, 1, sqrt(denom))
    out <- cbind.data.frame(out, mcc = mcc)
  }
  if(any(is.na(out)) & rmNAs){out[is.na(out)] <- 0}
  rownames(out) <- inds
  if(getVals){
    out2 <- cbind.data.frame(tp, tn, fp, fn)
    rownames(out2) <- inds
    out <- list(performance = out, indices = out2, 
                estVals = estVals, trueVals = trueVals)
  }
  if(tt){out <- t(out)}
  return(out)
}

##### getNets
getNets <- function(obj, interactions = TRUE, dat = NULL){
  if(!is.null(dat)){obj <- list(call = obj$call, dat = dat, fit = obj)}
  x <- switch(2 - all(c("dat", "fit") %in% names(obj)), 1, length(obj))
  if(x == 1){obj <- list(obj)}
  allNets <- c("temporal", "contemporaneous", "between")
  calls <- simnets <- fitnets <- list()
  for(i in seq_len(x)){
    out <- obj[[i]]
    calls[[i]] <- switch(2 - ("call" %in% names(out)), out$call, list())
    simnets[[i]] <- setNames(vector("list", 1), "nets0")
    simnets[[i]][[1]] <- setNames(lapply(allNets, function(z) net(out$dat, z)), allNets)
    if(is.null(out$call$m)){interactions <- FALSE}
    if(interactions){simnets[[i]][[1]]$interactions <- netInts(out$dat)}
    fitnets[[i]] <- setNames(vector("list", 3), paste0("nets", 0:2))
    for(j in 1:3){
      thresh <- c(FALSE, TRUE, TRUE); rule <- c("or", "or", "and")
      fitnets[[i]][[j]] <- setNames(lapply(allNets, function(z) net(out$fit, z, thresh[j], rule[j])), allNets)
      if(interactions){fitnets[[i]][[j]]$interactions <- netInts(out$fit, threshold = thresh[j])}
    }
  }
  if(!is.null(names(obj))){names(calls) <- names(simnets) <- names(fitnets) <- names(obj)}
  #if(x == 1){calls <- calls[[1]]; simnets <- simnets[[1]]; fitnets <- fitnets[[1]]}
  output <- list(calls = calls, simnets = simnets, fitnets = fitnets)
  return(output)
}

##### results
results <- function(obj, nets = 1, ind = "perf", X = NULL, abbv = 10){
  nn <- c("temporal", "contemporaneous", "between", "interactions")
  if(any(sapply(lapply(obj$calls, '[[', "m"), is.null))){nn <- nn[1:3]}
  if(isTRUE(attr(obj, "SURnet"))){nn <- setdiff(nn, "between")}
  if(ind %in% c("cosine", "cor", "mse", "mae")){
    Y1 <- lapply(1:3, function(z){
      yy <- do.call(rbind, lapply(seq_along(obj$calls), function(i){
        sapply(seq_along(nn), function(j){matrixDist(
          obj$fitnets[[i]][[z]][[j]], obj$simnets[[i]]$nets0[[j]], 
          ind, isTRUE(nn[j] %in% c("temporal", "interactions")))})
      }))
      colnames(yy) <- nn
      if(!is.null(names(obj$calls))){
        rownames(yy) <- abbreviate(names(obj$calls), minlength = abbv)}
      return(yy)
    })
  } else {
    Y1 <- lapply(1:3, function(z){
      setNames(lapply(seq_along(nn), function(j){
        yy <- do.call(rbind, lapply(seq_along(obj$calls), function(i){
          performance(obj$fitnets[[i]][[z]][[j]], obj$simnets[[i]]$nets0[[j]])}))
        if(!is.null(names(obj$calls))){
          rownames(yy) <- abbreviate(names(obj$calls), minlength = abbv)}
        return(yy)
      }), nn)
    })
  }
  if(is.null(X)){
    return(Y1[[nets]])
  } else {
    k <- gsub("_.*", "", names(obj$calls)[1])
    X <- subset(X, FUN == k)
    if(ind %in% c("cosine", "cor", "mse", "mae")){
      Y1 <- cbind(X, Y1[[nets]])
    } else {
      Y1 <- setNames(lapply(seq_along(nn), function(z) cbind(X, Y1[[nets]][[z]])), nn)
    }
    return(Y1)
  }
}


### ======================================================================== ###
### ======================= EXTRACT CLUSTER RESULTS ======================== ###
### ======================================================================== ###
##### GET RESULTS
getData <- function(fits, data, ind = "cosine", threshold = FALSE, ints = TRUE, 
                    exclude = "PDC", NAs = FALSE, rule = "OR"){
  if(ind %in% c("perf", "performance")){
    mets <- c("sens", "spec", "prec", "acc")
    out <- lapply(seq_along(fits), function(i){
      data.frame(t(performance(fits[[i]], data[[i]], threshold, rule = rule)))})
    out2 <- setNames(lapply(1:4, function(i){
      i <- do.call(rbind, lapply(out, '[[', i))
      colnames(i) <- mets
      return(i)
    }), c("kappa", "beta", "PCC", "PDC"))
    if(!is.null(exclude)){if(all(exclude != "none")){out2 <- out2[-which(names(out2) %in% exclude)]}}
    if(ints){
      out3 <- do.call(rbind, lapply(seq_along(fits), function(i){
        int <- tryCatch({netInts(fits[[i]], threshold = threshold)}, 
                        error = function(e){diag(0, ncol(net(fits[[i]])))})
        return(performance(int, data[[i]]$mm$mb2, threshold, rule = rule))
      }))
      out2$interactions <- as.matrix(out3)
    }
    out <- do.call(rbind, lapply(out2, function(z) data.frame(value = c(z), ind = rep(mets, each = length(fits)))))
    out$net <- factor(rep(names(out2), each = nrow(out)/length(names(out2))))
    out$ind <- factor(out$ind)
    rownames(out) <- 1:nrow(out)
  } else {
    out <- t(sapply(seq_along(fits), function(i){
      compareMods(fits[[i]], data[[i]], threshold = threshold, ind = ind, rule = rule)[[1]]}))
    out <- data.frame(x = c(out), net = rep(c("kappa", "beta", "PCC", "PDC"), each = length(fits)))
    if(!is.null(exclude)){if(all(exclude != "none")){out <- out[-which(out$net %in% exclude), ]}}
    if(ints){
      out2 <- sapply(seq_along(fits), function(i){
        int <- tryCatch({netInts(fits[[i]], threshold = threshold)}, 
                        error = function(e){diag(0, ncol(net(fits[[i]])))})
        return(matrixDist(int, data[[i]]$mm$mb2, ind = ind))
      })
      out <- rbind(out, data.frame(x = out2, net = rep("interactions", length(fits))))
    }
    names(out)[1] <- "value"
  }
  if(!NAs){if(any(is.na(out))){out[is.na(out[, 1]), 1] <- 0}}
  attr(out, "type") <- ind
  return(out)
}

##### GET DESCRIPTIVES
getMeans <- function(x, se = TRUE){
  nn <- as.character(unique(x$net))
  if("ind" %in% names(x)){
    met <- as.character(unique(x$ind))
    out1 <- do.call(rbind, lapply(seq_along(nn), function(i){
      sapply(seq_along(met), function(j){
        mean(subset(x, net == nn[i] & ind == met[j])[, 1], na.rm = TRUE)})}))
    dimnames(out1) <- list(nn, paste0(met, ".M"))
    out2 <- do.call(rbind, lapply(seq_along(nn), function(i){
      sapply(seq_along(met), function(j){
        j <- subset(x, net == nn[i] & ind == met[j])
        if(any(is.na(j[, 1]))){j <- j[!is.na(j[, 1]), ]}
        jj <- sd(j[, 1], na.rm = TRUE)
        if(se){jj <- jj/sqrt(nrow(j))}
        return(jj)
      })}))
    dimnames(out2) <- list(nn, paste0(met, ifelse(se, ".SE", ".SD")))
    out <- cbind(out1, out2)
  } else {
    ms <- sapply(seq_along(nn), function(i) mean(subset(x, net == nn[i])[, 1], na.rm = TRUE))
    ses <- sapply(seq_along(nn), function(i){
      i <- subset(x, net == nn[i])
      if(any(is.na(i[, 1]))){i <- i[!is.na(i[, 1]), ]}
      ii <- sd(i[, 1], na.rm = TRUE)
      if(se){ii <- ii/sqrt(nrow(i))}
      return(ii)
    })
    out <- rbind(ms, ses)
    dimnames(out) <- list(c("mean", ifelse(se, "SE", "SD")), nn)
    if("type" %in% names(attributes(x))){
      rownames(out) <- paste0(attr(x, "type"), c(".M", ifelse(se, ".SE", ".SD")))
    }
  }
  return(out)
}

##### PLOT RESULTS
plotData <- function(x, breaks = NULL){
  require(ggplot2)
  if("ind" %in% names(x)){
    g <- ggplot(x, aes(x = ind, y = value)) + facet_grid(~ net) + xlab("Metric")
  } else {
    g <- ggplot(x, aes(x = net, y = value)) + xlab("Network")
  }
  g <- g + geom_boxplot(outlier.size = .3) + theme_bw()
  if(!is.null(attr(x, "type"))){
    tt <- attr(x, "type")
    if(tt == "cor"){tt <- "Correlation"}
    if(tt == "perf"){tt <- "Performance"}
    if(tt == "cosine"){tt <- "Cosine"}
    if(tt == "mse"){tt <- "MSE"}
    if(tt == "mae"){tt <- "MAE"}
    g <- g + ylab(tt)
  }
  if(!is.null(breaks)){g <- g + scale_y_continuous(breaks = breaks)}
  return(g)
}

##### Extension of 'getData' for multiple folders
makeData <- function(x, ind = 'cosine', exclude = 'PDC', useFits = TRUE, 
                     ints = TRUE, alt = NULL, rule = 'OR'){
  if(!useFits){ints <- FALSE}
  if(!is.null(alt)){
    if(!grepl("^/", alt)){alt <- paste0("/", alt)}
    if(!grepl("[.]RDS$", alt)){alt <- paste0(alt, ".RDS")}
    useFits <- FALSE
  }
  dirs <- dir()[grep(paste0("^", x), dir())]
  stopifnot(length(dirs) > 0)
  dirs <- dirs[order(as.numeric(gsub(paste0(letters, collapse = "|"), "", dirs)))]
  pb <- txtProgressBar(max = length(dirs), style = 3)
  X <- setNames(lapply(seq_along(dirs), function(z){
    intsZ <- ints
    data <- readRDS(paste0(dirs[z], "/data.RDS"))
    if(intsZ & !"mm" %in% unique(lapply(data, names))[[1]]){intsZ <- FALSE}
    fits <- readRDS(paste0(dirs[z], ifelse(
      useFits, "/fits.RDS", ifelse(!is.null(alt), alt, "/vars.RDS"))))
    errs <- which(sapply(fits, length) == 0)
    if(any(errs)){fits <- fits[-errs]; data <- data[-errs]}
    if(length(fits) > 0){
      x1 <- getData(fits, data, ind, exclude = exclude, ints = intsZ, rule = rule)
      x2 <- getData(fits, data, ind, T, exclude = exclude, ints = intsZ, rule = rule)
      x3 <- getData(fits, data, "perf", exclude = exclude, ints = intsZ, rule = rule)
      x4 <- getData(fits, data, "perf", T, exclude = exclude, ints = intsZ, rule = rule)
      out <- list(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    } else {
      out <- list()
    }
    setTxtProgressBar(pb, z)
    if(z == length(dirs)){close(pb)}
    if(any(errs)){attr(out, 'errs') <- unname(errs)}
    return(out)
  }), dirs)
  if(any(sapply(X, function(z) 'errs' %in% names(attributes(z))))){
    attr(X, 'errs') <- lapply(X, attr, 'errs')
  }
  attr(X, "type") <- ind
  return(X)
}

##### Connected to 'makeData'
makePlots <- function(outs, x = 1, ms = NULL, measure = 'cosine', plot = TRUE, 
                      se = TRUE, critNames = NULL, box = FALSE, 
                      vertical = TRUE, legPos = 'right'){
  require(ggplot2)
  measure <- match.arg(tolower(measure), c('cosine', 'correlation', 'mae', 'mse'))
  if(sum(grepl(measure, names(outs))) == 1){outs <- outs[[grep(measure, names(outs))]]}
  if(!is.null(ms)){
    if(isTRUE(ms)){
      ms <- c(.05, .1, .25)
      outs <- appd(outs)
      if(!x %in% 1:2){stop('Not prepared to plot this many factors at once')}
    } else if(is.numeric(ms) & length(ms) == 1){
      outs <- outs[[switch(as.character(ms), '0.05' = 1, '0.1' = 2, '0.25' = 3, ms)]]
      ms <- NULL
    }
  }
  if("type" %in% names(attributes(outs))){outs <- list(outs); one <- TRUE} else {one <- FALSE}
  if("type" %in% names(attributes(outs[[1]]))){measure <- attr(outs[[1]], "type")}
  tt <- as.numeric(gsub(paste0(letters, collapse = "|"), "", names(outs[[1]])))
  Z <- lapply(outs, function(z) lapply(z, function(zz) lapply(zz, getMeans, se = se)))
  Z <- lapply(Z, function(z){
    if(any(sapply(z, length) == 0)){
      kk <- which(sapply(z, length) == 0)
      z <- z[-kk]
    }
    zout <- do.call(rbind, lapply(z, '[[', x))
    attr(zout, 'tt') <- switch(2 - (length(z) == length(tt)), tt, tt[-kk])
    return(zout)
  })
  if(any(grepl('[.]M$', rownames(Z[[1]])))){
    Z <- lapply(Z, function(z){
      tz <- attr(z, 'tt')
      data.frame(tt = rep(tz, ncol(z)), 
                 mean = c(z[seq(1, nrow(z), by = 2), ]),
                 se = c(z[seq(2, nrow(z), by = 2), ]), 
                 Network = rep(Hmisc::capitalize(colnames(z)), each = length(tz)))
    })
    if(is.null(ms) & length(unique(sapply(Z, nrow))) == 1 & (length(unique(gsub("[0-9]", "", c(sapply(outs, names))))) > 1 | !is.null(critNames))){
      dat <- data.frame(do.call(rbind, Z), crit = rep(
        toupper(unique(gsub("[0-9]", "", c(sapply(outs, names))))), 
        each = unique(sapply(Z, nrow))))
      if(!is.null(critNames)){dat$crit <- rep(critNames, each = unique(sapply(Z, nrow)))}
    } else if(!one){
      mm <- sapply(Z, nrow)
      xx <- rep(c('Moderator', 'No Moderator'), length(mm)/2)
      #xx <- c("Moderator", "No Moderator")[order(sapply(Z, nrow), decreasing = TRUE)]
      #mm <- sapply(Z, nrow)[order(sapply(Z, nrow), decreasing = TRUE)]
      #if(length(Z) > 2 & FALSE){
      #  if(length(unique(mm)) == 1){
      #    mm <- rep(2:1, length(Z)/2)
      #  } else {
      #    mm <- sapply(Z, nrow)
      #  }
      #  #stopifnot(all(sapply(Z, nrow)[seq(1, length(Z), by = 2)] > sapply(Z, nrow)[seq(2, length(Z), by = 2)]))
      #  xx[which(mm == max(unique(mm)))] <- "Moderator"
      #  xx[which(mm == min(unique(mm)))] <- "No Moderator"
      #  mm <- sapply(Z, nrow)
      #  crit <- c(unlist(sapply(seq_along(mm), function(z) rep(xx[z], mm[z]))))
      #} else {
      #  crit <- rep(xx, mm)
      #}
      dat <- data.frame(do.call(rbind, Z), crit = rep(xx, mm))
      if(!is.null(ms)){
        Ms0 <- apply(rbind(seq(1, length(mm), by = 2), seq(2, length(mm), by = 2)), 2, function(z) sum(mm[z]))
        dat$M <- rep(ms, Ms0)
      }
    } else if(one){
      dat <- Z[[1]]
    }
    g <- ggplot(dat, aes(x = tt, y = mean, colour = Network, fill = Network)) + theme_bw()
    if(!one & is.null(ms)){if(vertical){g <- g + facet_grid(~ crit)} else {g <- g + facet_grid(crit ~ .)}}
    if(!one & !is.null(ms)){if(vertical){g <- g + facet_grid(M ~ crit)} else {g <- g + facet_grid(crit ~ M)}}
    if(!box){
      g <- g + geom_line() + geom_point(size = 1) + 
        geom_ribbon(aes(ymin = mean - se, ymax = mean + se, colour = NULL), alpha = .3)
    } else { # WORKS, NOT PROPERLY
      g <- g + geom_boxplot(aes(colour = NULL), outlier.size = 1, alpha = .7)
    }
    if(measure %in% c("cor", "correlation")){measure <- "Correlation"}
    if(measure == "cosine"){measure <- "Cosine Similarity"}
    if(measure == "mse"){measure <- "Mean Squared Error (MSE)"}
    if(measure == "mae"){measure <- "Mean Absolute Error (MAE)"}
    g <- g + ylab(measure) + theme(legend.position = legPos)
  } else {
    Z <- lapply(Z, function(z){
      tz <- attr(z, 'tt')
      xx <- unique(rownames(z))
      data.frame(tt = rep(rep(tz, each = length(xx)), 4),
                 ind = rep(c("Sensitivity", "Specificity", "Precision", "Accuracy"), each = nrow(z)),
                 mean = c(z[, 1:4]), se = c(z[, 5:8]), Network = rep(Hmisc::capitalize(xx), 4))
    })
    if(length(unique(sapply(Z, nrow))) == 1 & (length(unique(gsub("[0-9]", "", c(sapply(outs, names))))) > 1 | !is.null(critNames))){
      dat <- data.frame(do.call(rbind, Z), crit = rep(
        toupper(unique(gsub("[0-9]", "", c(sapply(outs, names))))), 
        each = unique(sapply(Z, nrow))))
      if(!is.null(critNames)){dat$crit <- rep(critNames, each = unique(sapply(Z, nrow)))}
    } else if(!one){
      mm <- sapply(Z, nrow)
      xx <- rep(c('Moderator', 'No Moderator'), length(mm)/2)
      dat <- data.frame(do.call(rbind, Z), crit = rep(xx, mm))
    } else if(one){
      dat <- Z[[1]]
    }
    g <- ggplot(dat, aes(x = tt, y = mean, colour = Network, fill = Network)) + 
      theme_bw() + ylab("Average Performance")
    if(!one){
      g <- g + facet_grid(ind ~ crit)
    } else {
      if(vertical){g <- g + facet_grid(~ ind)} else {g <- g + facet_grid(ind ~ .)}
    }
    if(!box){
      g <- g + geom_line() + geom_point(size = 1) + 
        geom_ribbon(aes(ymin = mean - se, ymax = mean + se, colour = NULL), alpha = .3)
    } else { # WORKS, NOT PROPERLY
      g <- g + geom_boxplot(aes(colour = NULL), outlier.size = 1, alpha = .7)
    }
  }
  g <- g + xlab("Number of Time Points") + theme(legend.position = legPos)
  if(x %in% 1:2){attr(dat, 'measure') <- measure}
  if(plot){return(g)} else {return(dat)}
}

##### appd: take a list and append all sublists; or collect certain lists
appd <- function(x, recursive = FALSE){
  if(length(x) == 1){
    inds <- c('correlation', 'cosine', 'mse', 'mae')
    inds <- ifelse(is.character(x), list(inds), list(paste0(inds, x)))[[1]]
    inds <- paste(paste0('^', inds, '$'), collapse = '|')
    out <- nlist(inds, TRUE)
  } else {
    out <- list()
    for(i in seq_along(x)){out <- append(out, x[[i]])}
    if(recursive){
      while(all(sapply(out, is, 'list'))){
        out <- appd(out, recursive = FALSE)
      }
    }
  }
  return(out)
}


########## GATHERING DATA
fullData <- function(x, crit = NULL, ints = TRUE, full = FALSE, rmInts = FALSE){
  fullData2 <- function(x, rmInts = FALSE){
    out <- setNames(lapply(1:4, function(z){
      zz <- do.call(rbind, lapply(x, '[[', z))
      zz$net <- factor(zz$net)
      if(is.character(zz$crit)){zz$crit <- factor(zz$crit)}
      if(z %in% 3:4){zz$ind <- factor(zz$ind)}
      return(zz)
    }), paste0('x', 1:4))
    if(all(sapply(1:4, function(z) 'interactions' %in% out[[z]]$net))){
      if(all(sapply(1:4, function(z) 'M' %in% colnames(out[[z]])))){
        if(length(unique(sapply(out, function(z) length(unique(subset(z, net == 'interactions')$M))))) == 1){rmInts <- TRUE}
      }
      if(rmInts){
        out <- lapply(out, function(z){
          zz <- subset(z, net != 'interactions')
          zz$net <- factor(zz$net)
          return(zz)
        })
      }
    }
    return(out)
  }
  if(!is.null(crit)){
    if(all(crit == 'm')){
      crit <- c(.05, .1, .25)
      out <- lapply(1:3, function(z) fullData(x[[z]], crit[z], ints = ints, full = FALSE))
      if(full){out <- fullData2(out, rmInts = rmInts)}
      return(out)
    }
  }
  out <- setNames(lapply(1:4, function(z){
    do.call(rbind, lapply(seq_along(x), function(zz){
      nn <- names(x[[zz]])
      tt <- as.numeric(gsub(paste0(letters, collapse = "|"), "", nn))
      crit <- ifelse(is.null(crit), toupper(unique(gsub("[0-9]", "", nn))), crit)
      out <- lapply(x[[zz]], '[[', z)
      out <- cbind.data.frame(do.call(rbind, out), tt = rep(tt, sapply(out, nrow)), crit = crit)
      if(ints){out$M <- as.numeric(isTRUE("interactions" %in% unique(x[[zz]][[1]][[1]]$net)))}
      rownames(out) <- 1:nrow(out)
      return(out)
    }))
  }), paste0("x", 1:4))
  if(full){out <- fullData2(out, rmInts = rmInts)}
  return(out)
}
mms <- function(x, search = FALSE, double = FALSE, ind = 'net', y = NULL, ...){
  ff <- as.formula(paste0('value ~ .', ifelse(double, '^2', '')))
  if(!is.null(y)){x <- lapply(seq_along(x), function(z) x[[z]][x[[z]][[ind]] == y, -which(colnames(x[[z]]) == ind)])}
  if(!search){
    m <- lapply(seq_along(x), function(z) lm(ff, x[[z]]))
  } else {
    m <- lapply(seq_along(x), function(z) step(lm(ff, x[[z]]), ...))
  }
  return(m)
}
s2pool <- function(x, ind = NULL, vals = NULL, d = TRUE, all = FALSE){
  if(is(x, 'lm')){x <- x$model}
  if(is.null(vals)){vals <- c(0, 1)}
  if(is.null(ind)){ind <- 'M'}
  x1 <- x[x[, ind] == vals[1], 'value']
  x2 <- x[x[, ind] == vals[2], 'value']
  m1 <- mean(x1, na.rm = T)
  m2 <- mean(x2, na.rm = T)
  s1 <- var(x1, na.rm = T)
  s2 <- var(x2, na.rm = T)
  df1 <- length(x1) - 1
  df2 <- length(x2) - 1
  df0 <- df1 + df2
  out <- spool <- ((df1/df0) * s1) + ((df2/df0) * s2)
  if(d){out <- tryCatch({abs(m1 - m2)/sqrt(out)}, error = function(e){NA})}
  if(d){if(m1 == m2){out <- NA} else {names(out) <- paste0("M", 0:1)[which.max(c(m1, m2))]}}
  if(!all){return(out)} else {return(list(out = out, means = c(m1, m2), vars = c(s1, s2), s2pool = spool))}
}
ssx <- function(x, ind = 'net', v = 'PCC'){
  if(all(v == 'all') | isTRUE(v)){v <- unique(x[[1]][, ind])}
  if(is.factor(v)){v <- as.character(v)}
  if(length(v) == 1){
    out <- lapply(x, function(z) z[z[, ind] == v, -which(colnames(z) == ind)])
  } else {
    out <- lapply(seq_along(v), function(z) ssx(x, ind = ind, v = v[z]))
    names(out) <- paste0(ifelse(is.numeric(v), ind, ''), v)
    olen <- unique(sapply(out, length))
    stopifnot(length(olen) == 1)
    stopifnot(length(unique(lapply(out, names))) == 1)
    n2 <- unique(lapply(out, names))[[1]]
    out <- lapply(1:olen, function(z) lapply(out, '[[', z))
    names(out) <- n2
  }
  return(out)
}
s2poolBy <- function(x, ind = 'net'){
  ind <- match.arg(ind, setdiff(unique(unlist(lapply(x, colnames))), 'value'))
  if(ind == 'ind'){x <- x[which(sapply(lapply(x, colnames), function(z) 'ind' %in% z))]}
  xlen <- length(x)
  inds <- unique(x[[1]][, ind])
  if(is.factor(inds)){inds <- as.character(inds)}
  x <- setNames(lapply(inds, function(z){
    ssx(x, ind = ind, v = z)
  }), paste0(ifelse(is.numeric(inds), ind, ''), inds))
  stopifnot(all(unique(sapply(x, length)) == xlen))
  xn <- names(x[[1]])
  x <- setNames(lapply(seq_len(xlen), function(z) lapply(x, '[[', z)), xn)
  out <- setNames(lapply(seq_along(x), function(z){
    z1 <- sapply(x[[z]], s2pool, d = TRUE)
    z2 <- list()
    if(!all(is.na(z1))){
      z2 <- data.frame(do.call(rbind, strsplit(names(z1), '[.]')), d = z1, stringsAsFactors = FALSE)
      z3 <- any(sapply(names(z1), function(k) sum(strsplit(k, '')[[1]] == '.')) > 1)
      if(isTRUE(z3)){z2 <- z2[, -which(apply(z2, 2, function(z) any(grepl(ind, z))))]}
      rownames(z2) <- 1:nrow(z2)
      colnames(z2) <- c(ind, 'M', 'd')
      if(ind == 'crit'){z2$crit <- as.numeric(paste0('.', z2$crit))}
      z2$M <- as.numeric(gsub('M', '', z2$M))
      if(is.numeric(inds)){z2[, ind] <- as.numeric(gsub(paste0(letters, collapse = '|'), '', z2[, ind]))}
    }
    return(z2)
  }), names(x))
  if(any(sapply(out, length) == 0)){out <- out[-which(sapply(out, length) == 0)]}
  return(out)
}
s22 <- function(x, ind1 = 'net', ind2 = 'tt', neg = FALSE){
  inds <- unique(x[[1]][, ind1])
  if(is.factor(inds)){inds <- as.character(inds)}
  X <- lapply(inds, function(z) ssx(x, ind1, z))
  XX <- lapply(X, function(z) s2poolBy(z, ind = ind2))
  poop <- lapply(seq_along(XX), function(z) lapply(XX[[z]], function(zz) cbind(zz, ind = inds[z])))
  poop <- lapply(seq_len(unique(sapply(XX, length))), function(z) lapply(poop, '[[', z))
  poop <- lapply(poop, function(z) do.call(rbind, z))
  names(poop) <- names(XX[[1]])
  if(neg & all(sapply(poop, function(z) 'M' %in% colnames(z)))){
    if(all(sapply(poop, function(z) all(z$M %in% c(0, 1))))){
      poop <- lapply(poop, function(z){
        z[z$M == 0, 'd'] <- -z[z$M == 0, 'd']
        return(z)
      })
    }
  }
  return(poop)
}


### ======================================================================== ###
### ======================================================================== ###
##### getCommandArgs: only really works if submission scripts have been created by makeScripts
getCommandArgs <- function(x = commandArgs(trailingOnly = TRUE)){
  invisible(lapply(strsplit(x, '=', fixed = TRUE), function(a){
    if(length(a) == 2){
      assign(substr(a[1], 1, (nchar(a[1]) - 1)), switch(
        substr(a[1], nchar(a[1]), nchar(a[1])), 
        C = as.character, L = as.logical, 
        N = as.numeric)(a[2]), envir = .GlobalEnv)
    }
  }))
}

##### evalFits: summarize iterated fitting from cluster
evalFits <- function(out1, compare = FALSE, combine = TRUE, threshold = TRUE, 
                     avg = TRUE, getOriginal = FALSE, rmFailed = TRUE, 
                     multi = FALSE, getNets = FALSE){
  mconds <- isTRUE((combine %in% c(TRUE, 'all')) & isFALSE(getNets) & !getOriginal & rmFailed)
  if(!identical(multi, FALSE) & mconds){
    gg <- function(x, nums = TRUE){
      x <- gsub(ifelse(nums, '[^0-9]', '[0-9]'), '', x)
      if(nums){x <- as.numeric(x)}
      return(x)
    }
    if(any(grepl('_', names(out1)))){
      nn <- do.call(rbind, strsplit(names(out1), '_'))
      nn1 <- apply(nn, 2, gg, nums = TRUE)
      if(any(is.na(nn1))){
        nn0 <- which(apply(nn1, 2, function(z) any(is.na(z))))
        nn1 <- nn1[, -nn0, drop = FALSE]
        if(ncol(nn1) > 1){
          nn2 <- apply(nn[, -nn0], 2, function(z) toupper(unique(gg(z, FALSE))))
          stopifnot(any(nn2 == 'N'))
          N <- nn1[, which(nn2 == 'N')]
          if('P' %in% nn2){P <- nn1[, which(nn2 == 'P')]}
          if('M' %in% nn2){M <- nn1[, which(nn2 == 'M')]}
        } else {
          N <- as.vector(nn1)
        }
      } else {
        N <- as.vector(nn1)
      }
    } else {
      N <- switch(2 - isTRUE(multi), gg(names(out1), nums = TRUE), multi)
    }
    stopifnot(is.numeric(N) & length(N) == length(out1))
    call <- replace(as.list(match.call())[-1], 'multi', list(multi = FALSE))
    out <- lapply(seq_along(out1), function(i){
      ii <- do.call(evalFits, replace(call, 'out1', list(out1 = out1[[i]])))
      if(length(ii) == 0){stop(paste0('Error in list element ', i))}
      return(ii)
    })
    nn3 <- sapply(out, nrow)
    out <- cbind.data.frame(do.call(rbind, out), N = rep(N, nn3))
    if(exists('P', inherits = FALSE)){out <- cbind.data.frame(out, P = rep(P, nn3))}
    if(exists('M', inherits = FALSE)){out <- cbind.data.frame(out, M = rep(M/10, nn3))}
    out <- out[order(out$N), ]
    rownames(out) <- 1:nrow(out)
    fullorder1 <- c('Prune', 'noMod', 'covariate', 'AIC', 'BIC', 'EBIC25', 'EBIC50', 'CVmin', 'CV1se',
                    'splitAIC', 'bootstrapAIC', 'stabilityAIC')
    if(!'noMod' %in% levels(out$fit)){fullorder1 <- setdiff(fullorder1, 'noMod')}
    if(!'covariate' %in% levels(out$fit)){fullorder1 <- setdiff(fullorder1, 'covariate')}
    fullorder2 <- c(fullorder1, 'NCT', 'FGL')
    if(all(fullorder1 %in% levels(out$fit))){
      fullorder <- list(fullorder1, fullorder2)
      v <- which(sapply(fullorder, length) == length(levels(out$fit)))
      fullorder <- fullorder[[v]]
      if(!identical(levels(out$fit), fullorder)){
        out$fit <- as.character(out$fit)
        out$fit <- factor(out$fit, levels = fullorder)
      }
    }
    return(out)
  }
  stopifnot(!isTRUE(multi))
  if(is.character(combine)){
    combine <- match.arg(tolower(combine), c('none', 'rules', 'type', 'all'))
    combine <- switch(combine, none = FALSE, all = TRUE, 
                      type = ifelse(!identical(compare, FALSE), 'type', FALSE),
                      rules = ifelse(!identical(compare, FALSE), 'rules', TRUE))
  }
  if(is.null(names(out1))){names(out1) <- paste0('iter', 1:length(out1))}
  errs <- which(sapply(out1, function(z) any(sapply(z, class) == 'try-error')))
  if(length(unique(sapply(out1, length))) > 1 | length(errs) > 0){
    if(rmFailed){
      if(length(errs) > 0){
        for(i in seq_along(errs)){
          out1[[errs[i]]] <- out1[[errs[i]]][-which(sapply(out1[[errs[i]]], class) == 'try-error')]
        }
      }
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
  rules <- c('OR', 'AND')
  if(!identical(compare, FALSE) | !isFALSE(getNets)){
    if(!'criterion' %in% names(attributes(out1))){
      warning('"criterion" not provided as an attribute to the main list')
    }
    outs1 <- setNames(lapply(c('netInts', 'net'), function(FUN){
      setNames(lapply(rules, function(rule){
        args <- list(threshold = threshold, rule = rule)
        if(FUN == 'netInts'){args <- append(args, list(avg = avg, empty = FALSE))}
        lapply(out1, function(iter){
          lapply(iter[setdiff(names(iter), c('fitNCT', 'fitFGL'))], function(z){
            if(is.character(threshold) & 'call' %in% names(z)){
              args$threshold <- !any(c('criterion', 'which.lam') %in% names(z$call))
            } else {
              args$threshold <- threshold
            }
            do.call(match.fun(FUN), append(args, list(fit = z)))
          })
        })
      }), rules)
    }), c('interactions', 'pairwise'))
    if(!isFALSE(getNets)){
      if(is.character(getNets)){
        getNets <- match.arg(tolower(getNets), c('or', 'and', 'pairwise', 'interactions'), several.ok = TRUE)
        type <- c('pairwise', 'interactions')
        rules <- c('or', 'and')
        if(any(type %in% getNets)){
          outs1 <- outs1[[getNets[getNets %in% type]]]
          if(any(rules %in% getNets)){
            outs1 <- outs1[[toupper(getNets[getNets %in% rules])]]
          }
        } else if(any(rules %in% getNets)){
          rr <- toupper(getNets[getNets %in% rules])
          outs1 <- lapply(outs1, '[[', rr)
        }
      }
      return(outs1)
    }
    if(isTRUE(compare) | !is.character(compare)){
      compare0 <- c('Correlation', 'Cosine', 'MSE', 'MAE')
      compare <- ifelse(isTRUE(compare), 4, ifelse(all(compare %in% 3:4), compare[1], 4))
      compare <- switch(compare - 2, setdiff(compare0, 'Cosine'), compare0)
    } else {
      inds <- c('cosine', 'correlation', 'mse', 'rmse', 'ssd', 'mae', 'msd', 'performance')
      if('all' %in% compare){compare <- setdiff(inds, 'performance')}
      compare <- Hmisc::capitalize(match.arg(tolower(compare), inds, several.ok = TRUE))
      if(any(!startsWith(compare, 'C'))){
        compare[!startsWith(compare, 'C')] <- toupper(compare[!startsWith(compare, 'C')])
      }
    }
    if('PERFORMANCE' %in% compare){
      compare <- c('Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'MCC')
      outs <- list(lapply(outs1$pairwise, function(z){
        nn <- switch(2 - !is.null(names(z)), names(z), seq_along(z))
        setNames(do.call(rbind.data.frame, lapply(nn, function(xn){
          fits <- z[[xn]]
          names(fits) <- c('x', 'Prune', attr(out1, 'criterion'))
          if(is.character(xn)){xn <- as.numeric(gsub('[^0-9]', '', xn))}
          structure(data.frame(do.call(rbind.data.frame, lapply(
            2:length(fits), function(i){
              performance(fits[[i]], fits[[1]], inds = 'between')
            })), fit = factor(names(fits)[-1], levels = names(fits)[-1]), iter = xn), 
            row.names = 1:(length(fits) - 1))
        })), c(compare, 'fit', 'iter'))
      }))
    } else {
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
    }
    if(!identical(combine, FALSE)){
      outs <- lapply(outs, function(dat){
        d1 <- cbind.data.frame(do.call(rbind, dat), rule = factor(
          rep(names(dat), each = nrow(dat[[1]])), levels = names(dat)))
        data.frame(fit = rep(d1$fit, length(compare)), rule = rep(d1$rule, length(compare)),
                   measure = factor(rep(compare, each = nrow(d1)), levels = compare),
                   value = unlist(c(d1[, compare]), use.names = FALSE), 
                   iter = rep(d1$iter, length(compare)))
      })
      if(length(outs) == 1){outs <- outs[[1]]}
      if((isTRUE(combine) | combine == 'type') & is(outs, 'list')){
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
    if(getOriginal){outs <- list(out1 = out1, out2 = outs)}
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
      if(identical(w2, 'min')){
        message('That weird thing with criterion == min?')
        w2 <- 'CV'
      }
      if(identical(w2, 'CV')){
        w2 <- ifelse('which.lam' %in% names(mod$call), paste0(w2, mod$call$which.lam), 'CVmin')
      } 
      if(identical(w2, 'EBIC')){
        w2 <- ifelse('gamma' %in% names(mod$call), paste0(w2, mod$call$gamma * 100), 'EBIC50')
      }
      thresh <- ifelse(!is.character(threshold), threshold, isTRUE(w2 == 'Prune'))
      est1 <- netInts(mod, rule = xx$rule[i], avg = avg, threshold = thresh, empty = FALSE)
      est2 <- structure(performance(est1, x, inds = ifelse(avg, 'between', 'beta')), criterion = w2)
      return(est2)
    }), apply(xx, 1, paste, collapse = '_'))
    xxx <- append(f12, list(
      f3 = structure(z$fitNCT$performance, criterion = 'NCT'), 
      f4 = structure(z$fitFGL$performance, criterion = 'FGL')))
    return(xxx)
  })}, error = function(e){list()})
  out3 <- tryCatch({lapply(seq_len(unique(sapply(out2, length))), function(z){
    z1 <- data.frame(do.call(rbind, lapply(out2, '[[', z)))
    nn <- c('Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'MCC')
    rn <- as.numeric(gsub('[^0-9]', '', rownames(z1)))
    z2 <- unlist(c(z1), use.names = FALSE)
    z3 <- cbind.data.frame(measure = factor(rep(nn, each = nrow(z1)), levels = nn), 
                           value = z2, iter = rep(rn, ncol(z1)))
    rownames(z3) <- 1:nrow(z3)
    return(z3)
  })}, error = function(e){list()})
  out4 <- tryCatch({setNames(lapply(seq_along(rules), function(z){
    n <- length(out3)
    ind <- c(seq(z, (n - 2), length(rules)), (n - 1):n)
    nn <- unique(unique(lapply(out2, function(zz) unname(sapply(zz, attr, 'criterion'))))[[1]])
    z1 <- do.call(rbind.data.frame, out3[ind])
    z1 <- cbind.data.frame(fit = factor(rep(nn, each = nrow(out3[[1]])), levels = nn), z1)
    attr(z1, 'box') <- TRUE
    return(z1)
  }), rules)}, error = function(e){list()})
  if(combine & length(out4) > 0){
    out4 <- structure(cbind.data.frame(do.call(rbind, out4), rule = factor(
      rep(names(out4), each = nrow(out4[[1]])), levels = names(out4))), 
      row.names = 1:sum(sapply(out4, nrow)))
    out4 <- data.frame(fit = out4$fit, rule = out4$rule, 
                       out4[, setdiff(colnames(out4), c('fit', 'rule'))])
    attr(out4, 'box') <- TRUE
  }
  if(getOriginal){out4 <- list(out1 = out1, out2 = out4)}
  return(out4)
}

##### plotit: relevant to 'template' procedures
plotit <- function(dat, rule = 'or', err = 'se', legend = TRUE, plot = TRUE, 
                   ylab = FALSE, title = NULL, xlab = FALSE, ticks = FALSE, 
                   facet = 'measure', include = NULL, exclude = NULL){
  require(ggplot2)
  rules <- c('or', 'and')
  if(is.numeric(rule)){rule <- rules[rule]}
  if(is(dat, 'list')){
    if(length(facet) == 1 & is.character(facet)){
      dat <- dat[[match.arg(tolower(rule), rules)]]
    } else {
      atts <- setdiff(unique(sapply(dat, function(z) names(attributes(z)))), c('names', 'class', 'row.names'))
      dat <- structure(cbind.data.frame(do.call(rbind, dat), rule = factor(
        rep(names(dat), each = nrow(dat[[1]])), levels = names(dat))), 
        row.names = 1:sum(sapply(dat, nrow)))
      attributes(dat)[atts] <- TRUE
    }
  }
  atts <- names(attributes(dat))
  dat <- switch(2 - !is.null(exclude), subset(dat, !measure %in% Hmisc::capitalize(exclude)), switch(
    2 - !is.null(include), subset(dat, measure %in% Hmisc::capitalize(include)), dat))
  if(is.numeric(facet)){facet <- switch(facet, c('measure', 'rule'), c('rule', 'measure'))}
  facet <- as.formula(switch(length(facet), paste('~', facet), paste(facet, collapse = ' ~ ')))
  if('box' %in% atts){
    if(is.logical(rule)){legend <- rule}
    g <- ggplot(dat, aes_string(x = 'type', y = 'value', fill = 'type')) +
      geom_boxplot() + facet_grid(facet) + theme_bw()
  } else {
    errs <- colnames(dat)[startsWith(colnames(dat), 's')]
    rule <- match.arg(tolower(rule), c(rules, errs))
    err <- ifelse(rule %in% errs, rule, match.arg(tolower(err), errs))
    yMin <- paste0('mean - ', err)
    yMax <- paste0('mean + ', err)
    g <- ggplot(dat, aes_string(x = 'type', y = 'mean', colour = 'type')) +
      geom_point() + geom_errorbar(aes_string(ymin = yMin, ymax = yMax), width = .2, lwd = .7) +
      facet_grid(facet) + theme_bw()
  }
  if(!is.null(title)){g <- g + ggtitle(title)}
  if(identical(ylab, FALSE)){g <- g + ylab('')}
  if(identical(xlab, FALSE)){g <- g + xlab('')}
  if(identical(ticks, FALSE) & !identical(legend, FALSE)){g <- g + theme(axis.text.x = element_blank())}
  if(is.character(legend)){g <- g + theme(legend.position = legend)}
  if(identical(legend, FALSE)){g <- g + theme(legend.position = 'none')}
  if(plot){return(g)} else {return(invisible(g))}
}

##### multiplot (JUST PLAYING AROUND)
multiplot <- function(x, err = 'se', legend = NULL, ylab = FALSE, nrow = 1,
                      ncol = 2, xlab = FALSE, ticks = FALSE, title = TRUE,
                      getLegend = FALSE){
  rules <- c('OR', 'AND')
  P <- lapply(seq_along(rules), function(z){
    plotit(dat = x, rule = z, legend = 'bottom', plot = FALSE, 
           title = switch(2 - isTRUE(title), rules[z], NULL), 
           err = err, ylab = ylab, xlab = xlab, ticks = ticks)
  })
  if(is.null(legend) | getLegend){legend <- g_legend(P[[1]])}
  if(getLegend){return(legend)}
  P <- lapply(P, function(z) z + theme(legend.position = 'none'))
  gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = P, nrow = nrow, ncol = ncol), 
                          legend, nrow = 2, heights = c(10, 1))
}

##### checkit
checkit <- function(x1, x2, one = TRUE){
  x3 <- lapply(c('out3', 'out4'), function(z){
    n <- lapply(x1[[z]], naby, 1, T)
    sapply(seq_along(n), function(i){
      if(length(n[[i]]) != 0){
        identical(x1[[z]][[i]][-n[[i]], ], x2[[z]][[i]][-n[[i]], ])
      } else {
        identical(x1[[z]][[i]], x2[[z]][[i]])
      }
    })
  })
  if(one){all(unlist(x3))} else {x3}
}

##### agg: the aggregate function is so annoying
agg <- function(data, y = 'value', by = NULL, metrics = c('mean', 'se', 'ci'), 
                ci = .95, otherFUN = NULL, ex = NULL, naZero = FALSE){
  dn <- which(apply(data, 2, function(z) length(unique(z))) == 1)
  if(length(dn) > 0){data <- data[, -dn, drop = FALSE]}
  if(is.null(by)){by <- setdiff(colnames(data), ex)}
  if(any(y %in% by)){by <- intersect(setdiff(by, y), colnames(data))}
  mets <- c('mean', 'sd', 'se', 'se2', 'ci')
  if('all' %in% metrics){metrics <- setdiff(mets, 'ci')}
  metrics <- unique(gsub('ci', 'se2', match.arg(tolower(metrics), mets, several.ok = TRUE)))
  mean1 <- function(x){mean(x, na.rm = TRUE)}
  sd1 <- function(x){sd(x, na.rm = TRUE)}
  se1 <- function(x){sd(x, na.rm = TRUE)/sqrt(length(x))}
  se21 <- function(x, alpha = (1 - ci)){se1(x) * qnorm(1 - alpha/2)}
  fun <- function(x, picks = metrics){
    c(mean = mean1(x), sd = sd1(x), se = se1(x), se2 = se21(x))[picks]
  }
  if(!is.null(otherFUN) & is.function(otherFUN)){fun <- otherFUN}
  X <- switch(2 - isTRUE(length(by) == 1), setNames(list(data[, by]), by), as.list(data[, by]))
  out <- aggregate(data[, y], X, fun)
  nn <- colnames(out)[-ncol(out)]
  out <- cbind.data.frame(out[, 1:length(nn)], out[, ncol(out)])
  colnames(out) <- c(nn, gsub('se2', paste0('ci', ci * 100), metrics))
  if(naZero & any(is.na(out))){out[is.na(out)] <- 0}
  out
}

##### multiNets: trying to solve a problem...
multiNets <- function(x, nn = c('x', 'fit1', 'fit2', 'fit3', 'fit4', 'fit5'), ...){
  xx <- lapply(x, function(z) evalFits(z, getNets = TRUE, ...))
  K <- lapply(xx, function(X) setNames(lapply(nn, function(z){
    z0 <- lapply(X, lapply, lapply, '[[', z)
    z1 <- appd(z0$pairwise, recursive = TRUE)
    z2 <- appd(z0$interactions, recursive = TRUE)
    out <- list(pairwise = z1, interactions = z2)
    return(out)
  }), nn))
  kcount <- function(x){sum(x != 0)/2}
  out <- lapply(K, lapply, function(z) do.call(
    cbind.data.frame, lapply(z, sapply, kcount)))
  out <- lapply(out, function(z) cbind.data.frame(
    do.call(rbind, z), fit = factor(rep(nn, sapply(z, nrow)), levels = nn)))
  out <- cbind.data.frame(do.call(rbind, out), N = rep(
    as.numeric(gsub('[^0-9]', '', names(x))), sapply(out, nrow)))
  rownames(out) <- 1:nrow(out)
  out
}

##### plotting: load ggplot and gridExtra
plotting <- function(){
  library(ggplot2)
  library(gridExtra)
  fixfactors <- function(x, use = 'all', except = NULL, addm = FALSE, 
                         addtype = TRUE){
    realcands <- c('P', 'M', 'N')
    cands <- names(which(sapply(x, is.factor)))
    if(!is.null(cands)){
      if(identical(tolower(use), 'all')){use <- cands}
      use <- intersect(use, cands)
      if(!is.null(except)){use <- setdiff(use, except)}
      use <- intersect(use, realcands)
      stopifnot(length(use) > 0)
      if('P' %in% use){levels(x$P) <- paste0('P = ', levels(x$P))}
      if('N' %in% use){levels(x$N) <- paste0('N = ', levels(x$N))}
      if('M' %in% use){
        levels(x$M) <- paste0(as.numeric(as.character(levels(x$M))) * 100, '%')
        if(addm){levels(x$M) <- paste0('M = ', levels(x$M))}
      }
      if(addtype & !'type' %in% colnames(x) & 'fit' %in% colnames(x)){
        ty <- ifelse(any(grepl('^NCT$|^FGL$', levels(x$fit))), 'Interactions', 'Pairwise')
        x$type <- factor(rep(ty, nrow(x)))
      }
      return(x)
    } else {
      message('No factors to fix!')
    }
  }
  plots <- function(x, measures = NULL, means = FALSE, lines = TRUE, r = 'AND',
                    p = 5, pfactor = TRUE, nfactor = FALSE, m = NULL, 
                    fix = FALSE, ...){
    elgnis <- !isTRUE(length(unique(as.character(x$fit))) == 1)
    if(!elgnis){elgnis <- list(fit = unique(as.character(x$fit)), levs = levels(x$fit))}
    if('iter' %in% colnames(x)){x <- x[, setdiff(colnames(x), 'iter')]}
    if('M' %in% colnames(x)){x$M <- factor(x$M, levels = sort(unique(x$M)))}
    if('rule' %in% colnames(x)){
      x <- subset(x, rule == r)
      x$rule <- factor(x$rule)
    }
    if(means){x <- agg(x, 'value')}
    if(!is.null(p)){x <- subset(x, P == p)}
    if(!is.null(m)){x <- subset(x, M == m)}
    if(pfactor & 'P' %in% colnames(x)){x$P <- factor(x$P, levels = sort(unique(x$P)))}
    if(!is.null(measures)){
      measures <- match.arg(tolower(measures), tolower(unique(x$measure)), several.ok = TRUE)
      if(any(startsWith(measures, 'm'))){
        measures[startsWith(measures, 'm')] <- toupper(measures[startsWith(measures, 'm')])
      }
      if(all(measures %in% c('ll', 'df', 'aic', 'bic'))){
        measures[measures != 'df'] <- toupper(measures[measures != 'df'])
      } else {
        measures <- Hmisc::capitalize(measures)
      }
      x <- subset(x, measure %in% measures)
      x$measure <- factor(x$measure, levels = measures)
    }
    if(nfactor){x$N <- factor(x$N, levels = sort(unique(x$N)))}
    rownames(x) <- 1:nrow(x)
    if(!isTRUE(elgnis) & !'fit' %in% colnames(x)){
      x <- cbind.data.frame(fit = factor(rep(elgnis$fit, nrow(x)), levels = elgnis$levs), x)
    }
    if(isTRUE(fix)){
      args <- tryCatch({list(...)}, error = function(e){list()})
      args$x <- x
      args <- args[intersect(names(args), formalArgs('fixfactors'))]
      x <- do.call(fixfactors, args)
    }
    x
  }
  studies <- function(out = 1, study = 1, beyond = FALSE, newRes = FALSE,
                      folder = '~/Desktop/ultimate', fit = FALSE, m5 = FALSE){
    stopifnot(out %in% 1:6 & study %in% 1:2)
    if(grepl('/$', folder)){folder <- gsub('/$', '', folder)}
    newRes <- switch(study, FALSE, !isTRUE(beyond) & newRes)
    path <- paste0(folder, '/study', study, ifelse(beyond, '/beyond', ifelse(newRes, '/newRes', '')), '/out', out, '.RDS')
    if(fit){path <- paste0(folder, '/fit', study, '.RDS')}
    file <- readRDS(path)
    if(study == 1 & m5){
      file2 <- readRDS(paste0('~/Desktop/NEWDATA/FINAL/new/out', out, '.RDS'))
      file <- rbind.data.frame(file, file2)
      file <- file[order(file$N), ]
    }
    file
  }
  grrange <- function(g, leg = NULL, method = 1, nrow = NULL, ncol = NULL, 
                      nrow2 = NULL, ncol2 = NULL, heights = c(10, 1)){
    if(!is(g, 'list')){return(plot(g))}
    method <- ifelse(is.null(leg), 2, ifelse(is.numeric(method), method, 1))
    method <- pmax(1, pmin(method, 2))
    if(method == 1){
      if(is.null(nrow2)){nrow2 <- 2}
      stopifnot(length(heights) == nrow2)
      grid.arrange(arrangeGrob(grobs = g, nrow = nrow, ncol = ncol), 
                   leg, nrow = nrow2, ncol = ncol2, heights = heights)
    } else {
      if(!is.null(leg)){g[[length(g) + 1]] <- leg}
      grid.arrange(arrangeGrob(grobs = g, nrow = nrow, ncol = ncol))
    }
  }
  funs <- list(plots = plots, fixfactors = fixfactors, 
               studies = studies, grrange = grrange)
  list2env(funs, .GlobalEnv)
}


##### finalFiles
finalFiles <- function(study = 3, type = 3, newRes = FALSE){
  s1_1 <- '~/Desktop/study1/results/results'
  s1_2 <- '~/Desktop/study1/resampling/results'
  s2_1 <- '~/Desktop/study2/final/results/results'
  s2_2 <- '~/Desktop/study2/fullOutput/finalRes/results'
  out <- list(study1 = list(fits = s1_1, resample = s1_2), 
              study2 = list(fits = s2_1, resample = s2_2))
  if(study < 3){
    out <- out[[study]]
    if(type < 3){out <- out[[type]]}
  }
  if(!identical(newRes, FALSE)){
    out <- '~/Desktop/X/FINAL/'
    if(isTRUE(newRes)){newRes <- 'files'}
    newRes <- match.arg(tolower(newRes), c('files', 'path'))
    if(newRes == 'files'){out <- dir(out)}
  }
  return(out)
}


##### count0
count0 <- function(x, equal = TRUE, find = FALSE){
  if(all(grepl('iter', names(x)))){
    `%0%` <- Negate(ifelse(equal, '!=', '=='))
    xx <- unlist(sapply(x, sapply, length))
    if(find){
      out <- which(apply(xx, 2, function(z) any(z == 0)))
    } else {
      out <- sum(c(xx) %0% 0)
    }
    return(out)
  } else if(all(grepl('m[1-3]$', names(x)))){
    out <- sapply(x, function(z){
      do.call(count0, list(x = z, equal = equal))
    })
    return(out)
  }
}

##### isFALSE
isFALSE <- function(x){identical(x, FALSE)}


##### relmets
relmets <- function(fit, true = NULL){
  # Relative parameter estimate: (avg(ests) - true)/true
  # Relative standard error: (avg(ses) - sd(ests))/sd(ests)
  print('Not created yet')
}