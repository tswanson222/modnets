##### simNet: cross-sectional simulations
simNet <- function(N = 100, p = 5, m = FALSE, m2 = .1, b1 = NULL, b2 = NULL, 
                   sparsity = .5, intercepts = NULL, nIter = 250, msym = FALSE, 
                   onlyDat = FALSE, pbar = TRUE, div = 1000, gibbs = TRUE,
                   ordinal = FALSE, nLevels = 5, mord = FALSE, time = TRUE,
                   mbinary = FALSE, minOrd = 3, m1 = 0, m1_range = NULL,
                   m2_range = c(.1, .3), modType = 'none', V = 2, 
                   lags = NULL, nCores = 1, cluster = 'SOCK', ...){
  t1 <- Sys.time()
  args <- tryCatch({list(...)}, error = function(e){list()})
  if(!is.null(lags) & !identical(lags, 0)){
    args1 <- append(list(nTime = N, nPerson = 1, nNode = p, lag = 1, GGMsparsity = sparsity), as.list(match.call())[-1])
    if('nPerson' %in% names(args)){args1 <- replace(args1, 'nPerson', args$nPerson)}
    args1 <- append(args1, args[setdiff(names(args), names(args1))])
    out <- do.call(mlGVARsim, args1[intersect(names(args1), formalArgs('mlGVARsim'))])
    return(out)
  }
  if(!is.null(m1_range)){if(!is.numeric(m1_range) | length(m1_range) != 2){m1_range <- NULL}}
  if(!is.numeric(m2_range) | length(m2_range) != 2){m2_range <- c(.1, .3)}
  if(is.numeric(modType)){modType <- switch(modType + 1, 'none', 'full', 'partial')}
  modType <- match.arg(tolower(modType), c('none', 'full', 'partial'))
  if(is.null(intercepts)){
    intercepts <- rep(0, p)
  } else if(length(intercepts) != p){
    intercepts <- rnorm(p)
  }
  mfun <- switch(
    2 - mbinary, 
    function(size = 1){sample(x = 0:1, size = size)}, 
    function(mean = 0){rnorm(n = 1, mean = mean)}
  )
  if(isTRUE(m)){m <- ifelse(mbinary, 1, mfun())}
  m <- ifelse(identical(m, FALSE), 0, m)
  if(!identical(m1, 0)){
    m1 <- runif(p, min(m2_range), max(m2_range))
    m1 <- m1 * sample(0:1, p, TRUE, prob = c(sparsity, 1 - sparsity))
  } else {
    m1 <- rep(0, p)
  }
  if(is.null(b1)){
    b1 <- simPcor(p, sparsity, finalRange = m1_range)
  }
  if(is.null(b2)){
    mat <- diag(0, p)
    pn <- (p * (p - 1))/2
    b2 <- sample(0:1, pn, TRUE, prob = c(1 - m2, m2))
    b2 <- b2 * sample(c(-1, 1), pn, TRUE, prob = c(.5, .5))
    mat[upper.tri(mat)] <- b2 * runif(pn, min(m2_range), max(m2_range))
    b2 <- as.matrix(Matrix::forceSymmetric(mat))
    if(msym){b2 <- simPcor(p, sparsity = 1 - m2, finalRange = m2_range)}
    if(all(m == 0)){b2 <- diag(0, p)}
  }
  if(modType != 'none' & !all(b2 == 0)){
    while(ifelse(modType == 'full', !all(b1[b2 != 0] == 0), any(b1[b2 != 0] == 0))){
      b1 <- simPcor(p, sparsity, finalRange = m1_range)
    }
  }
  diag(b1) <- diag(b2) <- 0
  if(nCores > 1 | isTRUE(nCores)){
    if(isTRUE(nCores)){nCores <- parallel::detectCores()}
    if(grepl('Windows', sessionInfo()$running)){cluster <- 'SOCK'}
    if(tolower(cluster) != 'mclapply'){
      cluster <- match.arg(toupper(cluster), c('SOCK', 'FORK'))
      cl <- parallel::makeCluster(nCores, type = cluster)
    } else {
      cl <- nCores
    }
    sampFun <- function(case, nIter, p, m, b1, b2, intercepts, div, V, m1){
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
          sampling[iter, v] <- rnorm(1, v_mu, 1)
          if(any(abs(sampling[!is.na(sampling)]) > div)){stop('Sampler diverged')}
        }
      }
      out <- list(data = sampling[nIter, ], M = m0[nIter])
      return(out)
    }
    if(cluster == 'SOCK'){
      parallel::clusterExport(cl, c('mfun', 'sampFun'), envir = environment())
    }
    out <- tryCatch({
      if(pbar){
        pbapply::pboptions(type = 'timer', char = '-')
        pbapply::pblapply(1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1, 
                          b2 = b2, intercept = intercepts, div = div, V = V, 
                          m1 = m1, cl = cl)
      } else if(tolower(cluster) != 'mclapply'){
        parallel::parLapply(
          cl, 1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1, b2 = b2, 
          intercepts = intercepts, div = div, V = V, m1 = m1)
      } else {
        parallel::mclapply(
          1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1, b2 = b2, 
          intercepts = intercepts, div = div, V = V, m1 = m1, mc.cores = nCores)
      }}, 
      error = function(e){TRUE}
    )
    if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
    if(isTRUE(out)){stop('Sampler diverged')}
    data <- do.call(rbind, lapply(out, '[[', 'data'))
    M <- unlist(lapply(out, '[[', 'M'))
    rm(cl)
  } else if(isTRUE(gibbs) | !all(m == 0)){
    data <- matrix(NA, nrow = N, ncol = p); M <- c()
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
          sampling[iter, v] <- rnorm(1, v_mu, 1)
          if(any(abs(sampling[!is.na(sampling)]) > div)){stop('Sampler diverged')}
        }
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
    data <- mvtnorm::rmvnorm(n = N, mean = intercepts, sigma = Sigma)
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
    if(mord){
      while(length(unique(mord)) < minOrd){
        mord <- as.numeric(cut(M, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
      }
      M <- mord
    }
    data <- cbind(data, M)
  }
  out <- append(list(data = data.frame(data)), out)
  if(onlyDat){out <- data.frame(data)}
  if(!all(m == 0)){attributes(out)[c('m', 'm1', 'm2')] <- list(m, m1, m2)}
  class(out) <- c(ifelse(onlyDat, 'data.frame', 'list'), 'mgmSim')
  attr(out, 'time') <- t2 <- Sys.time() - t1
  if(time){print(Sys.time() - t1)}
  return(out)
}


shit2 <- function(N, p, m, ...){
  out <- TRUE
  while(isTRUE(out)){
    out <- tryCatch({sim(N = N, p = p, m = m, ...)}, error = function(e){TRUE})
  }
  fit <- fitNetwork(out$data, p + 1)
  out <- list(out = list(b1 = out$b1, b2 = out$b2), fit = fit)
  return(out)
}

shit3 <- function(n = 10, N = 200, p = 5, V = 2, nCores = 32){
  t1 <- Sys.time()
  pb <- txtProgressBar(max = n, style = 3)
  out <- lapply(1:n, function(z){
    pp <- shit2(N = N, p = p, m = TRUE, nCores = nCores, V = V, time = FALSE)
    setTxtProgressBar(pb, z)
    if(z == n){close(pb)}
    return(pp)
  })
  print(Sys.time() - t1)
  return(out)
}

results <- function(x, threshold = FALSE, rule = 'or', avg = FALSE, ind = 'cor'){
  if('fit' %in% names(x)){x <- list(x)}
  y1 <- lapply(x, function(z) z$out$b1)
  y2 <- lapply(x, function(z) z$out$b2)
  x1 <- lapply(x, function(z) net(z$fit, threshold = threshold, rule = rule))
  x2 <- lapply(x, function(z) netInts(z$fit, threshold = threshold, rule = rule, avg = avg))
  out1 <- sapply(seq_along(x), function(z) matrixDist(x1[[z]], y1[[z]], ind = ind, directed = FALSE))
  out2 <- sapply(seq_along(x), function(z) matrixDist(x2[[z]], y2[[z]], ind = ind, directed = !avg))
  out <- cbind(net = out1, netInts = out2)
  return(out)
}
