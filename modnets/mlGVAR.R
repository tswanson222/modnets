### ======================================================================== ###
### ========================= Simulate MLVAR data ========================== ###
### ======================================================================== ###
##### VARsim: simulate VAR data for single subject
VARsim <- function(N, p, b1, b2 = NULL, c1 = NULL, c2 = NULL, covariate = NULL, 
                   b0 = 0, eps = 0.2, burnin = 100, getDat = TRUE, 
                   fullDat = FALSE, tvInit = FALSE){
  N <- N + burnin
  y <- matrix(0, N, p)
  if(!is(eps, "matrix")){
    sigma2 <- matrix(eps, p, p)
    diag(sigma2) <- 1
  } else {
    sigma2 <- eps
  }
  if(all(b0 == 0)){b0 <- sapply(rep(0, p), rep, N)}
  else if(length(b0) == 1){b0 <- sapply(rep(b0, p), rep, N)}
  else if(length(b0) == p){b0 <- sapply(b0, rep, N)}
  else if(length(b0) == N){b0 <- matrix(b0, N, p)}
  if(is.null(b2)){b2 <- rep(0, p)} else {if(length(b2) != p * (p - 1)){stop("Length of b2 needs to equal p * (p - 1)")}}
  if(is.null(c1)){c1 <- rep(0, p)} else {if(length(c1) != p){stop("Length of c1 needs to be equal to p")}}
  if(is.null(c2)){c2 <- rep(0, p)} else {if(length(c2) != p * p){stop("Length of c2 needs to be equal to p^2")}}
  if(is.null(covariate)){covariate <- rep(0, N)} else if(length(covariate) == 1){
    if(covariate == "c"){covariate <- sample(c(0, 1), N, replace = T)} 
    else if(covariate == "g" | covariate == TRUE){covariate <- rnorm(N)}
    else if(covariate == "rmean"){covariate <- rnorm(n = N, mean = runif(1))}
    else if(covariate == "ar"){
      covariate <- rnorm(1)
      if(is(eps, "matrix")){eps <- 0.2}
      for(i in 2:N){covariate[i] <- .1 * covariate[i - 1] + 1.5 * rnorm(n = 1, sd = sqrt(eps))}
    }
  }
  covariate2 <- matrix(NA, N, p)
  covariate3 <- covariate
  for(i in 1:p){covariate2[, i] <- covariate}
  if(class(b1) != "matrix"){
    if(length(b1) != p^2){stop("Length of b1 needs to equal p^2")}
    b1 <- sapply(b1, rep, N)
  }
  b2 <- sapply(b2, rep, N)
  c2 <- sapply(c2, rep, N)
  if(tvInit){
    pphi <- kronecker(matrix(b1[1, ], p, p, byrow = T), matrix(b1[1, ], p, p, byrow = T))
    mu1 <- solve(diag(p) - matrix(b1[1, ], p, p, byrow = T)) %*% matrix(b0[1, ], p, 1)
    s1 <- matrix(solve(diag(p * p) - pphi) %*% matrix(sigma2, ncol = 1), p, p)
  } else {
    mu1 <- rep(0, p)
    s1 <- sigma2
  }
  y[1, ] <- rockchalk::mvrnorm(1, mu = mu1, Sigma = s1)
  for(t in 2:N){
    a0 <- matrix(b0[t, ], p, 1)
    a1 <- matrix(b1[t, ], p, p, T) %*% y[t - 1, ]
    a2 <- matrix(b2[t, ], p, (p - 1), T) %*% (y[t - 1, 1] * y[t - 1, -1])
    a3 <- matrix(c1 * covariate2[t - 1, ], p, 1)
    a4 <- matrix(c2[t, ], p, p, T) %*% (y[t - 1, ] * covariate2[t - 1, ])
    e <- matrix(rockchalk::mvrnorm(1, rep(0, p), sigma2), p, 1)
    y[t, ] <- a0 + a1 + a2 + a3 + a4 + e
    if(any(is.na(y) | is.infinite(y))){stop(paste0("Failed at t = ", t))}
  }
  colnames(y) <- paste0("y", 1:p)
  if(all(covariate == 0)){covariate <- NULL} else {covariate <- covariate[-N]}
  betas <- list(b0 = matrix(b0[1, ], p, 1), b1 = matrix(b1[1, ], p, p, T),
                b2 = matrix(b2[1, ], p, p - 1, T), c1 = matrix(c1, p, 1),
                c2 = matrix(c2[1, ], p, p, T))
  if(all(b0 == 0)){betas$b0 <- NULL}
  if(all(b2 == 0)){betas$b2 <- NULL}
  if(all(c1 == 0)){betas$c1 <- NULL}
  if(all(c2 == 0)){betas$c2 <- NULL}
  out <- list(Y = y[-(1:(burnin + 1)), ], X = data.frame(y[-c(1:burnin, N), ]), betas = betas)
  colnames(out$X) <- paste0(colnames(out$X), "L")
  if(p == 1){colnames(out$X) <- "y1L"}
  if(!is.null(covariate)){out$X <- data.frame(out$X, covariate[-(1:burnin)])}
  if(getDat){
    if(p == 1){
      out$dat <- data.frame(y1 = out$Y, out$X)
    } else {
      out$dat <- list()
      for(i in 1:p){
        out$dat[[i]] <- data.frame(y = out$Y[, i], out$X)
        colnames(out$dat[[i]])[1] <- paste0("y", i)
      }
      names(out$dat) <- paste0("dat", 1:p)
    }
  }
  if(fullDat){
    if(!is.null(covariate)){bottom <- c(out$Y[(N - burnin) - 1, ], covariate3[N])}
    else{bottom <- c(out$Y[(N - burnin) - 1, ])}
    full <- data.frame(rbind(out$X, bottom))
    colnames(full) <- paste0("V", 1:ncol(full))
    out$full <- full
  }
  out
}

##### MLVARsim: simulate VAR data for multiple subjects
MLVARsim <- function(n, t, covariate = "rmean", intercept = "random", p = 2, 
                     b1 = NULL, c1 = NULL, c2 = NULL, getBetas = FALSE, eps = 0.2){
  covariate <- match.arg(covariate, c("true", "rmean", "c", "g", "false", "ar"))
  if(covariate == "true"){covariate <- as.logical(covariate)}
  if(covariate == "false"){covariate <- as.logical(covariate)}
  if(class(intercept) != "character"){b0 <- intercept}
  if(is.null(b1)){b1 <- c(.3, .2, .1, .4)} else if(all(b1 == 0)){b1 <- rep(0, p * p)}
  if(is.null(c1)){c1 <- c(.2, .15)} else if(all(c1 == 0)){c1 <- rep(0, p)}
  if(is.null(c2)){c2 <- c(0, .2, 0, .2)} else if(all(c2 == 0)){c2 <- rep(0, p * p)}
  if(getBetas){
    B <- rbind(cbind(matrix(b1, p, p, T), matrix(c1, p, 1), matrix(c2, p, p, T)), 0)
    colnames(B) <- c("V1.1", "V2.1", "V3.1", "V1.1:V3.1", "V2.1:V3.1")
    rownames(B) <- c("V1", "V2", "V3")
    return(B)
  }
  stopifnot(length(b1) == p * p)
  stopifnot(length(c1) == p)
  stopifnot(length(c2) == p * p)
  subjectData <- list()
  for(i in 1:n){
    if(class(intercept) == "character"){b0 <- rnorm(2)}
    dat <- VARsim(N = t, p = p, b1 = b1, covariate = covariate, b0 = b0, 
                  c1 = c1, c2 = c2, eps = eps, fullDat = TRUE)
    subjectData[[i]] <- dat$full
  }
  if(n == 1){return(dat)}
  ids <- rep(1:n, each = t)
  finalDat <- cbind.data.frame(do.call(rbind, subjectData), ID = ids)
  vars <- colnames(finalDat)[!colnames(finalDat) %in% "ID"]
  out <- list(data = finalDat, vars = vars, idvar = "ID")
  out
}


### ======================================================================== ###
### =========================== Fit GVAR models ============================ ###
### ======================================================================== ###
##### setupVAR: prepare data for analysis
setupVAR <- function(data, idvar = NULL, method = c("gvar", "lmer", "all"),
                     center = TRUE, scale = TRUE, vars = NULL, 
                     centerWithin = TRUE, scaleWithin = FALSE){
  method <- match.arg(method)
  if(class(data) == "list"){
    if(!"data" %in% names(data)){stop("Must supply data frame")}
    data <- as.data.frame(data[["data"]])
  }
  if(is.null(idvar)){
    if(!any(grepl("ID", colnames(data)))){stop("Must supply 'idvar'")}
    idvar <- colnames(data)[grep("ID", colnames(data))]
  }
  data <- data.frame(data[, -which(colnames(data) == idvar)], ID = data[, idvar])
  if(is.null(vars)){vars <- colnames(data)[!colnames(data) %in% "ID"]}
  ids <- unique(data[, "ID"])
  binary <- apply(data[, -which(colnames(data) == "ID")], 2, function(z) length(table(z)) <= 2)
  binary <- ifelse(any(binary), list(names(which(binary))), list(NULL))[[1]]
  vars0 <- setdiff(vars, binary)
  if(center){data[, vars0] <- apply(data[, vars0], 2, scale, center, scale)}
  dataByID <- lapply(ids, function(z) data[data[, "ID"] == z, seq_along(vars)])
  N <- rep(seq_along(ids), (sapply(dataByID, nrow) - ifelse(method == "all", 0, 1)))
  dataMeans <- do.call(rbind, lapply(dataByID, colMeans))[N, ]
  colnames(dataMeans) <- paste0(colnames(dataMeans), ".m")
  if(method != "lmer"){
    data0 <- lapply(dataByID, function(z){if(centerWithin){
      z[, vars0] <- apply(z[, vars0], 2, scale, centerWithin, scaleWithin)}
      return(z)
    })
    if(method == "all"){return(data.frame(do.call(rbind, data0), dataMeans, ID = ids[N]))}
    Y <- do.call(rbind, lapply(data0, function(z) z[-1, ]))
    X <- do.call(rbind, lapply(data0, function(z) z[-nrow(z), ]))
  } else {
    Y <- do.call(rbind, lapply(dataByID, function(z) z[-1, ]))
    X <- do.call(rbind, lapply(dataByID, function(z){
      z <- z[-nrow(z), ]
      if(centerWithin){z[, vars0] <- apply(
        z[, vars0], 2, scale, centerWithin, scaleWithin)}
      return(z)
    }))
  }
  dat <- data.frame(Y, X, dataMeans, ID = ids[N])
  dat
}

##### mlGVAR: fit GVAR models with multilevel data
mlGVAR <- function(data, m = NULL, selectFUN = NULL, subjectNets = FALSE, idvar = 'ID',
                   exogenous = TRUE, center = TRUE, scale = TRUE, fixedType = 'g', 
                   betweenType = 'g', centerWithin = TRUE, scaleWithin = FALSE,
                   rule = 'OR', threshold = 'none', verbose = TRUE, pcor = FALSE, 
                   fixedArgs = NULL, betweenArgs = NULL, bm = FALSE, beepno = NULL,
                   dayno = NULL, deleteMissing = TRUE, ...){
  t1 <- Sys.time()
  mnames <- mi <- m
  args0 <- tryCatch({list(...)}, error = function(e){list()})
  if(any(is.na(data))){
    if(deleteMissing){
      ww <- which(apply(data, 1, function(z) any(is.na(z))))
      data <- na.omit(data)
      warning(paste0(length(ww), ' rows deleted due to missingness'))
    } else {
      stop(paste0(length(ww), ' rows contain missing values'))
    }
  }
  if(!bm){
    bm <- list(moderators = NULL)
    betweenArgs <- switch(
      2 - is.null(betweenArgs), bm, 
      if(!'moderators' %in% names(betweenArgs)){append(betweenArgs, bm)} else {betweenArgs})
  }
  if(is.character(threshold)){
    threshold <- sapply(match.arg(tolower(
      threshold), c('none', 'pcc', 'between', 'fixed', 'all'), several.ok = TRUE), 
      function(z) switch(z, none = FALSE, all = TRUE, z), USE.NAMES = FALSE)
  }
  if(!idvar %in% colnames(data)){stop('Must supply idvar')}
  if(!is.null(m)){
    mnames <- switch(2 - is.character(m), m, colnames(data)[m])
  } else if(!is.null(selectFUN)){
    a0 <- names(args0)
    if(!'method' %in% a0){
      args0$method <- 'glmnet'
    } else if(args0$method == 'glinternet'){
      mnames <- m <- setdiff(colnames(data), idvar)
    }
    if(!'criterion' %in% a0){args0$criterion <- 'AIC'}
    if(any(c('resample', 'split') %in% selectFUN) & !'sampMethod' %in% a0){
      args0$sampMethod <- 'split'
    }
  }
  data <- data.frame(data[, -which(colnames(data) == idvar)], ID = data[, idvar])
  if(!is.null(beepno) | !is.null(dayno)){
    stopifnot(!is.null(beepno) & !is.null(dayno))
    stopifnot(length(beepno) == 1 & length(dayno) == 1)
    if(is.numeric(dayno)){dayno <- colnames(data)[dayno]}
    if(is.numeric(beepno)){beepno <- colnames(data)[beepno]}
    data0 <- data
    data <- data[, setdiff(colnames(data), c(beepno, dayno))]
  }
  vars <- setdiff(colnames(data), 'ID')
  dat <- setupVAR(data = data, idvar = 'ID', method = 'all', center = center,
                  scale = scale, vars = vars, centerWithin = centerWithin,
                  scaleWithin = scaleWithin)
  fixedDat <- dat[, vars]
  samp_ind <- as.numeric(cumsum(table(dat[, 'ID'])))
  if(!is.null(beepno) & !is.null(dayno)){
    dat <- cbind.data.frame(dat, data0[, c(beepno, dayno)])
    dat0 <- split(dat, dat$ID)
    consec <- vector('list', length(dat0))
    for(i in seq_along(dat0)){
      consec[[i]] <- which(!getConsec(data = dat0[[i]], beepno = beepno, dayno = dayno, makeAtt = FALSE))
      if(i > 1){
        consec[[i]] <- consec[[i]] + sum(sapply(1:(i - 1), function(z) nrow(dat0[[z]])))
      }
    }
    samp_ind <- union(unlist(consec), samp_ind)
  }
  attr(fixedDat, 'samp_ind') <- (1:nrow(dat))[-samp_ind]
  if(!is.null(m)){
    m <- which(vars %in% mnames)
    if(length(m) >= ncol(fixedDat) - 1){exogenous <- FALSE}
    if(!is.null(selectFUN)){args0$method <- 'glinternet'}
  }
  if(verbose){message('Estimating fixed networks')}
  fixedThresh <- ifelse(!is.character(threshold), threshold, ifelse(
    'fixed' %in% threshold, TRUE, ifelse('pcc' %in% threshold, 'PCC', FALSE)))
  fitNetArgs <- setdiff(formalArgs('fitNetwork'), '...')
  args1 <- list(data = fixedDat, moderators = m, type = fixedType, lags = 1,
                exogenous = exogenous, center = FALSE, scale = FALSE, pcor = pcor,
                rule = rule, threshold = fixedThresh, verbose = verbose)
  if(length(args0) > 0){args0 <- args0[setdiff(names(args0), names(args1))]}
  if(!is.null(fixedArgs)){
    fix1 <- intersect(names(fixedArgs), names(args1))
    if(length(fix1) > 0){args1 <- replace(args1, fix1, fixedArgs[fix1])}
    fixedArgs <- fixedArgs[setdiff(names(fixedArgs), names(args1))]
    args1 <- append(args1, fixedArgs[intersect(fitNetArgs, names(fixedArgs))])
  }
  if(is.null(selectFUN)){
    args1 <- append(args1, args0[intersect(fitNetArgs, names(args0))])
    fixedNets <- do.call(fitNetwork, args1)
  } else {
    if(isTRUE(selectFUN)){selectFUN <- 'varSelect'}
    if(length(selectFUN) == 2){betweenType <- isTRUE(selectFUN[2] == FALSE)}
    selectFUN <- match.arg(selectFUN[1], c('varSelect', 'resample', 'stability', 'bootstrap', 'split'))
    if(!selectFUN %in% c('varSelect', 'resample')){
      args0$sampMethod <- selectFUN
      selectFUN <- 'resample'
    }
    FUNargs <- setdiff(formalArgs(selectFUN), '...')
    FUN <- match.fun(selectFUN)
    args1 <- setNames(args1, gsub('moderators', 'm', names(args1)))
    args1.1 <- append(args1, args0)
    args1.1 <- args1.1[intersect(FUNargs, names(args1.1))]
    if('criterion' %in% names(args1.1)){
      if(length(args1.1$criterion) > 1){
        args1.1$criterion <- args1.1$criterion[1]}}
    if('gamma' %in% names(args1.1)){
      if(length(args1.1$gamma) > 1){args1.1$gamma <- args1.1$gamma[1]}
    }
    fixedType <- tryCatch({do.call(FUN, args1.1)}, error = function(e){TRUE})
    if(selectFUN == 'varSelect' | isTRUE(fixedType)){
      args1.2 <- append(replace(args1, c('type', 'verbose'), list(
        type = 'g', verbose = FALSE)), args0)
      if(!isTRUE(fixedType)){args1.2$type <- fixedType}
      names(args1.2)[names(args1.2) == 'm'] <- 'moderators'
      args1.2 <- args1.2[intersect(fitNetArgs, names(args1.2))]
      fixedNets <- tryCatch({do.call(fitNetwork, replace(args1.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, replace(args1.2, c('type', 'saveMods'), list(type = 'g', saveMods = FALSE)))})
    } else {
      if('fit0' %in% names(fixedType)){fixedType$fit0 <- NULL}
      if(is(fixedType, 'list')){attr(fixedType[[2]], 'threshold') <- fixedThresh}
      args1.2 <- append(list(obj = fixedType, data = args1.1$data, fit = TRUE), 
                        args0[intersect(c('select', 'thresh'), names(args0))])
      fixedNets <- tryCatch({do.call(modSelect, replace(args1.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, setNames(args1, gsub('^m$', 'moderators', names(args1))))
      })
      if(verbose){cat('\n')}
    }
  }
  ids <- unique(dat[, 'ID'])
  meansDat <- data.frame(do.call(rbind, lapply(ids, function(z){
    z <- dat[dat[, 'ID'] == z, paste0(vars, '.m')]
    return(z[1, ])
  })))
  colnames(meansDat) <- vars
  rownames(meansDat) <- 1:nrow(meansDat)
  if(verbose){message('Estimating between-subjects network')}
  if(nrow(meansDat) <= ncol(meansDat)){
    warning('Not enough subjects to fit unregularized between-subjects network')
  }
  betThresh <- ifelse(!is.character(threshold), threshold, 'between' %in% threshold)
  args2 <- list(data = meansDat, moderators = m, type = betweenType, center = FALSE,
                scale = FALSE, exogenous = exogenous, lags = NULL, rule = rule, 
                threshold = betThresh, pcor = pcor, verbose = verbose)
  if(!is.null(betweenArgs)){
    if('covariates' %in% names(betweenArgs) & !'moderators' %in% names(betweenArgs)){
      betweenArgs <- append(betweenArgs, list(moderators = NULL))
    }
    fix2 <- intersect(names(betweenArgs), names(args2))
    if(length(fix2) > 0){
      if('moderators' %in% fix2){
        if(is.null(betweenArgs$moderators) & !is.null(args2$moderators)){
          if(!'covariates' %in% names(betweenArgs) & exogenous){
            args2$data <- args2$data[, -args2$moderators]
          }
        }
      }
      args2 <- replace(args2, fix2, betweenArgs[fix2])
    }
    betweenArgs <- betweenArgs[setdiff(names(betweenArgs), names(args2))]
    args2 <- append(args2, betweenArgs[intersect(fitNetArgs, names(betweenArgs))])
  }
  if(is.null(selectFUN) | isTRUE(betweenType)){
    if(isTRUE(betweenType)){args2$type <- betweenType <- 'g'}
    args2 <- append(args2, args0[intersect(fitNetArgs, names(args0))])
    betNet <- do.call(fitNetwork, args2)
  } else {
    if(is.logical(betweenType)){
      betweenType <- 'g'
      selectFUN <- 'varSelect'
      FUNargs <- formalArgs('varSelect')
      FUN <- match.fun('varSelect')
    }
    args2 <- setNames(args2, gsub('moderators', 'm', names(args2)))
    args2.1 <- append(args2, args0)
    if(is.null(args2.1$m) & args2.1$method == 'glinternet'){args2.1$method <- 'glmnet'}
    args2.1 <- args2.1[intersect(FUNargs, names(args2.1))]
    if('criterion' %in% names(args2.1)){
      if(length(args2.1$criterion) > 1){
        args2.1$criterion <- args2.1$criterion[2]}}
    if('gamma' %in% names(args2.1)){
      if(length(args2.1$gamma) > 1){args2.1$gamma <- args2.1$gamma[2]}
    }
    betweenType <- tryCatch({do.call(FUN, args2.1)}, error = function(e){TRUE})
    if(selectFUN == 'varSelect' | isTRUE(betweenType)){
      args2.2 <- append(replace(args2, c('type', 'verbose'), list(
        type = 'g', verbose = FALSE)), args0)
      if(!isTRUE(betweenType)){args2.2$type <- betweenType}
      names(args2.2)[names(args2.2) == 'm'] <- 'moderators'
      args2.2 <- args2.2[intersect(fitNetArgs, names(args2.2))]
      betNet <- tryCatch({do.call(fitNetwork, replace(args2.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, replace(args2.2, c('type', 'saveMods'), list(type = 'g', saveMods = FALSE)))})
    } else {
      if('fit0' %in% names(betweenType)){betweenType$fit0 <- NULL}
      if(is(betweenType, 'list')){attr(betweenType[[2]], 'threshold') <- betThresh}
      args2.2 <- append(list(obj = betweenType, data = args2.1$data, fit = TRUE), 
                        args0[intersect(c('select', 'thresh'), names(args0))])
      betNet <- tryCatch({do.call(modSelect, replace(args2.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, setNames(args2, gsub('^m$', 'moderators', names(args2))))
      })
    }
  }
  if(all(subjectNets != FALSE)){
    subjectSelect <- ifelse(all(subjectNets == 'select'), TRUE, FALSE)
    if(isTRUE(subjectNets) | isTRUE(subjectSelect)){subjectNets <- ids}
    if(any(!subjectNets %in% ids)){subjectNets <- intersect(subjectNets, ids)}
    if(verbose){
      message('Estimating subject-specific networks')
      pb <- txtProgressBar(max = length(subjectNets), style = 3)
    }
    indNets <- lapply(seq_along(subjectNets), function(i){
      dati <- dat[dat[, 'ID'] == subjectNets[i], vars]
      if(!is.null(m)){
        mi <- sapply(seq_along(m), function(j) length(unique(dati[, m[j]])) != 1)
        mi <- ifelse(all(!mi), list(NULL), list(m[mi]))[[1]]
        if(!identical(m, mi)){dati <- dati[, -setdiff(m, mi)]}
      }
      p0 <- ncol(dati)
      fiti <- nrow(dati) > (2 * p0) + 1
      if(!is.null(mi)){
        p1 <- sum(sapply(seq_along(mi), function(p) p0 - p))
        fiti <- ifelse(exogenous & length(mi) < ncol(dati) - 1, 
                       nrow(dati) > p1 + (2 * p0) - (length(mi) - 1), 
                       nrow(dati) > p1 + (2 * p0) + 1)
      }
      if(fiti){
        indi <- fitNetwork(data = dati, moderators = mi, threshold = fixedThresh,
                           exogenous = exogenous, center = center, type = 'g', 
                           lags = 1, scale = ifelse(center, scale, FALSE), 
                           saveMods = FALSE, ...)
      }
      if(verbose){setTxtProgressBar(pb, i)}
      if(fiti){return(indi)} else {return(list())}
    })
    names(indNets) <- paste0('subject', subjectNets)
    if(any(sapply(indNets, length) == 0)){
      inderrs <- subjectNets[sapply(indNets, length) == 0]
      indNets <- indNets[which(sapply(indNets, length) > 0)]
      if(length(indNets) > 0){
        message(paste0('Too few time points to estimate networks for subject',
                       ifelse(length(inderrs) == 1, ': ', 's: '),
                       paste(inderrs, collapse = ', ')))
      } else {
        subjectNets <- FALSE
        message('Too few time points to estimate subject-specific networks')
      }
    }
  }
  outcall <- list(fixedType = fixedType, betweenType = betweenType, m = mnames, 
                  center = center, scale = scale, exogenous = exogenous, 
                  threshold = threshold, centerWithin = centerWithin, 
                  scaleWithin = scaleWithin, rule = rule)
  if(length(args0) > 0){outcall <- append(outcall, args0)}
  if(!is.null(fixedArgs) & length(fixedArgs) > 0){outcall <- append(outcall, list(fixedArgs = fixedArgs))}
  if(!is.null(betweenArgs) & length(betweenArgs) > 0){outcall <- append(outcall, list(betweenArgs = betweenArgs))}
  out <- list(call = outcall, fixedNets = fixedNets, betweenNet = betNet, ids = ids)
  if(!is.null(selectFUN)){
    out$ids <- NULL
    out$varMods <- list(fixedMods = fixedType, betweenMods = betweenType)
    out$call$fixedType <- fixedNets$call$type
    out$call$betweenType <- betNet$call$type
    out$call$selectFUN <- selectFUN
    out$ids <- ids
  }
  if(all(subjectNets != FALSE)){out$subjectNets <- indNets}
  out$netData <- list(samp_ind = samp_ind, fixedDat = fixedDat, meansDat = meansDat)
  out$dat <- dat
  attr(out, 'mlGVAR') <- TRUE
  class(out) <- c('list', 'mlGVAR')
  attr(out, 'time') <- t2 <- Sys.time() - t1
  if(verbose){cat('\nCOMPLETE:', round(t2, 2), attr(t2, 'units')); cat('\n')}
  return(out)
}

# NOTE: beepno & dayno only work for mlGVAR, not lmerVAR. Also, does not
# yet work with model selection procedures

### ======================================================================== ###
### ========================= Fit lmerVAR models =========================== ###
### ======================================================================== ###
##### lmerVAR: fit lmerVAR models
lmerVAR <- function(data, m = NULL, temporal = "default", contemp = "default",
                    idvar = "ID", intvars = NULL, center = TRUE, scale = TRUE, 
                    centerWithin = TRUE, scaleWithin = FALSE, exogenous = TRUE,
                    covariates = NULL, fix = NULL, warnings = FALSE, verbose = TRUE,
                    beepno = NULL, dayno = NULL, deleteMissing = TRUE){
  t1 <- Sys.time()
  suppressMessages(invisible(c(require(lme4), require(lmerTest))))
  if(!warnings){oldw <- getOption("warn"); options(warn = -1)}
  mnames <- m
  data <- data.frame(data)
  vars <- colnames(data)
  if(!idvar %in% vars){stop("Must supply 'idvar'")}
  if(!is.null(m)){mnames <- switch(2 - is.character(m), m, vars[m])}
  if(any(is.na(data))){
    if(deleteMissing){
      ww <- which(apply(data, 1, function(z) any(is.na(z))))
      data <- na.omit(data)
      warning(paste0(length(ww), ' rows deleted due to missingness'))
    } else {
      stop(paste0(length(ww), ' rows contain missing values'))
    }
  }
  data <- data.frame(data[, -which(vars == idvar)], ID = data[, idvar])
  if(!is.null(beepno) | !is.null(dayno)){
    stopifnot(!is.null(beepno) & !is.null(dayno))
    stopifnot(length(beepno) == 1 & length(dayno) == 1)
    if(is.numeric(dayno)){dayno <- colnames(data)[dayno]}
    if(is.numeric(beepno)){beepno <- colnames(data)[beepno]}
    data0 <- data
    data <- data[, setdiff(colnames(data), c(beepno, dayno))]
  }
  vars <- setdiff(colnames(data), idvar)
  dat <- setupVAR(data = data, idvar = "ID", method = "all", center = center,
                  scale = scale, centerWithin = FALSE, scaleWithin = FALSE)
  dat0 <- dat[, vars]
  samp_ind <- as.numeric(cumsum(table(dat[, "ID"])))
  if(!is.null(beepno) & !is.null(dayno)){
    dat <- cbind.data.frame(dat, data0[, c(beepno, dayno)])
    dat1 <- split(dat, dat$ID)
    consec <- vector('list', length(dat1))
    for(i in seq_along(dat1)){
      consec[[i]] <- which(!getConsec(data = dat1[[i]], beepno = beepno, dayno = dayno, makeAtt = FALSE))
      if(i > 1){
        consec[[i]] <- consec[[i]] + sum(sapply(1:(i - 1), function(z) nrow(dat1[[z]])))
      }
    }
    samp_ind <- union(unlist(consec), samp_ind)
  }
  nid <- attr(dat0, "samp_ind") <- setdiff(1:nrow(dat), samp_ind)
  ids0 <- unique(dat[, "ID"])
  ids <- dat[nid, "ID"]
  if(!is.null(m)){
    m <- match(mnames, vars)
    exogenous <- ifelse(length(m) >= length(vars) - 1, FALSE, exogenous)
  }
  dat0 <- lagMat(data = dat0, m = m, center = FALSE, scale = FALSE, 
                 exogenous = exogenous, lags = 1, checkType = TRUE, 
                 covariates = covariates)
  if(centerWithin){
    dat0$X[, vars] <- do.call(rbind, lapply(seq_along(ids0), function(i){
      apply(subset(dat0$X[, vars], ids == ids0[i]), 2, scale, TRUE, scaleWithin)
    }))
  }
  yvars <- colnames(dat0$Y)
  mvars <- paste0(vars, ".m")
  dat1 <- data.frame(dat0$Y, dat0$X, dat[nid, mvars], ID = ids, check.names = FALSE)
  if(is.null(intvars) & !is.null(m)){intvars <- colnames(dat1)[grep(":", colnames(dat1))]}
  temporal <- match.arg(temporal, c("default", "correlated", "orthogonal", "fixed", "intfixed"))
  if(temporal == "default"){temporal <- ifelse(length(yvars) > 6, "orthogonal", "correlated")}
  if(verbose){
    message(paste0("Estimating temporal and between-subject networks (", temporal, ")"))
    pb <- txtProgressBar(min = 0, max = length(yvars), style = 3)
  }
  tempForms <- lapply(yvars, function(y){
    x <- paste0(vars, collapse = " + ")
    means <- paste0(paste0(setdiff(
      gsub("[.]m$", "", mvars), 
      gsub("[.]y$", "", y)), ".m"), 
      collapse = " + ")
    ints <- paste0(intvars, collapse = " + ")
    ints <- switch(2 - isTRUE(ints == ""), NULL, ints)
    fixed <- paste0(c(x, means, ints), collapse = " + ")
    if(temporal == "correlated"){
      rands <- paste0("(", paste0(c(x, ints), collapse = " + "), " | ID)")
    } else {
      rands <- "(1 | ID)"
      randvars <- switch(2 - isTRUE(temporal == "intfixed"), vars, c(vars, intvars))
      if(!is.null(fix)){randvars <- setdiff(randvars, fix)}
      if(temporal != "fixed"){
        rands <- paste0("(", rands, " + ", paste0(
          paste0("(0 + ", randvars, " | ID)"), 
          collapse = " + "), ")")
      }
    }
    return(as.formula(paste0(y, " ~ ", fixed, " + ", rands)))
  })
  tempMods <- setNames(lapply(seq_along(yvars), function(z){
    tm <- suppressMessages(lmer(tempForms[[z]], data = dat1, REML = FALSE))
    if(verbose){setTxtProgressBar(pb, z)}
    return(tm)
  }), yvars)
  resDat <- structure(do.call(data.frame, lapply(tempMods, resid)), names = yvars)
  contemp <- match.arg(contemp, c("default", "correlated", "orthogonal"))
  if(contemp == "default"){contemp <- ifelse(length(yvars) > 6, "orthogonal", "correlated")}
  if(verbose){
    message(paste0("\nEstimating contemporaneous network (", contemp, ")"))
    pb <- txtProgressBar(min = 0, max = length(yvars), style = 3)
  }
  contempForms <- lapply(yvars, function(y){
    x <- setdiff(yvars, y)
    fixed <- paste0(c(0, x), collapse = " + ")
    if(contemp == "correlated"){
      rands <- paste0("(", fixed, " | ID)")
    } else {
      rands <- paste0("((1 | ID) + ", paste0(
        paste0("(0 + ", x, " | ID)"), collapse = " + "), ")")
    }
    return(as.formula(paste0(y, " ~ ", fixed, " + ", rands)))
  })
  resDat$ID <- ids
  contempMods <- setNames(lapply(seq_along(yvars), function(z){
    cm <- suppressMessages(lmer(contempForms[[z]], data = resDat, REML = FALSE))
    if(verbose){setTxtProgressBar(pb, z)}
    return(cm)
  }), yvars)
  fit <- do.call(rbind, lapply(tempMods, function(z) c(AIC(z), BIC(z))))
  fit <- structure(cbind.data.frame(yvars, fit), names = c("var", "aic", "bic"))
  inds <- list(yvars = yvars, vars = vars, mvars = mvars, intvars = intvars)
  outcall <- list(m = mnames, temporal = temporal, contemp = contemp, 
                  exogenous = exogenous, center = center, scale = scale, 
                  centerWithin = centerWithin, scaleWithin = scaleWithin)
  if(!is.null(covariates)){outcall <- append(outcall, list(covariates = covariates))}
  model <- list(tempMods = tempMods, contempMods = contempMods, fit = fit)
  out <- tryCatch({lmerNets(model = model, inds = inds, m = mnames)}, 
                  error = function(e){list(inds = inds)})
  for(i in seq_along(yvars)){
    attributes(model$tempMods[[i]])$formula <- tempForms[[i]]
    attributes(model$contempMods[[i]])$formula <- contempForms[[i]]
  }
  out <- append(list(call = outcall), append(out, list(mods = model, data = dat1)))
  attr(out, "temporal") <- paste0(temporal, ifelse(
    !is.null(intvars), " (interaction)", ifelse(
      !is.null(covariates), " (covariate)", "")))
  attr(out, "contemporaneous") <- contemp
  attr(out, "lmerVAR") <- TRUE
  class(out) <- c('list', 'lmerVAR')
  attr(out, "time") <- t2 <- Sys.time() - t1
  if(verbose){cat("\n"); print(Sys.time() - t1)}
  if(!warnings){options(warn = oldw)}
  return(out)
}

##### lmerNets: create networks from lmerVAR models
lmerNets <- function(model, inds, m = NULL, threshold = FALSE, 
                     rule = "OR", ggm = "pcor"){
  rule <- match.arg(tolower(rule), c("or", "and"))
  ggm <- match.arg(tolower(ggm), c("pcor", "cor", "cov", "prec"))
  forcePositive <- function(x){
    x <- (x + t(x))/2
    if(any(eigen(x)$values < 0)){
      x <- x - diag(nrow(x)) * min(eigen(x)$values) - 0.001
    }
    return(x)
  }
  y <- inds$yvars
  x <- inds$vars
  mvars <- inds$mvars
  intvars <- inds$intvars
  k <- length(y)
  y1 <- list(y, x)
  y2 <- rep(list(y), 2)
  ### BETA
  beta <- beta0 <- do.call(rbind, lapply(
    model$tempMods, function(z) fixef(z)[x]))
  betaSE <- betaSE0 <- do.call(rbind, lapply(
    model$tempMods, function(z) arm::se.fixef(z)[x]))
  beta_pvals <- betaPs0 <- (1 - pnorm(abs(beta/betaSE))) * 2
  dimnames(beta) <- dimnames(betaSE) <- dimnames(beta_pvals) <- y1
  dimnames(beta0) <- dimnames(betaSE0) <- dimnames(betaPs0) <- y1
  if(ncol(beta) != nrow(beta)){
    beta <- beta[y, match(y, paste0(x, ".y"))]
    betaSE <- betaSE[y, match(y, paste0(x, ".y"))]
    beta_pvals <- beta_pvals[y, match(y, paste0(x, ".y"))]
  }
  ### GAMMA THETA
  gammaTheta <- lapply(model$contempMods, fixef)
  gammaThetaSE <- lapply(model$contempMods, arm::se.fixef)
  gt1 <- gt2 <- matrix(NA, k, k)
  for(i in 1:k){
    gt1[i, match(names(gammaTheta[[i]]), y)] <- gammaTheta[[i]]
    gt2[i, match(names(gammaThetaSE[[i]]), y)] <- gammaThetaSE[[i]]
  }
  diag(gt1) <- diag(gt2) <- 0
  gammaTheta <- gt1
  gammaThetaSE <- gt2
  gammaTheta_pvals <- (1 - pnorm(abs(gammaTheta/gammaThetaSE))) * 2
  diag(gammaTheta_pvals) <- 1
  dimnames(gammaTheta) <- dimnames(gammaThetaSE) <- dimnames(gammaTheta_pvals) <- y2
  ### THETA AND PDC
  D <- diag(1/sapply(model$contempMods, sigma)^2)
  inv <- tryCatch({forcePositive(D %*% (diag(k) - gammaTheta))}, error = function(e){diag(k)})
  if(identical(inv, diag(k))){message("\nNon-convergent estimate of Theta")}
  Theta_prec <- inv
  Theta_cov <- corpcor::pseudoinverse(inv)
  Theta_cor <- cov2cor(Theta_cov)
  d <- 1/sqrt(diag(inv))
  Theta <- -t(d * inv) * d
  diag(Theta) <- 0
  dimnames(Theta) <- dimnames(Theta_cov) <- dimnames(Theta_cor) <- dimnames(Theta_prec) <- y2
  PDC <- beta/(sqrt(diag(Theta_cov) %o% diag(Theta_prec) + beta^2))
  colnames(beta) <- colnames(PDC) <- colnames(betaSE) <- colnames(beta_pvals) <- paste0(colnames(beta), ".lag1.")
  ### GAMMA OMEGA
  between <- TRUE
  if(length(y) != length(x)){
    if(!is.null(m)){
      if(length(setdiff(x, m)) < 3){between <- FALSE}
      if(between){mvars <- setdiff(mvars, paste0(m, ".m"))}
    } else {
      mvars <- mvars[gsub("[.]m$", ".y", mvars) %in% y]
    }
  }
  gammaOmega <- lapply(model$tempMods, function(z){
    z2 <- fixef(z)
    z2 <- z2[intersect(names(z2), mvars)]
    return(z2)
  })
  gammaOmegaSE <- lapply(model$tempMods, function(z){
    z2 <- arm::se.fixef(z)
    z2 <- z2[intersect(names(z2), mvars)]
    return(z2)
  })
  if(between){
    go1 <- go2 <- matrix(NA, k, k)
    for(i in 1:k){
      go1[i, match(names(gammaOmega[[i]]), mvars)] <- gammaOmega[[i]]
      go2[i, match(names(gammaOmegaSE[[i]]), mvars)] <- gammaOmegaSE[[i]]
    }
    diag(go1) <- diag(go2) <- 0
    gammaOmega <- go1
    gammaOmegaSE <- go2
    gammaOmega_pvals <- (1 - pnorm(abs(gammaOmega/gammaOmegaSE))) * 2
    diag(gammaOmega_pvals) <- 1
    dimnames(gammaOmega) <- dimnames(gammaOmegaSE) <- dimnames(gammaOmega_pvals) <- y2
    ### OMEGA
    mu_SD <- sapply(model$tempMods, function(z) attr(lme4::VarCorr(z)[[1]], "stddev")[1])
    D <- diag(1/mu_SD^2)
    inv <- tryCatch({forcePositive(D %*% (diag(k) - gammaOmega))}, error = function(e){diag(k)})
    if(identical(inv, diag(k))){message("\nNon-convergent estimate of Omega")}
    Omega_prec <- inv
    Omega_cov <- corpcor::pseudoinverse(inv)
    Omega_cor <- cov2cor(Omega_cov)
    d <- 1/sqrt(diag(inv))
    Omega <- -t(d * inv) * d
    diag(Omega) <- 0
  } else {
    gammaOmega <- gammaOmegaSE <- gammaOmega_pvals <- Omega <- Omega_cov <- Omega_cor <- Omega_prec <- diag(0, k)
    dimnames(gammaOmega) <- dimnames(gammaOmegaSE) <- dimnames(gammaOmega_pvals) <- y2
  }
  dimnames(Omega) <- dimnames(Omega_cov) <- dimnames(Omega_cor) <- dimnames(Omega_prec) <- y2
  ### COLLECT RESULTS
  out <- list(Beta = list(mu = beta, SE = betaSE, Pvals = beta_pvals), 
              Theta = list(cor = Theta_cor, cov = Theta_cov, pcor = Theta, prec = Theta_prec), 
              Omega = list(cor = Omega_cor, cov = Omega_cov, pcor = Omega, prec = Omega_prec), 
              gammaTheta = list(mu = gammaTheta, SE = gammaThetaSE, Pvals = gammaTheta_pvals),
              gammaOmega = list(mu = gammaOmega, SE = gammaOmegaSE, Pvals = gammaOmega_pvals))
  if(!is.null(m)){
    ints <- do.call(rbind, lapply(model$tempMods, function(z){
      fixef(z)[grepl(":", names(fixef(z)))]}))
    intsSE <- do.call(rbind, lapply(model$tempMods, function(z){
      arm::se.fixef(z)[grepl(":", names(arm::se.fixef(z)))]}))
    intsPvals <- (1 - pnorm(abs(ints/intsSE))) * 2
    out$ints <- list(coefs = ints, SE = intsSE, Pvals = intsPvals) 
  }
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    dimnames(colMat) <- dimnames(adjMat)
    colMat
  }
  ### THRESHOLDING
  if(threshold != FALSE){
    if(threshold == TRUE){threshold <- 0.05}
    thresh <- function(obj, alpha = 0.05, rule = "or", ggm = "pcor"){
      beta <- obj$Beta$mu
      theta <- obj$Theta[[ggm]]
      omega <- obj$Omega[[ggm]]
      bp <- obj$Beta$Pvals
      tp <- obj$gammaTheta$Pvals
      op <- obj$gammaOmega$Pvals
      diag(tp) <- diag(op) <- 1
      beta2 <- beta * ifelse(bp <= alpha, 1, 0)
      if("PDC" %in% names(obj)){pdc2 <- obj$PDC * ifelse(bp <= alpha, 1, 0)}
      if(rule == "or"){
        theta2 <- theta * ifelse(tp <= alpha | t(tp) <= alpha, 1, 0)
        omega2 <- omega * ifelse(op <= alpha | t(op) <= alpha, 1, 0)
      } else if(rule == "and"){
        theta2 <- theta * ifelse(tp <= alpha & t(tp) <= alpha, 1, 0)
        omega2 <- omega * ifelse(op <= alpha & t(op) <= alpha, 1, 0)
      }
      res <- list(beta = beta2, theta = theta2, omega = omega2)
      if("PDC" %in% names(obj)){res$PDC <- pdc2}
      return(res)
    }
    out3 <- thresh(obj = append(out, list(PDC = PDC)), 
                   alpha = threshold, rule = rule, ggm = ggm)
    beta <- out3$beta
    Theta <- out3$theta
    Omega <- out3$omega
    PDC <- out3$PDC
  }
  ### OUTPUT
  out2 <- list(temporal = list(adjMat = beta, edgeColors = getEdgeColors(beta)),
               contemporaneous = list(adjMat = Theta, edgeColors = getEdgeColors(Theta)),
               between = list(adjMat = Omega, edgeColors = getEdgeColors(Omega)))
  out2$temporal <- append(out2$temporal, list(PDC = list(adjMat = PDC, edgeColors = getEdgeColors(PDC))))
  out2$contemporaneous$kappa <- Theta_prec
  output <- append(out2, list(coefs = out))
  #output <- append(out2, list(coefs = out, mods = model[grepl("Mods|fit", names(model))]))
  attr(output, "lmerVAR") <- TRUE
  return(output)
}

##### compareVAR: compare lmerVAR models
compareVAR <- function(m1, m2, m3 = NULL, anova = NULL, type = "tempMods"){
  if(is.null(m3)){
    if(all(unlist(lapply(list(m1, m2), function(z) "temporal" %in% names(attributes(z)))))){
      cat(paste0("Model 1: ", attributes(m1)$temporal), "\n")
      cat(paste0("Model 2: ", attributes(m2)$temporal), "\n\n")
    }
    m1 <- m1$mods; m2 <- m2$mods
    aic <- apply(cbind(m1$fit$aic, m2$fit$aic), 1, which.min)
    bic <- apply(cbind(m1$fit$bic, m2$fit$bic), 1, which.min)
    if(all.equal(m1$fit$aic, m2$fit$aic) == TRUE){aic <- rep(0, nrow(m1$fit))}
    if(all.equal(m1$fit$bic, m2$fit$bic) == TRUE){bic <- rep(0, nrow(m1$fit))}
    return(data.frame(aic, bic))
  } else if(is.null(anova)){
    M <- list(m1, m2, m3)
    if(all(unlist(lapply(M, function(z) "temporal" %in% names(attributes(z)))))){
      names(M) <- c("Model 1", "Model 2", "Model 3")
      for(i in 1:3){cat(paste0(names(M)[i], ": ", lapply(M, attr, "temporal")[[i]]), "\n")}
      cat("\n")
    }
    m1 <- m1$mods; m2 <- m2$mods; m3 <- m3$mods
    vars <- m1$fit$var
    out <- list()
    for(i in 1:length(vars)){
      out[[i]] <- rbind(m1$fit[m1$fit$var == vars[i], ], m2$fit[m2$fit$var == vars[i], ],
                        m3$fit[m3$fit$var == vars[i], ])
    }
    return(cbind.data.frame(var = vars, do.call(rbind, lapply(out, function(z) apply(z[,-1], 2, which.min)))))
  } else if(is.numeric(anova) & anova <= length(m1$mods$tempMods)){
    m <- anova
    type <- match.arg(type, c("tempMods", "contempMods"))
    M1 <- m1$mods[[type]][[m]]
    M2 <- m2$mods[[type]][[m]]
    M3 <- m3$mods[[type]][[m]]
    return(anova(M1, M2, M3))
  }
}

compareVAR2 <- function(mods, p){
  vs <- paste0(names(mods), "$mods$tempMods[[", p, "]]")
  obj <- paste0("anova(", paste0(vs, collapse = ", "), ")")
  for(i in 1:length(vs)){assign(names(mods)[i], mods[[i]], pos = 1)}
  out <- eval(parse(text = obj))
  out
}


### ======================================================================== ###
### ============================= PLOTS ==================================== ###
### ======================================================================== ###
##### makePlot:
makePlot <- function(model, data, zxy = NULL, res = 50, radius = .03, 
                     points = TRUE, ppoints = TRUE, save = FALSE){
  if(is.null(zxy)){zxy <- as.list(colnames(data))}
  z <- zxy[[1]]
  x <- zxy[[2]]
  y <- zxy[[3]]
  if(is.null(names(zxy))){
    zlab <- z
    xlab <- x
    ylab <- y
  } else {
    zlab <- names(zxy)[1]
    xlab <- names(zxy)[2]
    ylab <- names(zxy)[3]
  }
  require(rgl)
  predictgrid <- function(model, xvar, yvar, zvar, res, type = NULL){
    xr <- range(model$model[[xvar]])
    yr <- range(model$model[[yvar]])
    newdata <- expand.grid(x = seq(xr[1], xr[2], length = res), y = seq(yr[1], yr[2], length = res))
    names(newdata) <- c(xvar, yvar)
    newdata[[zvar]] <- predict(model, newdata = newdata, type = type)
    newdata
  }
  df2mat <- function(p, xvar = NULL, yvar = NULL, zvar = NULL){
    if(is.null(xvar)){xvar <- names(p)[1]}
    if(is.null(yvar)){yvar <- names(p)[2]}
    if(is.null(zvar)){zvar <- names(p)[3]}
    x <- unique(p[[xvar]])
    y <- unique(p[[yvar]])
    z <- unique(p[[zvar]], nrow = length(y), ncol = length(x))
    m <- list(x, y, z)
    names(m) <- c(xvar, yvar, zvar)
    m
  }
  interleave <- function(v1, v2){as.vector(rbind(v1, v2))}
  data$pred <- predict(model)
  grid_df <- predictgrid(model, x, y, z, res = res)
  grid_list <- df2mat(grid_df)
  plot3d(data[,x], data[,y], data[,z], xlab = "", ylab = "", zlab = "", 
         axes = F, type = ifelse(points, "s", "n"), size = .4, lit = F)
  if(ppoints & points){
    spheres3d(data[,x], data[,y], data$pred, radius = radius, alpha = .4, type = "s", size = .5, lit = F)
    segments3d(interleave(data[,x], data[,x]), interleave(data[,y], data[,y]), 
               interleave(data[,z], data$pred), alpha = .4, col = "red")
  }
  surface3d(grid_list[[x]], grid_list[[y]], grid_list[[z]], 
            alpha = .4, front = "lines", back = "lines")
  rgl.bbox(color = "grey50", emission = "grey50", xlen = 0, ylen = 0, zlen = 0)
  rgl.material(color = "black")
  axes3d(edges = c("x--", "y+-", "z--"), ntick = 6, cex = .75)
  mtext3d(xlab, edge = "x--", line = 2)
  mtext3d(ylab, edge = "y+-", line = 3)
  mtext3d(zlab, edge = "z--", line = 3)
  if(save == TRUE){rgl.snapshot("~/Desktop/fuckery.png", fmt = "png")}
}

##### makePlot2:
makePlot2 <- function(data, zxy = NULL, center = TRUE, res = 50, effect = TRUE){
  if(is.null(zxy)){zxy <- as.list(colnames(data))}
  require(ggplot2)
  require(margins)
  if(is.null(names(zxy))){
    zlab <- zxy[[1]]
    xlab <- zxy[[2]]
    ylab <- zxy[[3]]
  } else {
    zlab <- names(zxy)[1]
    xlab <- names(zxy)[2]
    ylab <- names(zxy)[3]
  }
  z <- data[,zxy[[1]]]
  x <- data[,zxy[[2]]]
  y <- data[,zxy[[3]]]
  if(center){
    x <- x - mean(x)
    y <- y - mean(y)
  }
  dat <- data.frame(z, x, y)
  mod <- lm(z ~ x * y, data = dat)
  xvals <- seq(min(x), max(x), length = res)
  lowy <- mean(y) - sd(y)
  highy <- mean(y) + sd(y)
  lowdf <- data.frame(x = xvals, y = lowy)
  highdf <- data.frame(x = xvals, y = highy)
  pred_low <- predict(mod, data.frame(x = lowdf$x, y = lowdf$y), se.fit = T)
  yhat_low <- pred_low$fit
  se_low <- pred_low$se.fit * qnorm(.975)
  pred_low <- data.frame(xvals = xvals, zvals = yhat_low, upper = yhat_low + se_low, 
                         lower = yhat_low - se_low, y = -1)
  pred_high <- predict(mod, data.frame(x = highdf$x, y = highdf$y), se.fit = T)
  yhat_high <- pred_high$fit
  se_high <- pred_high$se.fit * qnorm(.975)
  pred_high <- data.frame(xvals = xvals, zvals = yhat_high, upper = yhat_high + se_high,
                          lower = yhat_high - se_high, y = 1)
  pred_mod <- rbind(pred_low, pred_high)
  pred_mod$y <- factor(pred_mod$y)
  colnames(pred_mod)[2] <- "zvals"
  levels(pred_mod$y) <- c(paste0("Low ", ylab), paste0("High ", ylab))
  pred_plot <- ggplot(data = pred_mod, aes(x = xvals, group = y)) +
    geom_line(aes(y = zvals, color = y)) + 
    geom_ribbon(alpha = .2, aes(ymin = lower, ymax = upper)) + 
    labs(title = paste0("Effect of ", xlab, " on ", zlab, " Moderated by ", ylab),
         subtitle = paste0("Simple Slopes at +/- 1 SD of Mean-Centered ", ylab),
         y = paste0("Predicted Values of ", zlab),
         x = paste0("Mean Centered Value of ", xlab),
         colour = "") + 
    theme_bw() +
    theme(legend.position = c(.1, .85))
  marg_mod <- margins::cplot(mod, dx = "x", x = "y", what = "effect", data = dat, draw = FALSE)
  marg_plot <- ggplot(data = marg_mod, aes(x = xvals, y = yvals)) +
    geom_line(color = "red") +
    geom_ribbon(alpha = .2, aes(ymin = lower, ymax = upper)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste0("Marginal Effect of ", xlab, " on ", zlab, " as a Function of ", ylab),
         subtitle = paste0("Marginal Effect of ", xlab, " Across Range of Mean-Centered ", ylab),
         y = paste0("Estimated Effect of ", xlab, " on ", zlab),
         x = paste0("Mean Centered Values of ", ylab),
         color = "") +
    theme_bw()
  if(effect){
    plot(marg_plot)
  } else {
    cowplot::plot_grid(pred_plot, marg_plot, ncol = 2, align = "h")
  }
}

##### moderationPlot:
moderationPlot <- function(data, vars = list("Y", "X", "M"), plot = TRUE, 
                           n = 20, mtype = "g", mean = TRUE){
  if(!is(vars, 'list')){vars <- as.list(vars)}
  if(length(vars) != 3){stop("Must provide 3 variable names")}
  if(is.null(names(vars))){names(vars) <- c("Y", "X", "M"); message("Assumed vars order: Y, X, M")}
  varCols <- sapply(vars, function(x) which(colnames(data) == x))
  colnames(data)[varCols] <- names(vars)
  X <- seq(min(data[,"X"]), max(data[,"X"]), length = n)
  form <- "Y ~ X * M"
  if(all(data$M %in% c(0, 1))){mtype <- "c"}
  if(mtype == "c"){
    dat <- list(datZero = data.frame(X = X, M = 0),
                datOne = data.frame(X = X, M = 1))
  } else {
    Mm <- mean(data[,"M"])
    Mlo <- Mm - sd(data[,"M"])
    Mhi <- Mm + sd(data[,"M"])
    dat <- list(datLo = data.frame(X = X, M = Mlo), 
                datM = data.frame(X = X, M = Mm), 
                datHi = data.frame(X = X, M = Mhi))
  }
  if(ncol(data) > 3){
    colnames(data)[-varCols] <- paste0("Z", 1:(ncol(data) - 3))
    zs <- grep("Z", colnames(data))
    for(i in 1:length(zs)){
      zs[i] <- mean(data[, zs[i]])
      dat <- lapply(dat, function(z){
        zz <- cbind(z, zs[i])
        colnames(zz)[i + 2] <- paste0("Z", i)
        return(zz)
      })
      form <- paste0(form, " + Z", i)
    }
  }
  mod <- lm(form, data = data)
  preds <- lapply(dat, function(z){
    zz <- predict(mod, z, se.fit = TRUE)
    zz$se.fit <- zz$se.fit * qnorm(.975)
    zz <- data.frame(X = X, Y = zz$fit, upper = zz$fit + zz$se.fit, lower = zz$fit - zz$se.fit)
    return(zz)
  })
  if(mtype == "c"){
    predMod <- cbind(do.call(rbind, preds), M = c(rep(0, n), rep(1, n)))
  } else {
    predMod <- cbind(do.call(rbind, preds), M = c(rep("Low", n), rep("Mean", n), rep("High", n)))
    if(mean == FALSE){predMod <- predMod[predMod$M != "Mean",]}
  }
  predMod$M <- factor(predMod$M)
  rownames(predMod) <- 1:nrow(predMod)
  if(plot == TRUE){
    require(ggplot2)
    predPlot <- ggplot(data = predMod, aes(x = X, group = M)) +
      geom_line(aes(y = Y, color = M)) +
      geom_ribbon(alpha = .2, aes(ymin = lower, ymax = upper)) +
      labs(title = "Effect of X on Y Moderated by M",
           subtitle = ifelse(mtype == "c", "Simple slopes at levels of M", "Simple slopes at +/- 1 SD of M"),
           y = "Predicted values of Y",
           x = "X",
           color = ifelse(mtype == "c", "M", "")) +
      theme_bw() +
      theme(legend.position = c(.1, .85))
    predPlot
  } else {
    predMod
  }
}

