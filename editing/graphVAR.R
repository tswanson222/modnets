trevGVAR <- function(data, nLambda = 50, verbose = TRUE, gamma = 0.5, scale = TRUE, 
                     lambda_beta = NULL, lambda_kappa = NULL, maxit.in = 100, maxit.out = 100, 
                     deleteMissings = TRUE, penalize.diagonal = TRUE, lambda_min_kappa = 0.05, 
                     lambda_min_beta = lambda_min_kappa, mimic = c("current", "0.1.2", "0.1.4", "0.1.5", "0.2"), 
                     vars, beepvar, dayvar, idvar, lags = 1, centerWithin = TRUE, trevor = FALSE,
                     likelihood = c("unpenalized", "penalized"), manual = FALSE){
  if(manual == TRUE){
    nLambda = 50
    verbose = TRUE
    gamma = .5
    scale = T
    lambda_min_beta <- lambda_min_kappa <- .001
    lambda_beta <- lambda_kappa <- NULL
    mimic = "current"
    if(class(data) == "list"){vars = data$vars}
    else{vars = colnames(data)[!colnames(data) %in% "ID"]}
    centerWithin = T
    likelihood = "unpenalized"
    manual = T
    penalize.diagonal = T
    maxit.in <- maxit.out <- 100
    deleteMissings = T
    lags = 1
  } else {
    mimic <- match.arg(mimic)
  }
  if(mimic == "0.1.2"){
    if(lambda_min_beta != lambda_min_kappa){warning("mimic = 0.1.2 only uses lambda_min_kappa, not lambda_min_beta")}
    if(lambda_min_kappa != 0.01){warning("Set lambda_min_kappa = 0.01 to mimic 0.1.2 default behavior")}
  }
  if(is.list(data) && !is.data.frame(data)){
    if(!("data_c" %in% names(data) & "data_l" %in% names(data))){stop("'data_c' and 'data_l' must be contained in 'data'")}
    data_c <- data$data_c
    data_l <- data$data_l
  } else {
    if(mimic == "0.1.5"){
      if(is.data.frame(data)){data <- as.matrix(data)}
      stopifnot(is.matrix(data))
      data <- scale(data, TRUE, scale)
      data_c <- data[-1, , drop = FALSE]
      data_l <- cbind(1, data[-nrow(data), , drop = FALSE])
    } else {
      if(manual == FALSE & !trevor){
        data <- graphicalVAR:::tsData(as.data.frame(data), vars = vars, 
                                      beepvar = beepvar, dayvar = dayvar, idvar = idvar, 
                                      scale = scale, centerWithin = centerWithin, lags = lags)
      } else if(trevor){
        vars0 <- ifelse(missing(vars), list(NULL), list(vars))[[1]]
        data <- trevTS(as.data.frame(data), vars = vars0, scale = scale, centerWithin = centerWithin, lags = lags)
      } else {
        data <- graphicalVAR:::tsData(as.data.frame(data), vars = vars, 
                                      scale = scale, centerWithin = centerWithin, lags = lags)
      }
      data_c <- data$data_c
      data_l <- data$data_l
    }
  }
  data_c <- as.matrix(data_c)
  data_l <- as.matrix(data_l)
  Nvar <- ncol(data_c)
  Ntime <- nrow(data_c)
  if(any(is.na(data_c)) || any(is.na(data_l))){
    if(deleteMissings){
      warnings("Data with missings deleted")
      missing <- rowSums(is.na(data_c)) > 0 | rowSums(is.na(data_l)) > 0
      data_c <- data_c[!missing, ]
      data_l <- data_l[!missing, ]
    } else {
      stop("Missing data not supported")
    }
  }
  if(is.null(lambda_beta) | is.null(lambda_kappa)){
    if(mimic == "0.1.2"){
      if(!trevor){
        lams <- graphicalVAR:::SparseTSCGM_lambdas(data_l, data_c, nLambda, lambda.min.ratio = lambda_min_kappa)
      } else {
        lams <- trev_sparseTSCGM_lambdas(data_l, data_c, nLambda, lambda.min.ratio = lambda_min_kappa)
      }
    } else {
      if(!trevor){
        lams <- graphicalVAR:::generate_lambdas(data_l, data_c, nLambda, nLambda, 
                                                lambda_min_kappa = lambda_min_kappa, 
                                                lambda_min_beta = lambda_min_beta, 
                                                penalize.diagonal = penalize.diagonal, 
                                                version0.1.4 = mimic == "0.1.4")
      } else {
        lams <- trev_generateLambdas(data_l, data_c, nLambda, nLambda, 
                                     lambda_min_kappa = lambda_min_kappa, 
                                     lambda_min_beta = lambda_min_beta, 
                                     penalize.diagonal = penalize.diagonal, 
                                     version0.1.4 = mimic == "0.1.4")
      }
    }
    if(is.null(lambda_beta)){lambda_beta <- lams$lambda_beta}
    if(is.null(lambda_kappa)){lambda_kappa <- lams$lambda_kappa}
  }
  Nlambda_beta <- length(lambda_beta)
  Nlambda_kappa <- length(lambda_kappa)
  lambdas <- expand.grid(kappa = lambda_kappa, beta = lambda_beta)
  Estimates <- vector("list", nrow(lambdas))
  if(verbose){pb <- txtProgressBar(0, nrow(lambdas), style = 3)}
  for(i in seq_len(nrow(lambdas))){
    if(lambdas$beta[i] == 0 & lambdas$kappa[i] == 0){
      X <- data_l
      Y <- data_c
      nY <- ncol(Y)
      nX <- ncol(X)
      n <- nrow(X)
      beta <- t(Y) %*% X %*% solve(t(X) %*% X)
      S <- 1/(nrow(Y) - 1) * (t(Y) %*% Y - t(Y) %*% X %*% t(beta) - beta %*% t(X) %*% Y + beta %*% t(X) %*% X %*% t(beta))
      S <- (S + t(S))/2
      if(any(eigen(S)$value < 0)){stop("Residual covariances not postive definite")}
      kappa <- solve(S)
      kappa <- (kappa + t(kappa))/2
      lik1 = determinant(kappa)$modulus[1]
      lik2 <- sum(diag(kappa %*% S))
      pdO = sum(sum(kappa[upper.tri(kappa, diag = FALSE)] != 0))
      pdB = sum(sum(beta != 0))
      LLk <- (n/2) * (lik1 - lik2)
      LLk0 <- (n/2) * (-lik2)
      EBIC <- -2 * LLk + (log(n)) * (pdO + pdB) + (pdO + pdB) * 4 * gamma * log(2 * nY)
      Estimates[[i]] <- list(beta = beta, kappa = kappa, EBIC = EBIC)
    } else {
      tryres <- try(Rothmana(data_l, data_c, lambdas$beta[i], 
                             lambdas$kappa[i], gamma = gamma, maxit.in = maxit.in, 
                             maxit.out = maxit.out, penalize.diagonal = penalize.diagonal, 
                             mimic = mimic, likelihood = likelihood))
      if(is(tryres, "try-error")){
        Estimates[[i]] <- list(beta = matrix(NA, Nvar, Nvar + 1), 
                               kappa = matrix(NA, Nvar, Nvar), 
                               EBIC = Inf, error = tryres)
      } else {
        Estimates[[i]] <- tryres
      }
    }
    if(verbose){setTxtProgressBar(pb, i)}
  }
  if(verbose){close(pb)}
  lambdas$ebic <- sapply(Estimates, "[[", "EBIC")
  if(all(lambdas$ebic == Inf)){stop("No model estimated without error")}
  min <- which.min(lambdas$ebic)
  Results <- Estimates[[min]]
  if(length(lambda_beta) > 1){
    if(lambdas$beta[[min]] == min(lambda_beta)){message("Minimal tuning parameter for beta selected.")}
  }
  if(length(lambda_kappa) > 1){
    if(lambdas$kappa[[min]] == min(lambda_kappa)){message("Minimal tuning parameter for kappa selected.")}
  }
  colnames(Results$beta) <- colnames(data_l)
  rownames(Results$beta) <- colnames(data_c)
  Results$PCC <- trevPCC(Results$kappa)
  if(1 %in% lags){
    Results$PDC <- trevPDC(Results$beta[, c("1", paste0(data$vars, "_lag1"))], Results$kappa)
    if(length(lags) > 1){warning("Partial directed correlations only computed for lag 1 network.")}
  }
  Results$path <- lambdas
  Results$labels <- colnames(data_c)
  if(is.null(Results$labels)){Results$labels <- paste0("V", seq_len(ncol(data_c)))}
  rownames(Results$beta) <- colnames(Results$kappa) <- rownames(Results$kappa) <- colnames(Results$PCC) <- rownames(Results$PCC) <- colnames(Results$PDC) <- rownames(Results$PDC) <- Results$labels
  Results$gamma <- gamma
  Results$allResults <- Estimates
  Results$N <- nrow(data_c)
  Results$data <- data
  class(Results) <- "graphicalVAR"
  return(Results)
}

################################################################################
Rothmana <- function(X, Y, lambda_beta, lambda_kappa, convergence = 0.0001, 
                     gamma = 0.5, maxit.in = 100, maxit.out = 100, penalize.diagonal, 
                     interceptColumn = 1, mimic = "current", 
                     likelihood = c("unpenalized", "penalized")){
  likelihood <- match.arg(likelihood)
  nY <- ncol(Y)
  nX <- ncol(X)
  if(missing(penalize.diagonal)){
    penalize.diagonal <- ifelse(mimic == "0.1.2", nY != nX, (nY != nX - 1) & (nY != nX))
  }
  lambda_mat <- matrix(lambda_beta, nX, nY)
  if(!penalize.diagonal){
    if(nY != nX & nY != nX - 1){stop("Beta is not P x P or P x P+1, cannot detect diagonal.")}
    add <- ifelse(nY == nX, 0, ifelse(nY == nX - 1, 1, NA))
    for(i in 1:min(c(nY, nX))){lambda_mat[i + add, i] <- 0}
  }
  if(!is.null(interceptColumn) && !is.na(interceptColumn)){lambda_mat[interceptColumn, ] <- 0}
  n <- nrow(X)
  beta_ridge <- graphicalVAR:::beta_ridge_C(X, Y, lambda_beta)
  beta <- matrix(0, nX, nY)
  it <- 0
  repeat{
    it <- it + 1
    kappa <- graphicalVAR:::Kappa(beta, X, Y, lambda_kappa)
    beta_old <- beta
    beta <- graphicalVAR:::Beta_C(kappa, beta, X, Y, lambda_beta, lambda_mat, convergence, maxit.in)
    if(sum(abs(beta - beta_old)) < (convergence * sum(abs(beta_ridge)))){break}
    if(it > maxit.out){warning("Model did NOT converge in outer loop"); break}
  }
  ZeroIndex <- which(kappa == 0, arr.ind = TRUE)
  WS <- (t(Y) %*% Y - t(Y) %*% X %*% beta - t(beta) %*% t(X) %*% Y + t(beta) %*% t(X) %*% X %*% beta)/(nrow(X))
  if(any(eigen(WS, only.values = TRUE)$values < -sqrt(.Machine$double.eps))){stop("Residual covariance matrix is not non-negative definite")}
  if(likelihood == "unpenalized"){
    if(nrow(ZeroIndex) == 0){
      out4 <- suppressWarnings(glasso::glasso(WS, rho = 0, trace = FALSE))
    } else {
      out4 <- suppressWarnings(glasso::glasso(WS, rho = 0, zero = ZeroIndex, trace = FALSE))
    }
    lik1 <- determinant(out4$wi)$modulus[1]
    lik2 <- sum(diag(out4$wi %*% WS))
  } else {
    lik1 <- determinant(kappa)$modulus[1]
    lik2 <- sum(diag(kappa %*% WS))
  }
  pdO = sum(sum(kappa[upper.tri(kappa, diag = FALSE)] != 0))
  pdB <- ifelse(mimic == "0.1.2", sum(sum(beta != 0)), sum(sum(beta[lambda_mat != 0] != 0)))
  LLk <- (n/2) * (lik1 - lik2)
  LLk0 <- (n/2) * (-lik2)
  EBIC <- -2 * LLk + (log(n)) * (pdO + pdB) + (pdO + pdB) * 4 * gamma * log(2 * nY)
  return(list(beta = t(beta), kappa = kappa, EBIC = EBIC))
}

################################################################################
trevMLGVAR <- function(data, vars = NULL, idvar = NULL, subjectNetworks = TRUE, useTrevor = FALSE, 
                       beepvar = NULL, dayvar = NULL, scale = TRUE, centerWithin = TRUE, 
                       gamma = 0.5, verbose = TRUE, lambda_min_kappa_fixed = 0.001,
                       lambda_min_beta_fixed = 0.001, lambda_min_kappa = 0.05, 
                       lambda_min_beta = lambda_min_kappa, lambda_min_glasso = 0.01, ...){
  if(is.null(idvar)){stop("'idvar' must be assigned")}
  if(is.null(beepvar) & is.null(dayvar)){
    dataPrepped <- trevTS(data, vars = vars, beepvar = beepvar, dayvar = dayvar, idvar = idvar, scale = scale, centerWithin = centerWithin)
  } else {
    dataPrepped <- graphicalVAR:::tsData(data, vars = vars, beepvar = beepvar, dayvar = dayvar, idvar = idvar, scale = scale, centerWithin = centerWithin)
  }
  if(verbose){message("Estimating fixed networks")}
  if(useTrevor){
    ResFixed <- trevGVAR(dataPrepped, lambda_min_kappa = lambda_min_kappa_fixed, lambda_min_beta = lambda_min_beta_fixed, gamma = gamma, ...)
  } else{
    ResFixed <- graphicalVAR::graphicalVAR(dataPrepped, lambda_min_kappa = lambda_min_kappa_fixed, lambda_min_beta = lambda_min_beta_fixed, gamma = gamma, ...)
  }
  if(verbose){message("Estimating between-subjects network")}
  meansData <- dataPrepped$data_means
  meansData <- meansData[, names(meansData) != idvar]
  meansData <- meansData[rowMeans(is.na(meansData)) != 1, ]
  ResBetween <- qgraph::EBICglasso(cov(meansData), nrow(meansData), gamma, returnAllResults = TRUE, lambda.min.ratio = lambda_min_glasso)
  IDs <- unique(dataPrepped$data[[idvar]])
  idResults <- list()
  if(!identical(subjectNetworks, FALSE)){
    if(isTRUE(subjectNetworks)){subjectNetworks <- IDs}
    if(verbose){
      message("Estimating subject-specific networks")
      pb <- txtProgressBar(max = length(subjectNetworks), style = 3)
    }
    for(i in seq_along(subjectNetworks)){
      if(useTrevor == FALSE){
        capture.output({
          idResults[[i]] <- try(suppressWarnings(graphicalVAR::graphicalVAR(dataPrepped$data[dataPrepped$data[[idvar]] == subjectNetworks[i], ], 
                                                                            vars = dataPrepped$vars, beepvar = dataPrepped$beepvar, 
                                                                            dayvar = dataPrepped$dayvar, idvar = dataPrepped$idvar, 
                                                                            scale = scale, lambda_min_kappa = lambda_min_kappa, 
                                                                            lambda_min_beta = lambda_min_beta, gamma = gamma, 
                                                                            centerWithin = centerWithin, ..., verbose = FALSE)))
        })
      } else {
        capture.output({
          idResults[[i]] <- try(suppressWarnings(trevGVAR(dataPrepped$data[dataPrepped$data[[idvar]] == subjectNetworks[i], ], 
                                                          vars = dataPrepped$vars, beepvar = dataPrepped$beepvar, 
                                                          dayvar = dataPrepped$dayvar, idvar = dataPrepped$idvar, 
                                                          scale = scale, lambda_min_kappa = lambda_min_kappa, 
                                                          lambda_min_beta = lambda_min_beta, gamma = gamma, 
                                                          centerWithin = centerWithin, ..., verbose = FALSE)))
        })
      }
      if(verbose){setTxtProgressBar(pb, i)}
      if(is(idResults[[i]], "try-error")){idResults[[i]] <- list()}
    }
    if(verbose){close(pb)}
  } else {
    idResults <- lapply(seq_along(IDs), function(x) list())
  }
  Results <- list(fixedPCC = ResFixed$PCC, fixedPDC = ResFixed$PDC, 
                  fixedResults = ResFixed, betweenNet = ResBetween$optnet, 
                  betweenResults = ResBetween, ids = IDs, 
                  subjectPCC = lapply(idResults, "[[", "PCC"), 
                  subjectPDC = lapply(idResults, "[[", "PDC"), 
                  subjecResults = idResults)
  class(Results) <- "mlGraphicalVAR"
  return(Results)
}

################################################################################
trevTS <- function(data, vars = NULL, beepvar = NULL, dayvar = NULL, idvar = NULL, 
                   lags = 1, scale = TRUE, centerWithin = TRUE, deleteMissings = TRUE){
  . <- NULL
  data <- as.data.frame(data)
  invisible(suppressMessages(require(dplyr)))
  if(is.null(idvar)){
    idvar <- "ID"
    data[[idvar]] <- 1
  }
  if(is.null(dayvar)){
    dayvar <- "DAY"
    data[[dayvar]] <- 1
  }
  if(is.null(beepvar)){
    beepvar <- "BEEP"
    data <- data %>% dplyr::group_by_(dayvar, idvar) %>% dplyr::mutate_(BEEP = ~seq_len(n()))
  }
  if(is.null(vars)){vars <- names(data[!names(data) %in% c(idvar, dayvar, beepvar)])}
  data <- data[, c(vars, idvar, dayvar, beepvar)]
  data[, vars] <- scale(data[, vars], TRUE, scale)
  MeansData <- data %>% dplyr::group_by_(idvar) %>% dplyr::summarise_at(funs(mean(., na.rm = TRUE)), .vars = vars)
  if(centerWithin){
    if(length(unique(data[[idvar]])) > 1){
      data <- data %>% dplyr::group_by_(idvar) %>% dplyr::mutate_at(funs(scale(., center = TRUE, scale = FALSE)), .vars = vars)
    }
  }
  augData <- data
  beepsPerDay <- eval(substitute(dplyr::summarize_(data %>% group_by_(idvar, dayvar), first = ~min(beepvar, na.rm = TRUE), 
                                                   last = ~max(beepvar, na.rm = TRUE)), list(beepvar = as.name(beepvar))))
  allBeeps <- expand.grid(unique(data[[idvar]]), unique(data[[dayvar]]), seq(min(data[[beepvar]], na.rm = TRUE), max(data[[beepvar]], na.rm = TRUE)))
  names(allBeeps) <- c(idvar, dayvar, beepvar)
  allBeeps <- eval(substitute({
    allBeeps %>% dplyr::left_join(beepsPerDay, by = c(idvar, dayvar)) %>% dplyr::group_by_(idvar, dayvar) %>% 
      dplyr::filter_(~BEEP >= first, ~BEEP <= last) %>% 
      dplyr::arrange_(idvar, dayvar, beepvar)
  }, list(BEEP = as.name(beepvar))))
  augData <- augData %>% dplyr::right_join(allBeeps, by = c(idvar, dayvar, beepvar))
  data_c <- augData %>% ungroup %>% dplyr::select_(.dots = vars)
  shift <- function(x, n){length <- length(x); c(rep(NA, n), x)[1:length]}
  data_l <- do.call(cbind, lapply(lags, function(l){
    data_lagged <- augData %>% dplyr::group_by_(idvar, dayvar) %>% 
      dplyr::mutate_at(funs(shift), .vars = vars) %>% ungroup %>% 
      dplyr::select_(.dots = vars)
    names(data_lagged) <- paste0(vars, "_lag", l)
    data_lagged
  }))
  if(deleteMissings){
    isNA <- rowSums(is.na(data_c)) > 0 | rowSums(is.na(data_l)) > 0
    data_c <- data_c[!isNA, ]
    data_l <- data_l[!isNA, ]
  }
  Results <- list(data = augData, data_c = data_c[, vars], 
                  data_l = cbind(1, data_l), data_means = MeansData, vars = vars, 
                  idvar = idvar, dayvar = dayvar, beepvar = beepvar, lags = lags)
  class(Results) <- "tsData"
  return(Results)
}

################################################################################
trevTS2 <- function(data, vars = NULL, idvar = NULL, scale = TRUE, centerWithin = TRUE, lags = 1){
  data <- as.data.frame(data)
  if(is.null(idvar)){idvar <- "ID"; if(!idvar %in% colnames(data)){data[[idvar]] <- 1}}
  if(is.null(vars)){vars <- colnames(data)[!colnames(data) %in% idvar]}
  data[,vars] <- apply(data[,vars], 2, scale, center = TRUE, scale = scale)
  ids <- sort(unique(data[,idvar]))
  dataByID <- lapply(ids, function(z) data[data[,idvar] == z, vars])
  dataMeans <- cbind.data.frame(ID = ids, do.call(rbind, lapply(dataByID, function(z) colMeans(z))))
  if(centerWithin){
    dataByID <- lapply(dataByID, function(z) apply(z, 2, scale, center = TRUE, scale = FALSE))
    data[,vars] <- do.call(rbind, dataByID)
  }
  shift <- function(x, n){c(rep(NA, n), x)[1:length(x)]}
  data_l <- do.call(rbind, lapply(dataByID, function(z) apply(z, 2, shift, lags)))
  data_c <- data[!rowSums(is.na(data_l)) > 0, vars]
  data_l <- cbind.data.frame(1, data_l[!rowSums(is.na(data_l)) > 0, ])
  colnames(data_l)[colnames(data_l) %in% vars] <- paste0(vars, "_lag", lags)
  Results <- list(data = data, data_c = data_c, data_l = data_l, data_means = dataMeans,
                  vars = vars, idvar = idvar, lags = lags)
  Results
}

################################################################################
trevSim_MLGVAR <- function(nTime, nVar, nPerson, propPositive = 0.5, 
                           kappaRange = c(0.25, 0.5), betaRange = c(0.25, 0.5), 
                           betweenRange = c(0.25, 0.5), rewireWithin = 0, 
                           betweenVar = 1, withinVar = 0.25, temporalOffset = 2){
  repeat{
    trueKappa <- as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1, nVar, 1, 0)))
    trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] * sample(c(-1, 1), sum(upper.tri(trueKappa)), TRUE, prob = c(propPositive, 1 - propPositive))
    trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]
    trueBeta <- diag(1, nVar)
    for(i in 1:nVar){trueBeta[(i + (temporalOffset - 1))%%nVar + 1, i] <- sample(c(-1, 1), 1, propPositive)}
    trueBetween <- as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1, nVar, 1, 1)))
    trueBetween[upper.tri(trueBetween)] <- trueBetween[upper.tri(trueBetween)] * sample(c(-1, 1), sum(upper.tri(trueBetween)), TRUE, prob = c(propPositive, 1 - propPositive))
    trueBetween[upper.tri(trueBetween)] <- runif(sum(upper.tri(trueBetween)), betweenRange[1], betweenRange[2]) * trueBetween[upper.tri(trueBetween)]
    trueBetween[lower.tri(trueBetween)] <- t(trueBetween)[lower.tri(trueBetween)]
    diag(trueBetween) <- 1
    evK <- round(eigen(trueKappa)$values, 10)
    evB <- round(eigen(trueBeta)$values, 10)
    evBet <- round(eigen(trueBetween)$values, 10)
    if(all(evBet > 0)){break}
  }
  Sigma <- cov2cor(solve(trueBetween))
  D <- diag(sqrt(betweenVar), nVar)
  Means <- mvtnorm::rmvnorm(nPerson, sigma = D %*% Sigma %*% D)
  SubjectData <- lapply(1:nPerson, function(i){
    try <- 1
    maxtry <- 10
    repeat{
      kappa <- trueKappa
      kappa[upper.tri(kappa)] <- runif(sum(upper.tri(kappa)), kappaRange[1], kappaRange[2]) * kappa[upper.tri(kappa)]
      kappa[lower.tri(kappa)] <- t(kappa)[lower.tri(kappa)]
      diag(kappa) <- 1
      beta <- trueBeta * runif(nVar^2, betaRange[1], betaRange[2])
      kappa <- trevRewire(kappa, rewireWithin)
      beta <- trevRewire(beta, rewireWithin)
      evK <- eigen(kappa)$values
      evB <- eigen(beta)$values
      while(any(Re(evB)^2 + Im(evB)^2 > 1)){
        warning("Shrinking parameters")
        beta <- 0.95 * beta
        evB <- eigen(beta)$values
      }
      if(all(evK > 0) & all(Re(evB)^2 + Im(evB)^2 < 1)){break}
      try <- try + 1
      if(try > maxtry){stop("Maximum number of tries reached.")}
    }
    D <- diag(sqrt(withinVar), nVar)
    Delta <- diag(1/sqrt(diag(solve(kappa))))
    kappa <- solve(D) %*% solve(Delta) %*% kappa %*% solve(Delta) %*% solve(D)
    Data <- as.data.frame(trevSim_GVAR(nTime, beta, kappa, mean = Means[i, ]))
    Data$ID <- i
    return(list(kappa = kappa, beta = beta, PCC = trevPCC(kappa), 
                PDC = trevPDC(beta, kappa), data = Data))
  })
  fixedKappa <- Reduce("+", lapply(SubjectData, "[[", "kappa"))/nPerson
  fixedBeta <- Reduce("+", lapply(SubjectData, "[[", "beta"))/nPerson
  allData <- do.call(rbind, lapply(SubjectData, "[[", "data"))
  Results <- list(data = allData, fixedKappa = fixedKappa, 
                  fixedPCC = trevPCC(fixedKappa), fixedBeta = fixedBeta, 
                  fixedPDC = trevPDC(fixedBeta, fixedKappa), between = trueBetween, 
                  means = Means, personData = SubjectData, idvar = "ID", 
                  vars = names(allData)[names(allData) != "ID"])
  #class(Results) <- "simMLgvar"
  return(Results)
}

################################################################################
trevSim_GVAR <- function(nTime, beta = NULL, kappa = NULL, mean = rep(0, ncol(kappa)), 
                         init = mean, warmup = 100, lbound = rep(-Inf, ncol(kappa)), 
                         ubound = rep(Inf, ncol(kappa))){
  stopifnot(!is.null(beta))
  stopifnot(!is.null(kappa))
  Nvar <- ncol(kappa)
  init <- rep(init, length = Nvar)
  totTime <- nTime + warmup
  Data <- t(matrix(init, Nvar, totTime))
  Sigma <- solve(kappa)
  for(t in 2:totTime){
    Data[t, ] <- mean + t(beta %*% (Data[t - 1, ] - mean)) + mvtnorm::rmvnorm(1, rep(0, Nvar), Sigma)
    Data[t, ] <- ifelse(Data[t, ] < lbound, lbound, Data[t, ])
    Data[t, ] <- ifelse(Data[t, ] > ubound, ubound, Data[t, ])
  }
  return(Data[-seq_len(warmup), , drop = FALSE])
}


################################################################################
################################################################################
################################################################################
trevPCC <- function(x){
  x <- -cov2cor(x)
  diag(x) <- 0
  x <- as.matrix(Matrix::forceSymmetric(x))
  return(x)
}

################################################################################
trevPDC <- function(beta, kappa){
  if(ncol(beta) == nrow(beta) + 1){beta <- beta[, -1, drop = FALSE]}
  sigma <- solve(kappa)
  t(beta/sqrt(diag(sigma) %o% diag(kappa) + beta^2))
}

################################################################################
trevRewire <- function(x, p, directed){
  if(missing(directed)){directed <- !all(x == t(x))}
  if(directed){
    ind <- diag(1, ncol(x)) != 1
  } else {
    ind <- upper.tri(x)
  }
  curEdges <- which(x != 0 & ind, arr.ind = TRUE)
  toRewire <- which(runif(nrow(curEdges)) < p)
  for(i in seq_along(toRewire)){
    curZeros <- which(x == 0 & ind, arr.ind = TRUE)
    dest <- sample(seq_len(nrow(curZeros)), 1)
    x[curZeros[dest, 1], curZeros[dest, 2]] <- x[curEdges[toRewire[i], 1], curEdges[toRewire[i], 2]]
    x[curEdges[toRewire[i], 1], curEdges[toRewire[i], 2]] <- 0
    if(!directed){
      x[curZeros[dest, 2], curZeros[dest, 1]] <- x[curEdges[toRewire[i], 2], curEdges[toRewire[i], 1]]
      x[curEdges[toRewire[i], 2], curEdges[toRewire[i], 1]] <- 0
    }
  }
  return(x)
}

################################################################################
trev_sparseTSCGM_lambdas <- function(X, Y, nlambda = 100, lambda.min.ratio = 0.01){
  lambda.seq <- function(SS, SA, nlambda){
    if(length(nlambda) == 1){nlambda <- rep(nlambda, 2)}
    d = dim(SS)[2]
    lambda.max1 = max(max(SS - diag(d)), -min(SS - diag(d)))
    lambda.min1 = lambda.min.ratio * lambda.max1
    lambda1 = exp(seq(log(lambda.max1), log(lambda.min1), length = nlambda[1]))
    lambda.min.ratio2 = 0.15
    lambda.max2 = max(max(SA), -min(SA))
    lambda.min2 = lambda.min.ratio2 * lambda.max2
    lambda2 = exp(seq(log(lambda.max2), log(lambda.min2), length = nlambda[2]))
    return(list(lambda1 = lambda1, lambda2 = lambda2))
  }
  T <- dim(Y)[1]
  p <- dim(X)[2]
  n <- 1
  q <- dim(Y)[2]
  xtyi <- array(NA, c(p, q, n))
  xtxi <- array(NA, c(p, p, n))
  ytyi <- array(NA, c(q, q, n))
  XX <- X
  YY <- Y
  XX2 <- X^2
  YY2 <- Y^2
  xtyi <- crossprod(XX, YY)
  xtxi <- crossprod(XX)
  ytyi <- crossprod(YY)
  xty = apply(xtyi, c(1, 2), sum)
  xtx = apply(xtxi, c(1, 2), sum)
  yty = apply(ytyi, c(1, 2), sum)
  xtxt = apply(xtxi, c(1, 2), sum)/(n * T)
  xtx2 = (n * T) * colMeans(apply(XX2, c(1, 2), sum))
  yty2 = (n * T) * colMeans(apply(YY2, c(1, 2), sum))
  SX <- xtx/(n * T)
  mSX <- glasso::glasso(SX, 0.05, penalize.diagonal = FALSE)
  SX <- xtx/(n * T)
  mSX <- glasso::glasso(SX, 0.05, penalize.diagonal = FALSE)
  SXi <- mSX$wi
  SS = (yty)/(n * T)
  SS = cov2cor(SS)
  SAs = xty/(n * T)
  SA = t(SAs) %*% SXi
  lambda <- lambda.seq(SS = SS, SA = SA, nlambda = nlambda)
  lam1 <- round(lambda$lambda1, 3)
  lam2 <- round(lambda$lambda2, 3)
  lam2 <- round(lam2/max(lam2), 3)
  return(list(lambda_kappa = lam1, lambda_beta = lam2))
}

################################################################################
trev_generateLambdas <- function(X, Y, nLambda_kappa = 10, nLambda_beta = 10, 
                                 lambda_min_kappa = 0.05, lambda_min_beta = 0.05, 
                                 penalize.diagonal = TRUE, version0.1.4 = FALSE){
  N <- nrow(Y)
  P <- ncol(Y)
  corY <- cov2cor(t(Y) %*% Y/nrow(Y))
  if(version0.1.4){
    lam_K_max = max(abs(corY))
  } else {
    lam_K_max = max(abs(corY[upper.tri(corY)]))
  }
  lam_K_min = lambda_min_kappa * lam_K_max
  lam_K = exp(seq(log(lam_K_max), log(lam_K_min), length = nLambda_kappa))
  Yinv <- trevInvGlasso(t(Y) %*% Y)
  lam_B_max = max(abs(t(X) %*% Y %*% Yinv))
  lam_B_min = lambda_min_beta * lam_B_max
  lam_B = exp(seq(log(lam_B_max), log(lam_B_min), length = nLambda_beta))
  return(list(lambda_kappa = lam_K, lambda_beta = lam_B))
}

################################################################################
trevInvGlasso <- function(x){
  if(all(eigen(x)$values > sqrt(.Machine$double.eps))){
    Xinv <- solve(x)
  } else {
    Xglas <- glasso::glasso(x, 0.05, penalize.diagonal = FALSE)
    Xinv <- Xglas$wi
  }
  Xinv
}
