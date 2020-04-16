### ======================================================================== ###
### ================= Simulate mlGVAR data (mlVAR Version) ================= ###
### ======================================================================== ###
##### mlGVARsim: main workhorse for simulating VAR and mlGVAR data
mlGVARsim <- function(nTime = 50, nPerson = 10, nNode = 3, m = NULL, m2 = .25, m1 = .7, 
                      m0 = 1, lag = 1, thetaVar = NULL, mu_SD = NULL, init_beta_SD = NULL, 
                      fixedMuSD = 1, shrink_fixed = 0.9, propPos = .5, m1SD = .1, m2SD = .1,
                      m1_range = NULL, m2_range = NULL, shrink_deviation = 0.9, getM = FALSE,
                      contemporaneous = "wishart", GGMsparsity = .5, mcenter = TRUE, skew = FALSE, 
                      skewErr = FALSE, ordinal = FALSE, nLevels = 5, ordWithin = TRUE, minOrd = 3, 
                      thresholds = NULL, mseed = NULL, onlyNets = FALSE, modType = 'none'){
  modType <- match.arg(tolower(modType), c('none', 'full', 'partial', 'full2', 'partial2', 'zero'))
  if(identical(m, FALSE)){m <- NULL}
  if(minOrd < 2){minOrd <- 2}
  if(is.numeric(ordinal)){nLevels <- ordinal; ordinal <- TRUE}
  if(is.null(thetaVar)){thetaVar <- rep(1, nNode)}
  if(is.null(mu_SD)){mu_SD <- c(1, 1)}
  if(is.null(init_beta_SD)){init_beta_SD <- c(.1, 1)}
  if(is.null(m1_range)){m1_range <- c(.1, .4)}
  if(is.null(m2_range)){m2_range <- c(.1, .3)}
  if(length(nTime) == 1){nTime <- rep(nTime, nPerson)}
  DF_theta <- nNode * 2
  nTemporal <- nNode^2 * lag
  contemporaneous <- match.arg(contemporaneous, c("wishart", "randomGGM", "fixed"))
  Omega_mu <- simCor(nNode, GGMsparsity)
  Omega_Beta <- simCor(nTemporal, GGMsparsity)
  Omega <- rbind(cbind(Omega_mu, matrix(0, nNode, nTemporal)), 
                 cbind(matrix(0, nTemporal, nNode), Omega_Beta))
  SD <- runif(nNode + nTemporal, c(rep(mu_SD[1], nNode), rep(init_beta_SD[1], nNode)), 
              c(rep(mu_SD[2], nNode), rep(init_beta_SD[2], nNode)))
  Omega <- diag(SD) %*% Omega %*% diag(SD)
  if(contemporaneous == "wishart"){
    Theta_fixed <- simCor(nNode, GGMsparsity)
    Theta_fixed <- diag(sqrt(thetaVar)) %*% Theta_fixed %*% diag(sqrt(thetaVar))
    Theta <- rWishart(nPerson, DF_theta, Theta_fixed/DF_theta)
  } else if(contemporaneous == "randomGGM"){
    Theta <- lapply(1:nPerson, function(x) simCor(nNode, GGMsparsity))
    Theta <- do.call(abind::abind, c(Theta, along = 3))
    Theta_fixed <- apply(Theta, 1:2, mean)
  } else {
    Theta_fixed <- simCor(nNode, GGMsparsity)
    Theta <- lapply(1:nPerson, function(x) Theta_fixed)
    Theta <- do.call(abind::abind, c(Theta, along = 3))
  }
  if(!is.null(m)){
    if(isTRUE(m)){m <- "random"}
    m <- match.arg(tolower(m), c(
      "fixed", "random", "mixed1", "mixed2", "ar", "binary", 
      "skewed", "random0", "ordinal"), several.ok = TRUE)
    if(all(c("fixed", "random") %in% m)){m <- m[-which(m == "fixed")]}
    if("ar" %in% m){
      if(m0 >= 1){m0 <- .3}
      mm <- lapply(seq_len(nPerson), function(z){
        as.numeric(arima.sim(n = nTime[z] + 100, model = list(ar = m0)))})
    } else if("binary" %in% m){
      mm <- lapply(seq_len(nPerson), function(z) rbinom(nTime[z] + 100, 1, .5))
    } else if("skewed" %in% m){
      mm <- lapply(seq_len(nPerson), function(z) sn::rsn(nTime[z] + 100, 0, 1, m0))
    } else {
      mm <- lapply(seq_len(nPerson), function(z) rnorm(nTime[z] + 100, 0, m0))
    }
    if(!is.null(mseed)){set.seed(mseed)}
    if(m2 >= 0 & m2 < 1){
      m20 <- diag(0, nNode)
      if(modType != 'zero'){
        wcond <- function(m20, modType){
          all(m20 == 0) | ifelse(grepl('full', modType), any(diag(m20) != 0), FALSE)
        }
        while(wcond(m20, modType)){
          m20 <- matrix(sample(c(0, 1), nNode^2, TRUE, prob = c(1 - m2, m2)), nNode, nNode)
        }
      }
    } else if(m2 >= 1){
      m20 <- numeric(nNode^2)
      m20[sample(1:(nNode^2), ifelse(m2 > nNode^2, nNode^2, round(m2)))] <- 1
      if(modType == 'zero'){
        m20 <- diag(0, nNode)
      } else if(grepl('full', modType)){
        m20 <- matrix(m20, nNode, nNode)
        while(any(diag(m20) != 0)){
          if(m2 >= (nNode^2) - nNode){
            m20 <- matrix(1, nNode, nNode)
            diag(m20) <- 0
          } else {
            m20 <- numeric(nNode^2)
            m20[sample(1:(nNode^2), round(m2))] <- 1
            m20 <- matrix(m20, nNode, nNode)
          }
        }
      }
    }
    if(isTRUE(m1) | is.null(m1)){m1 <- .7}
    if(m1 >= 0 & m1 < 1){
      m10 <- sample(c(0, 1), nNode, TRUE, prob = c(1 - m1, m1))
    } else if(m1 >= 1){
      m10 <- numeric(nNode)
      m10[sample(1:nNode, ifelse(m1 > nNode, nNode, round(m1)))] <- 1
    }
    if(grepl('2', modType) & !all(m20 == 0)){
      mm12 <- apply(m20, 1, function(z) any(z != 0))
      while(ifelse(grepl('full', modType), !all(m10[mm12] == 0), any(m10[mm12] == 0))){
        if(m1 >= 0 & m1 < 1){
          m10 <- sample(0:1, nNode, TRUE, prob = c(1 - m1, m1))
        } else {
          m10 <- numeric(nNode)
          m10[sample(1:nNode, ifelse(m1 > nNode, nNode, round(m1)))] <- 1
        }
      }
    }
    m1 <- m10 * sample(c(-1, 1), nNode, TRUE, prob = c(1 - propPos, propPos))
    m2 <- m20 * sample(c(-1, 1), nNode^2, TRUE, prob = c(1 - propPos, propPos))
    if(any(c("fixed", "mixed1") %in% m)){
      m1 <- m1 * runif(nNode, min(m1_range), max(m1_range))
      mmb1 <- rep(list(m1), nPerson)
    } else if("random0" %in% m){
      mmb1 <- lapply(1:nPerson, function(z){
        m10 <- rnorm(nNode, 0, m1SD)
        while(any(abs(m10) < min(m1_range)) | any(abs(m10) > max(m1_range))){
          m10 <- rnorm(nNode, 0, m1SD)
        }
        m10 <- m1 * m10
        return(m10)
      })
    } else {
      mmb1 <- lapply(1:nPerson, function(z) m1 * runif(nNode, min(m1_range), max(m1_range)))
    }
    if(any(c("fixed", "mixed2") %in% m)){
      m2 <- matrix(m2 * runif(nNode^2, min(m2_range), max(m2_range)), ncol = nNode)
      mmb2 <- rep(list(m2), nPerson)
    } else if("random0" %in% m){
      mmb2 <- lapply(1:nPerson, function(z){
        m20 <- rnorm(nNode^2, 0, m2SD)
        while(any(abs(m20) < min(m2_range)) | any(abs(m20) > max(m2_range))){
          m20 <- rnorm(nNode^2, 0, m2SD)
        }
        m20 <- c(m2) * m20
        return(matrix(m20, ncol = nNode))
      })
    } else {
      mmb2 <- lapply(1:nPerson, function(z){
        matrix(c(m2) * runif(nNode^2, min(m2_range), max(m2_range)), ncol = nNode)})
    }
    if(getM){return(list(m = mm[[1]], mb1 = mmb1[[1]], mb2 = mmb2[[1]]))}
  } else {
    mm <- mmb1 <- mmb2 <- NULL
  }
  mu_fixed <- rnorm(nNode, 0, fixedMuSD)
  beta_fixed <- rnorm(nTemporal, 0)
  beta_fixed[order(abs(beta_fixed))[1:round(nTemporal * GGMsparsity)]] <- 0
  mat <- matrix(0, nNode, nNode * lag)
  diag(mat) <- 1
  t1 <- 0
  beta_fixed[c(mat) == 1] <- runif(sum(c(mat) == 1), 0, 1)
  if(!is.null(m) & !modType %in% c('none', 'zero')){
    beta_fixed2 <- matrix(beta_fixed, nNode, nNode)
    while(ifelse(grepl('full', modType), !all(beta_fixed2[m2 != 0] == 0), any(beta_fixed2[m2 != 0] == 0))){
      beta_fixed <- matrix(rnorm(nTemporal, 0), nNode, nNode)
      allzeros <- round(nTemporal * GGMsparsity)
      if(grepl('full', modType)){
        beta_fixed[m2 != 0] <- 0
        allzeros <- allzeros - sum(beta_fixed == 0)
        bfix <- beta_fixed[beta_fixed != 0]
        beta_fixed[beta_fixed != 0][order(abs(bfix))[1:allzeros]] <- 0
      } else {
        bfix <- beta_fixed[m2 == 0]
        beta_fixed[m2 == 0][order(abs(bfix))[1:allzeros]] <- 0
      }
      mat <- switch(2 - (lag == 1), diag(nNode), cbind(diag(nNode), matrix(0, nNode, nNode * (lag - 1))))
      beta_fixed[c(mat) == 1] <- runif(sum(c(mat) == 1), 0, 1)
      beta_fixed2 <- matrix(beta_fixed, nNode, nNode)
    }
  }
  if(lag > 0){
    repeat{
      t1 <- t1 + 1
      if(all(skew == FALSE)){
        Pars <- mvtnorm::rmvnorm(nPerson, c(mu_fixed, beta_fixed), sigma = Omega)
      } else {
        if(is.numeric(skew)){skew <- rep(skew, nNode + nNode^2)}
        if(is.character(skew)){if(all(skew == 'random')){skew <- runif(nNode + nNode^2, -10000, 10000)}}
        if(isTRUE(skew) | length(skew) != (nNode + nNode^2)){skew <- rep(3, nNode + nNode^2)}
        Omega <- as.matrix(forceSymmetric(Omega))
        Pars <- sn::rmsn(nPerson, c(mu_fixed, beta_fixed), Omega, skew)
      }
      Mus <- Pars[, 1:nNode]
      Betas <- array(c(t(Pars[, -(1:nNode)])), c(nNode, nNode * lag, nPerson))
      if(lag > 1){
        under <- cbind(diag(nNode * (lag - 1)), matrix(0, nNode * (lag - 1), nNode))
        ev <- sapply(seq_len(nPerson), function(i){
          mat <- rbind(Betas[, , i], under)
          eigen(mat)$values
        })
      } else {
        ev <- sapply(seq_len(nPerson), function(i){eigen(Betas[, , i])$values})
      }
      allEV <- c(ev)
      if(all(Re(ev)^2 + Im(ev)^2 < 1)){
        if(nPerson == 1){
          Mus <- matrix(mu_fixed, ncol = nNode, nrow = 1)
          Betas <- array(matrix(beta_fixed, nNode, nNode), c(nNode, nNode, 1))
          Theta <- array(Theta_fixed, c(nNode, nNode, 1))
        }
        DataList <- lapply(1:nPerson, function(p){
          parms <- lapply(seq_len(lag), function(l) array(c(Betas[, , p]), c(nNode, nNode, lag))[, , l])
          if(lag > 0){
            if(onlyNets){
              res <- list(parms = parms, means = Mus[p, ], lags = seq_len(lag),
                          Nt = nTime[p], init = Mus[p, ], burnin = 100, 
                          residuals = Theta[, , p], m = mm[[p]], mb1 = mmb1[[p]],
                          mb2 = mmb2[[p]], mcenter = mcenter, skewErr = skewErr)
              return(res)
            }
            res <- trevSimulateVAR(parms, means = Mus[p, ], lags = seq_len(lag), 
                                   Nt = nTime[p], init = Mus[p, ], burnin = 100, 
                                   residuals = Theta[, , p], m = mm[[p]], 
                                   mb1 = mmb1[[p]], mb2 = mmb2[[p]], 
                                   mcenter = mcenter, skewErr = skewErr)
          } else {
            res <- mvtnorm::rmvnorm(nTime[p], Mus[p, ], Theta[, , p])
          }
          colnames(res) <- paste0("V", 1:nNode)
          if(ordinal & ordWithin){
            for(vv in 1:ncol(res)){
              tick <- ord <- 0
              while(length(unique(ord)) < minOrd){
                thresh <- switch(2 - is.null(thresholds), rnorm(nLevels - 1), thresholds[[vv]])
                ord <- as.numeric(cut(res[, vv], sort(c(-Inf, thresh, Inf))))
                tick <- tick + 1
                if(tick == 10 & !is.null(thresholds)){thresholds <- NULL}
                if(tick == 20){break}
              }
              res[, vv] <- ord
            }
          }
          res$ID <- p
          res
        })
        if(onlyNets){return(DataList[[1]])}
        Data <- do.call(rbind, DataList)
        if(!any(abs(Data[, 1:nNode]) > 100)){break}
      }
      beta_fixed <- beta_fixed * shrink_fixed
      D <- diag(sqrt(diag(Omega)))
      D[-(1:nNode), -(1:nNode)] <- shrink_deviation * D[-(1:nNode), -(1:nNode)]
      Omega <- D %*% cov2cor(Omega) %*% D
    }
  } else {
    Pars <- mvtnorm::rmvnorm(nPerson, mu_fixed, sigma = Omega)
    Mus <- Pars[, 1:nNode]
    Betas <- array(dim = c(0, 0, nPerson))
    DataList <- lapply(1:nPerson, function(p){
      res <- as.data.frame(mvtnorm::rmvnorm(nTime[p], Mus[p, ], Theta[, , p]))
      colnames(res) <- paste0("V", 1:nNode)
      res$ID <- p
      res
    })
    Data <- do.call(rbind, DataList)
  }
  model <- list(mu = trevModelArray(mean = mu_fixed, SD = mu_SD, subject = lapply(1:nrow(Mus), function(i) Mus[i, ])), 
                Beta = trevModelArray(mean = array(beta_fixed, c(nNode, nNode, lag)), 
                                      SD = array(sqrt(diag(Omega[-(1:nNode), -(1:nNode)])), c(nNode, nNode, lag)), 
                                      subject = lapply(1:nPerson, function(p) array(Betas[, , p], c(nNode, nNode, lag)))), 
                Omega_mu = trevModelCov(cov = trevModelArray(mean = Omega[1:nNode, 1:nNode])), 
                Theta = trevModelCov(cov = trevModelArray(mean = Theta_fixed, subject = lapply(1:nPerson, function(p) Theta[, , p]))), 
                Omega = trevModelCov(cov = trevModelArray(mean = Omega)))
  kappa <- corpcor::pseudoinverse(Theta_fixed)
  beta <- matrix(beta_fixed, nNode, nNode)
  mod2 <- list(fixedKappa = kappa, fixedPCC = pcor2(Theta_fixed), fixedBeta = beta, 
               fixedPDC = trevPDC(beta, kappa), between = pcor2(Omega[1:nNode, 1:nNode]))
  if(ordinal & !ordWithin){
    for(vv in 1:(ncol(Data) - 1)){
      tick <- ord <- 0
      while(length(unique(ord)) < minOrd){
        thresh <- switch(2 - is.null(thresholds), rnorm(nLevels - 1), thresholds[[vv]])
        ord <- as.numeric(cut(Data[, vv], sort(c(-Inf, thresh, Inf))))
        tick <- tick + 1
        if(tick == 10 & !is.null(thresholds)){thresholds <- NULL}
        if(tick == 20){break}
      }
      Data[, vv] <- ord
    }
  }
  Results <- list(data = Data, vars = paste0("V", 1:nNode), idvar = "ID", lag = lag, model = model)
  Results <- append(list(data = Data), append(mod2, Results[-1]))
  if(!is.null(m)){
    if('ordinal' %in% m){
      mm <- lapply(mm, function(z){
        ord <- c()
        while(length(unique(ord)) < minOrd){
          ord <- as.numeric(cut(z, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
        }
        return(ord)
      })
    }
    mm <- lapply(mm, function(z) z[-(1:100)])
    dat2 <- data.frame(Data[, 1:nNode], M = unlist(mm), ID = Data[, "ID"])
    mb1 <- if(any(c("fixed", "mixed1") %in% m)){mmb1[[1]]} else {Reduce("+", mmb1)/nPerson}
    mb2 <- if(any(c("fixed", "mixed2") %in% m)){mmb2[[1]]} else {Reduce("+", mmb2)/nPerson}
    Results$mm <- list(data = dat2, m = m, subjects = list(mb1 = mmb1, mb2 = mmb2), mb1 = mb1, mb2 = mb2)
    if(nPerson == 1){Results$mm$data <- Results$mm$data[, -ncol(Results$mm$data)]}
  }
  if(nPerson == 1){
    Results$data <- Results$data[, -ncol(Results$data)]
    Results$between <- NULL
    attr(Results, "simMLgvar") <- TRUE
    class(Results) <- c('list', 'simMLgvar')
  } else {
    class(Results) <- attr(Results, "mlVARsim") <- "mlVARsim"
  }
  return(Results)
}

##### trevSimulateVAR: core single-subject VAR sampler
trevSimulateVAR <- function(parms, means = 0, lags = 1, Nt = 100, init, 
                            residuals = 0.1, burnin, m = NULL, mb1 = NULL, 
                            mb2 = NULL, mcenter = TRUE, skewErr = FALSE){
  if(is.matrix(parms)){parms <- list(parms)} ##### START 3
  if(any(sapply(parms, function(x) length(unique(dim(x))) > 1))){
    stop("non-square graph detected.")
  }
  if(missing(burnin)){burnin <- min(round(Nt/2), 100)}
  Ni <- ncol(parms[[1]])
  if(length(means) == 1){means <- rep(means, Ni)}
  maxLag <- max(lags)
  if(length(residuals) == 1){
    residuals <- diag(residuals, Ni)
  } else if(length(residuals) == Ni){
    residuals <- diag(residuals)
  }
  if(!is.matrix(residuals) && ncol(residuals) != Ni && nrow(residuals) != Ni){
    stop("'residuals' is not a square matrix")
  }
  totTime <- Nt + burnin
  if(missing(init)){init <- matrix(0, maxLag, Ni)}
  Res <- matrix(NA, totTime, Ni)
  if(!is.null(m)){
    if(!is(m, "matrix")){
      if(mcenter){m <- m - mean(m)}
      m <- matrix(rep(m, ncol(Res)), ncol = ncol(Res))
    }
  }
  skewErr <- ifelse(is.numeric(skewErr), skewErr, ifelse(isTRUE(skewErr), 3, FALSE))
  Res[1:maxLag, ] <- init ##### STOP 3
  for(t in (maxLag + 1):(totTime)){
    if(!is.null(m)){
      x1 <- rowSums(do.call(cbind, lapply(seq_along(lags), function(i){
        parms[[i]] %*% (Res[t - lags[i], ] - means)})))
      x2 <- mb1 * m[t - 1, ]
      x3 <- rowSums(mb2 %*% ((Res[t - 1, ] - means) * m[t - 1, ]))
      Res[t, ] <- means + x1 + x2 + x3
    } else {
      Res[t, ] <- means + rowSums(do.call(cbind, lapply(seq_along(lags),function(i){
        parms[[i]] %*% (Res[t - lags[i], ] - means)})))
    }
    e <- switch(
      2 - is.numeric(skewErr), 
      sn::rmsn(1, rep(0, Ni), residuals, rep(skewErr, Ni)), 
      mvtnorm::rmvnorm(1, rep(0, Ni), residuals)
    )
    Res[t, ] <- Res[t, ] + e
  }
  out <- as.data.frame(Res[-(1:burnin), ])
  return(out)
}

##### simPcor: creates network matrices
simPcor <- function(Nvar, sparsity = 0.5, parRange = c(0.5, 1), constant = 1.5, 
                    propPos = 0.5, precision = FALSE, finalRange = NULL){
  trueKappa <- matrix(0, Nvar, Nvar)
  kupper <- upper.tri(trueKappa)
  klower <- lower.tri(trueKappa)
  totEdges <- sum(kupper)
  trueKappa[kupper][sample(seq_len(totEdges), round((1 - sparsity) * totEdges))] <- 1
  vals <- sample(c(-1, 1), totEdges, TRUE, prob = c(propPos, 1 - propPos)) * runif(totEdges, min(parRange), max(parRange))
  trueKappa[kupper] <- trueKappa[kupper] * vals
  trueKappa[klower] <- t(trueKappa)[klower]
  diag(trueKappa) <- constant * rowSums(abs(trueKappa))
  diag(trueKappa) <- ifelse(diag(trueKappa) == 0, 1, diag(trueKappa))
  trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
  trueKappa <- (trueKappa + t(trueKappa))/2
  if(!precision){trueKappa <- as.matrix(qgraph::wi2net(trueKappa))}
  if(!is.null(finalRange)){
    if(max(abs(trueKappa)) > max(finalRange)){
      trueKappa <- trueKappa/(max(abs(abs(trueKappa)))/max(finalRange))
    }
    if(min(abs(trueKappa[trueKappa != 0])) < min(finalRange)){
      wmin <- abs(trueKappa[trueKappa != 0]) < min(finalRange)
      while(min(abs(trueKappa[trueKappa != 0])) < min(finalRange)){
        trueKappa[trueKappa != 0][wmin] <- trueKappa[trueKappa != 0][wmin] * constant
      }
    }
  }
  return(trueKappa)
}

##### simCor: generates correlation matrices based on sparse precision matrices
simCor <- function(x, sparsity = .5, maxiter = 100, ...){
  ow <- getOption("warn"); options(warn = 2)
  try <- 0
  while(try < maxiter){
    out <- tryCatch({
      cov2cor(solve(diag(x) - simPcor(x, sparsity, ...)))}, 
      error = function(e){TRUE})
    if(!isTRUE(out)){break}
    try <- try + 1
  }
  options(warn = ow)
  if(try == maxiter){stop("Estimate failed to converge")}
  return(out)
}

##### trevModelArray: straight from mlVAR
trevModelArray <- function(mean, subject, SE, P, lower, upper, SD){
  Results <- list()
  if(missing(subject)){subject <- NULL}
  if(is.null(subject) && missing(mean)){stop("Either 'subject' or 'mean' must be used")}
  if(missing(mean) && !is.null(subject)){
    N <- length(subject)
    mean <- Reduce("+", subject)/N
  }
  dim <- switch(2 - is.null(dim(mean)), length(mean), dim(mean))
  if(is.null(dim)){dim <- length(mean)}
  if(missing(SD)){
    if(!is.null(subject)){
      N <- length(subject)
      SD <- Reduce("+", lapply(subject, function(x) (x - mean)^2))/(N - 1)
    } else {
      SD <- array(NA, dim)
    }
  }
  if(missing(P)){
    if(!missing(SE)){
      P <- 2 * (1 - pnorm(abs(mean/SE)))
    } else {
      P <- array(NA, dim)
    }
  }
  if(missing(SE)){SE <- array(NA, dim)}
  if(missing(lower)){
    if(!missing(SE)){
      lower <- mean - 1.959964 * SE
    } else {
      lower <- array(NA, dim)
    }
  }
  if(missing(upper)){
    if(!missing(SE)){
      upper <- mean + 1.959964 * SE
    } else {
      upper <- array(NA, dim)
    }
  }
  if(is.null(subject)){subject <- NULL}
  Results[["mean"]] <- mean
  Results[["SD"]] <- SD
  Results[["lower"]] <- lower
  Results[["upper"]] <- upper
  Results[["SE"]] <- SE
  Results[["P"]] <- P
  Results[["subject"]] <- subject
  class(Results) <- c("mlVARarray", "list")
  return(Results)
}

##### trevModelCov: straight from mlVAR
trevModelCov <- function(cov, cor, prec, pcor){
  if(missing(cov)){stop("'cov' can not be missing")}
  if(!missing(cor) && !is(cor, "mlVARarray")){stop("'cor' must be missing or an object of class 'mlVARarray'")}
  if(!missing(prec) && !is(prec, "mlVARarray")){stop("'prec' must be missing or an object of class 'mlVARarray'")}
  if(!missing(pcor) && !is(pcor, "mlVARarray")){stop("'pcor' must be missing or an object of class 'mlVARarray'")}
  cov2corNA <- function(x){
    if(any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(cov2cor(x))
    }
  }
  solveNA <- function(x){
    if(any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(corpcor::pseudoinverse(x))
    }
  }
  cor2pcorNA <- function(x){
    if(any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(corpcor::cor2pcor(x))
    }
  }
  if(missing(cor)){
    cor <- trevModelArray(mean = cov2corNA(cov[["mean"]]), 
                          subject = lapply(cov[["subject"]], cov2corNA))
  }
  if(missing(prec)){
    prec <- trevModelArray(mean = solveNA(cov[["mean"]]), 
                           subject = lapply(cov[["subject"]], solveNA))
  }
  if(missing(pcor)){
    pcor <- trevModelArray(mean = cor2pcorNA(cov[["mean"]]), 
                           subject = lapply(cov[["subject"]], cor2pcorNA))
  }
  Results <- list(cov = cov, cor = cor, prec = prec, pcor = pcor)
  class(Results) <- c("mlVarCov", "list")
  return(Results)
}

##### pcor2: creates partial correlation matrix
pcor2 <- function(x){
  x <- -corpcor::pseudoinverse(x)
  diag(x) <- -diag(x)
  x <- cov2cor(x) - diag(ncol(x))
  return((x + t(x))/2)
}


### ======================================================================== ###
### ======================= graphicalVAR Version =========================== ###
### ======================================================================== ###
##### simMLGVAR: simulate mlGVAR data
simMLGVAR <- function(nTime = 50, nPerson = 20, nVar = 3, propPositive = 0.5, 
                      kappaRange = NULL, betaRange = NULL, betweenRange = NULL, 
                      rewireWithin = 0, betweenVar = 1, withinVar = 0.25, 
                      temporalOffset = 2, m = NULL, m0 = 1, m1 = .3, 
                      m2 = .3, m2pos = .7, m2max = .3, mcenter = TRUE){
  kbb <- list(kappaRange = kappaRange, betaRange = betaRange, betweenRange = betweenRange)
  if(all(sapply(kbb, is.null))){
    kappaRange <- betaRange <- betweenRange <- c(.25, .5)
  } else if(any(sapply(kbb, is.null))){
    kbb[which(sapply(kbb, is.null))] <- c(.25, .5)
    invisible(lapply(1:3, function(i) assign(names(kbb)[i], kbb[[i]], pos = 1)))
  }
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
  if(!is.null(m)){
    if(is(m, "list")){
      mm <- m
    } else {
      mm <- list()
      for(j in 1:nPerson){
        if(tolower(m) == "ar"){
          if(m0 >= 1){m0 <- .3}
          mm[[j]] <- as.numeric(arima.sim(n = nTime + 100, list(ar = m0)))
        } else {
          mm[[j]] <- rnorm(nTime + 100, mean = 0, sd = m0)
        }
      }
    }
    if(length(m1) != nVar){m1 <- runif(nVar, -m1, m1)}
    mmb1 <- rep(list(m1), nPerson)
    if(!is(m2, "matrix")){
      m2 <- matrix(sample(c(0, 1), nVar^2, TRUE, prob = c(1 - m2, m2)), nVar, nVar)
      m2 <- m2 * sample(c(-1, 1), nVar^2, TRUE, prob = c(1 - m2pos, m2pos)) * runif(nVar^2, 0, m2max)
    }
    mmb2 <- rep(list(m2), nPerson)
  } else {
    mm <- mmb1 <- mmb2 <- NULL
  }
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
    Data <- as.data.frame(
      trevSimGVAR(nTime, beta, kappa, mean = Means[i, ], 
                   m = mm[[i]], mb1 = mmb1[[i]], mb2 = mmb2[[i]], 
                   mcenter = mcenter))
    Data$ID <- i
    return(list(kappa = kappa, beta = beta, PCC = trevPCC(kappa), 
                PDC = trevPDC(beta, kappa), data = Data))
  })
  fixedKappa <- Reduce("+", lapply(SubjectData, "[[", "kappa"))/nPerson
  fixedBeta <- Reduce("+", lapply(SubjectData, "[[", "beta"))/nPerson
  allData <- do.call(rbind, lapply(SubjectData, "[[", "data"))
  Results <- list(data = allData, fixedKappa = fixedKappa, 
                  fixedPCC = trevPCC(fixedKappa), fixedBeta = fixedBeta, 
                  fixedPDC = trevPDC(fixedBeta, fixedKappa), between = Sigma,
                  Omega = trueBetween, means = Means, personData = SubjectData, 
                  idvar = "ID", vars = names(allData)[names(allData) != "ID"])
  if(!is.null(m)){
    mm <- lapply(mm, function(z) z[-(1:100)])
    dat2 <- data.frame(allData[, 1:nVar], M = unlist(mm), ID = allData[, "ID"])
    Results$mm <- list(data = dat2, mb1 = mmb1[[1]], mb2 = mmb2[[1]])
  }
  if(nPerson == 1){
    Results[c("personData", "idvar")] <- NULL
    Results$data <- Results$data[, -which(colnames(Results$data) == "ID")]
    if(!is.null(m)){Results$mm$data <- Results$mm$data[, -ncol(Results$mm$data)]}
  }
  attr(Results, "simMLgvar") <- TRUE
  return(Results)
}

##### trevSimGVAR
trevSimGVAR <- function(nTime = 100, beta = NULL, kappa = NULL, mean = 0, 
                        init = NULL, warmup = 100, lbound = -Inf, ubound = Inf, 
                        m = NULL, mb1 = NULL, mb2 = NULL, mcenter = TRUE, ...){
  ubound <- Inf; lbound <- -Inf; warmup <- 100; init <- NULL; mm <- NULL
  if(length(mean) == 1){mean <- rep(mean, ncol(kappa))}
  if(length(lbound) == 1){lbound <- rep(lbound, ncol(kappa))}
  if(length(ubound) == 1){ubound <- rep(ubound, ncol(kappa))}
  if(is.null(init)){init <- mean}
  stopifnot(!is.null(beta))
  stopifnot(!is.null(kappa))
  Nvar <- ncol(kappa)
  init <- rep(init, length = Nvar)
  totTime <- nTime + warmup
  Data <- t(matrix(init, Nvar, totTime))
  if(!is.null(m)){
    if(!is.numeric(m)){
      mm <- mlGVARsim(nTime = nTime, nPerson = 1, nNode = Nvar, 
                      m = m, mcenter = mcenter, getM = TRUE, ...)
      m <- mm$m
      mb1 <- mm$mb1
      mb2 <- mm$mb2
    }
    if(!is(m, "matrix")){
      if(mcenter){m <- m - mean(m)}
      m <- matrix(rep(m, ncol(Data)), ncol = ncol(Data))
    }
  }
  Sigma <- solve(kappa)
  for(t in 2:totTime){
    if(!is.null(m)){
      x1 <- t(beta %*% (Data[t - 1, ] - mean))
      x2 <- mb1 * m[t - 1, ]
      x3 <- t(mb2 %*% ((Data[t - 1, ] - mean) * m[t - 1, ]))
      x4 <- mvtnorm::rmvnorm(1, rep(0, Nvar), Sigma)
      Data[t, ] <- mean + x1 + x2 + x3 + x4
    } else {
      Data[t, ] <- mean + t(beta %*% (Data[t - 1, ] - mean)) + mvtnorm::rmvnorm(1, rep(0, Nvar), Sigma)
    }
    Data[t, ] <- ifelse(Data[t, ] < lbound, lbound, Data[t, ])
    Data[t, ] <- ifelse(Data[t, ] > ubound, ubound, Data[t, ])
  }
  if(!is.null(mm)){Data <- cbind(Data, mm$m)}
  output <- Data[-seq_len(warmup), , drop = FALSE]
  if(!is.null(mm)){output <- list(data = output, mb1 = mb1, mb2 = mb2)}
  return(output)
}

##### trevRewire
trevRewire <- function(x, p, directed){
  if(missing(directed)){directed <- !all(x == t(x))}
  ind <- if(directed){diag(1, ncol(x)) != 1} else {upper.tri(x)}
  curEdges <- which(x != 0 & ind, arr.ind = TRUE)
  toRewire <- which(runif(nrow(curEdges)) < p)
  if(any(x == 0)){
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
  }
  return(x)
}

##### trevPDC
trevPDC <- function(beta, kappa){
  if(ncol(beta) == nrow(beta) + 1){beta <- beta[, -1, drop = FALSE]}
  sigma <- solve(kappa)
  t(beta/sqrt(diag(sigma) %o% diag(kappa) + beta^2))
}

##### trevPCC
trevPCC <- function(x){
  x <- -cov2cor(x)
  diag(x) <- 0
  return(as.matrix(Matrix::forceSymmetric(x)))
}

##### randomGVAR: simulate the matrices only
randomGVAR <- function(Nvar, probKappaEdge = 0.1, probBetaEdge = 0.1, 
                       probKappaPositive = 0.5, probBetaPositive = 0.5, 
                       maxtry = 10, kappaConstant = 1.1){
  try <- 0
  repeat {
    kappaRange = c(0.5, 1)
    trueKappa <- matrix(0, Nvar, Nvar)
    trueKappa[upper.tri(trueKappa)] <- sample(c(0, 1), sum(upper.tri(trueKappa)), TRUE, prob = c(1 - probKappaEdge, probKappaEdge))
    trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] * sample(
      c(1, -1), sum(upper.tri(trueKappa)), TRUE, prob = c(1 - probKappaPositive, probKappaPositive)) * runif(
        sum(upper.tri(trueKappa)), min(kappaRange), max(kappaRange))
    trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]
    diag(trueKappa) <- kappaConstant * rowSums(abs(trueKappa))
    diag(trueKappa) <- ifelse(diag(trueKappa) == 0, 1, diag(trueKappa))
    trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
    trueKappa <- (trueKappa + t(trueKappa))/2
    Vmin <- min(abs(trueKappa[trueKappa != 0]))
    trueBeta <- matrix(sample(c(0, 1), Nvar^2, TRUE, prob = c(1 - probBetaEdge, probBetaEdge)), Nvar, Nvar)
    trueBeta <- trueBeta * sample(c(-1, 1), Nvar^2, TRUE, prob = c(1 - probBetaPositive, probBetaPositive)) * runif(Nvar^2, Vmin, 1)
    diag(trueBeta) <- Vmin
    evK <- eigen(trueKappa)$values
    evB <- eigen(trueBeta)$values
    while(any(Re(evB)^2 + Im(evB)^2 > 1)){
      trueBeta <- 0.95 * trueBeta
      evB <- eigen(trueBeta)$values
    }
    if(all(evK > 0) & all(Re(evB)^2 + Im(evB)^2 < 1)){break}
    try <- try + 1
    if(try > maxtry){stop("Maximum number of tries reached.")}
  }
  Res <- list(kappa = trueKappa, beta = trueBeta, PCC = trevPCC(trueKappa), 
              PDC = trevPDC(trueBeta, trueKappa))
  class(Res) <- "gVARmodel"
  return(Res)
}

##### cor0: deal with sd == 0
cor0 <- function(x, y){
  if(all(is.na(x)) | all(is.na(y))){return(NA)}
  if(sd(x, na.rm = TRUE) == 0 | sd(y, na.rm = TRUE) == 0){return(0)}
  return(cor(x, y, use = 'pairwise.complete.obs'))
}


### ======================================================================== ###
### ============================ CROSS-SECTIONAL =========================== ###
### ======================================================================== ###
##### mgmSim: simulates data with interactions from mgm
mgmSim <- function(N = 500, p = 5, prop = .2, bmin = .1, bmax = .2, exo = NULL,
                   propPos = .5, omega_mu = 1, muSD = 1, maxiter = 10, mgm2 = FALSE){
  if(mgm2){source('editing/extras/mgmsampler2.R')}
  ff <- t(combn(p, 2))
  n <- nrow(ff)
  k <- floor(n * prop)
  edges <- ff[sample(1:n, k, replace = FALSE), ]
  if(!is.null(exo) & mgm2){
    while(exo %in% edges){
      edges <- ff[sample(1:n, k, replace = FALSE), ]
    }
  }
  factors <- list(edges, matrix(c(2, 3, 4), ncol = 3, byrow = T))
  ints0 <- runif(k, bmin, bmax) * sample(
    c(-1, 1), k, TRUE, prob = c(propPos, 1 - propPos))
  interactions <- list(lapply(as.list(ints0), array, dim = c(1, 1)), 
                       list(array(.15, dim = c(1, 1, 1))))
  intercepts <- as.list(rnorm(p, 0, muSD))
  sds <- if(length(muSD) == 1){rep(muSD, p)} else {muSD[1:k]}
  if(any(is.na(sds))){sds[is.na(sds)] <- 1}
  oldw <- getOption("warn"); options(warn = 2)
  out <- iter <- TRUE
  while(isTRUE(out)){
    out <- tryCatch({
      if(mgm2){
        mgmsampler2(factors = factors, interactions = interactions, 
                    thresholds = intercepts, sds = sds, type = rep("g", p), 
                    level = rep(1, p), N = N, nIter = 100, pbar = TRUE,
                    exogenous = exo)
      } else {
        mgm::mgmsampler(factors = factors, interactions = interactions, 
                        thresholds = intercepts, sds = sds, type = rep("g", p), 
                        level = rep(1, p), N = N, nIter = 100, pbar = TRUE)
      }
    }, error = function(e){TRUE})
    if(isTRUE(out)){
      interactions[[2]] <- lapply(interactions[[2]], function(z){
        z[, , 1] <- z[, , 1] * .9
        return(z)
      })
    }
    iter <- as.numeric(iter) + 1
    if(iter > maxiter){break}
  }
  out$chains <- NULL
  out$trueNet <- diag(0, p)
  for(i in seq_len(k)){out$trueNet[edges[i, 1], edges[i, 2]] <- ints0[i]}
  out$trueNet <- out$trueNet + t(out$trueNet)
  class(out) <- c('list', 'mgmSim')
  options(warn = oldw)
  return(out)
}

##### mgm2: just quick 'mgm' wrapper
mgm2 <- function(data, m = NULL, sel = 'EBIC', gamma = .25, ...){
  if(identical(m, FALSE)){m <- NULL}
  type <- rep('g', ncol(data))
  level <- rep(1, ncol(data))
  mgm::mgm(data = data, type = type, level = level, lambdaSel = sel, 
           lambdaGam = gamma, moderators = m, ...)
}


### ======================================================================== ###
### ======================================================================== ###
setupModnets <- function(folders = '.', rm = FALSE){
  if(any(folders %in% c('.', ''))){folders <- '.'}
  fold <- "~/C:\\ Desktop/COMPS/METHODS/CODE/modnets/"
  mods <- paste0(c("functions", "ggm", "centrality", "sim", "mlGVAR", 
                   "simGVAR", "penalized", "power", "plots"), ".R")
  for(j in seq_along(folders)){
    if(!folders[j] %in% dir() & folders[j] != '.'){
      stopifnot(length(folders) == 1)
      folders[j] <- '.'
    }
    if(isTRUE(rm)){system(paste0('rm -rf ', folders[j], '/modnets'))}
    system(paste0('mkdir ', folders[j], '/modnets'))
    for(k in seq_along(mods)){system(paste0("cp ", fold, mods[k], " ./", folders[j], "/modnets/"))}
    f0 <- readLines(paste0(folders[j], "/modnets/functions.R"))
    f1 <- which(f0 == "setwd('~/C: Desktop/COMPS/METHODS/CODE/modnets')")
    f2 <- which(f0 == "files <- paste0('./', c('ggm', 'centrality', 'sim', 'mlGVAR', 'simGVAR', 'penalized', 'power', 'plots'),'.R')")
    f0[f1] <- paste0("#", f0[f1])
    f0[f2] <- gsub("'./'", "'modnets/'", f0[f2])
    f0 <- paste0(paste0(f0, collapse = "\n"), "\n")
    cat(f0, file = paste0(folders[j], "/modnets/functions.R"))
  }
}

