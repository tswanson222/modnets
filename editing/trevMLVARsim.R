################################################################################
mlGVARsim <- function(nTime = 50, nPerson = 10, nNode = 3, lag = 1, 
                      thetaVar = NULL, mu_SD = NULL, init_beta_SD = NULL, 
                      fixedMuSD = 1, manual = FALSE, shrink_fixed = 0.9, 
                      shrink_deviation = 0.9, contemporaneous = "wishart",
                      GGMsparsity = .5, m = NULL, m0 = 1, m1 = .7, m2 = .25, 
                      m1_range = NULL, m2_range = NULL, propPos = .5,
                      mcenter = TRUE, skew = FALSE, maxiter = 100, 
                      getM = FALSE, m1SD = .1, m2SD = .1){
  if(is.null(thetaVar)){thetaVar <- rep(1, nNode)} ##### START
  if(is.null(mu_SD)){mu_SD <- c(1, 1)}
  if(is.null(init_beta_SD)){init_beta_SD <- c(.1, 1)}
  if(is.null(m1_range)){m1_range <- c(.1, .4)}
  if(is.null(m2_range)){m2_range <- c(.1, .3)}
  if(length(nTime) == 1){nTime <- rep(nTime, nPerson)}
  DF_theta <- nNode * 2
  nTemporal <- nNode^2 * lag
  contemporaneous <- match.arg(contemporaneous, c("wishart", "randomGGM", "fixed"))
  Omega_mu <- simPcor(nNode, GGMsparsity, maxiter)
  Omega_Beta <- simPcor(nTemporal, GGMsparsity, maxiter)
  Omega <- rbind(cbind(Omega_mu, matrix(0, nNode, nTemporal)), 
                 cbind(matrix(0, nTemporal, nNode), Omega_Beta))
  SD <- runif(nNode + nTemporal, c(rep(mu_SD[1], nNode), rep(init_beta_SD[1], nNode)), 
              c(rep(mu_SD[2], nNode), rep(init_beta_SD[2], nNode)))
  Omega <- diag(SD) %*% Omega %*% diag(SD)
  if(contemporaneous == "wishart"){
    Theta_fixed <- simPcor(nNode, GGMsparsity, maxiter)
    Theta_fixed <- diag(sqrt(thetaVar)) %*% Theta_fixed %*% diag(sqrt(thetaVar))
    Theta <- rWishart(nPerson, DF_theta, Theta_fixed/DF_theta)
  } else {
    if(contemporaneous == "randomGGM"){
      Theta <- lapply(1:nPerson, function(x) simPcor(nNode, GGMsparsity, maxiter))
      Theta <- do.call(abind::abind, c(Theta, along = 3))
      Theta_fixed <- apply(Theta, 1:2, mean)
    } else {
      Theta_fixed <- simPcor(nNode, GGMsparsity, maxiter)
      Theta <- lapply(1:nPerson, function(x) Theta_fixed)
      Theta <- do.call(abind::abind, c(Theta, along = 3))
    }
  }
  if(!is.null(m)){
    if(isTRUE(m)){m <- "random"}
    m <- match.arg(tolower(m), c(
      "fixed", "random", "mixed1", "mixed2", "ar", "binary", 
      "skewed", "mlm", "random0"), several.ok = TRUE)
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
    if("mlm" %in% m){ # NOT WORKING YET
      m10 <- rnorm(nNode, 0, mean(m1_range))
      while(any(abs(m10) < min(m1_range)) | any(abs(m10) > max(m1_range))){
        m10 <- rnorm(nNode, 0, mean(m1_range))
      }
      m20 <- rnorm(nNode^2, 0, mean(m2_range))
      while(any(abs(m20) < min(m2_range)) | any(abs(m20) > max(m2_range))){
        m20 <- rnorm(nNode^2, 0, mean(m2_range))
      }
      m10[order(abs(m10))][1:round(nNode * (1 - m1))] <- 0
      m20[order(abs(m20))][1:round((nNode^2) * (1 - m2))] <- 0
      Omega_m1 <- simPcor(nNode, (1 - m1), maxiter)
      Omega_m2 <- simPcor(nNode^2, (1 - m2), maxiter)
      mmb1 <- mvtnorm::rmvnorm(nPerson, m10, Omega_m1)
      mmb2 <- mvtnorm::rmvnorm(nPerson, m20, Omega_m2) # NOPE......
      mmb1 <- lapply(apply(mmb1, 1, list), unlist)
      mmb2 <- lapply(lapply(apply(mmb2, 1, list), unlist), matrix, nNode, nNode)
    } else {
      m1 <- sample(c(0, 1), nNode, TRUE, prob = c(1 - m1, m1))
      m1 <- m1 * sample(c(-1, 1), nNode, TRUE, prob = c(1 - propPos, propPos))
      m2 <- matrix(sample(c(0, 1), nNode^2, TRUE, prob = c(1 - m2, m2)), nNode, nNode)
      m2 <- m2 * sample(c(-1, 1), nNode^2, TRUE, prob = c(1 - propPos, propPos))
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
    }
    if(getM){return(list(m = mm[[1]], mb1 = mmb1[[1]], mb2 = mmb2[[1]]))}
  } else {
    mm <- mmb1 <- mmb2 <- NULL
  }
  mu_fixed <- rnorm(nNode, 0, fixedMuSD)
  beta_fixed <- rnorm(nTemporal, 0)
  beta_fixed[order(abs(beta_fixed))[1:round(nTemporal/2)]] <- 0
  mat <- matrix(0, nNode, nNode * lag)
  diag(mat) <- 1
  t1 <- 0
  beta_fixed[c(mat) == 1] <- runif(sum(c(mat) == 1), 0, 1) ##### STOP
  if(lag > 0){
    repeat { ##### START 2
      t1 <- t1 + 1
      if(all(skew == FALSE)){
        Pars <- mvtnorm::rmvnorm(nPerson, c(mu_fixed, beta_fixed), sigma = Omega)
      } else {
        if(isTRUE(skew) | length(skew) != (nNode + nNode^2)){skew <- runif(nNode + nNode^2, -10000, 10000)}
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
        if(manual){break}
        if(nPerson == 1){
          Mus <- matrix(mu_fixed, ncol = nNode, nrow = 1)
          Betas <- array(matrix(beta_fixed, nNode, nNode), c(nNode, nNode, 1))
          Theta <- array(Theta_fixed, c(nNode, nNode, 1))
        }
        DataList <- lapply(1:nPerson, function(p){
          parms <- lapply(seq_len(lag), function(l) array(c(Betas[, , p]), c(nNode, nNode, lag))[, , l])
          if(lag > 0){
            if(manual){
              pars("trevSimulateVAR", parms = parms, means = Mus[p, ], 
                   lags = 1, Nt = nTime[p], init = Mus[p, ], burnin = 100, 
                   residuals = Theta[, , p], m = mm[[p]], 
                   mb1 = mmb1[[p]], mb2 = mmb2[[p]])
            }
            res <- trevSimulateVAR(parms, means = Mus[p, ], lags = seq_len(lag), 
                                   Nt = nTime[p], init = Mus[p, ], burnin = 100, 
                                   residuals = Theta[, , p], m = mm[[p]], 
                                   mb1 = mmb1[[p]], mb2 = mmb2[[p]], 
                                   mcenter = mcenter)
          } else {
            res <- mvtnorm::rmvnorm(nTime[p], Mus[p, ], Theta[, , p])
          }
          colnames(res) <- paste0("V", 1:nNode)
          res$ID <- p
          res
        })
        Data <- do.call(rbind, DataList)
        if(!any(abs(Data[, 1:nNode]) > 100)){break}
      }
      beta_fixed <- beta_fixed * shrink_fixed
      D <- diag(sqrt(diag(Omega)))
      D[-(1:nNode), -(1:nNode)] <- shrink_deviation * D[-(1:nNode), -(1:nNode)]
      Omega <- D %*% cov2cor(Omega) %*% D
    } ##### STOP 2
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
                Beta = trevModelArray(
                  mean = array(beta_fixed, c(nNode, nNode, lag)), 
                  SD = array(sqrt(diag(Omega[-(1:nNode), -(1:nNode)])), c(nNode, nNode, lag)), 
                  subject = lapply(1:nPerson, function(p) array(Betas[, , p], c(nNode, nNode, lag)))), 
                Omega_mu = trevModelCov(cov = trevModelArray(mean = Omega[1:nNode, 1:nNode])), 
                Theta = trevModelCov(cov = trevModelArray(mean = Theta_fixed, subject = lapply(1:nPerson, function(p) Theta[, , p]))), 
                Omega = trevModelCov(cov = trevModelArray(mean = Omega)))
  kappa <- corpcor::pseudoinverse(Theta_fixed)
  beta <- matrix(beta_fixed, nNode, nNode)
  mod2 <- list(fixedKappa = kappa, fixedPCC = pcor2(Theta_fixed), fixedBeta = beta, 
               fixedPDC = trevPDC(beta, kappa), between = pcor2(Omega[1:nNode, 1:nNode]))
  Results <- list(data = Data, vars = paste0("V", 1:nNode), idvar = "ID", lag = lag, model = model)
  Results <- append(list(data = Data), append(mod2, Results[-1]))
  if(!is.null(m)){
    mm <- lapply(mm, function(z) z[-(1:100)])
    dat2 <- data.frame(Data[, 1:nNode], M = unlist(mm), ID = Data[, "ID"])
    mb1 <- if(any(c("fixed", "mixed1") %in% m)){mmb1[[1]]} else {Reduce("+", mmb1)/nPerson}
    mb2 <- if(any(c("fixed", "mixed2") %in% m)){mmb2[[1]]} else {Reduce("+", mmb2)/nPerson}
    Results$mm <- list(data = dat2, m = m, subjects = list(mb1 = mmb1, mb2 = mmb2),
                       mb1 = mb1, mb2 = mb2)
    if(nPerson == 1){Results$mm$data <- Results$mm$data[, -ncol(Results$mm$data)]}
  }
  if(nPerson == 1){
    Results$data <- Results$data[, -ncol(Results$data)]
    Results$between <- NULL
    attr(Results, "simMLgvar") <- TRUE
  } else {
    class(Results) <- attr(Results, "mlVARsim") <- "mlVARsim"
  }
  return(Results)
}


################################################################################
trevSimulateVAR <- function(parms, means = 0, lags = 1, Nt = 100, init, 
                            residuals = 0.1, burnin, m = NULL, mb1 = NULL, 
                            mb2 = NULL, mcenter = TRUE){
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
  Res[1:maxLag, ] <- init ##### STOP 3
  for(t in (maxLag + 1):(totTime)){
    if(!is.null(m)){
      x1 <- rowSums(do.call(cbind, lapply(seq_along(lags), function(i){
        parms[[i]] %*% (Res[t - lags[i], ] - means)})))
      x2 <- mb1 * m[t - 1, ]
      x3 <- rowSums(mb2 %*% ((Res[t - 1, ] - means) * m[t - 1, ]))
      x4 <- mvtnorm::rmvnorm(1, rep(0, Ni), residuals)
      Res[t, ] <- means + x1 + x2 + x3 + x4
    } else {
      Res[t, ] <- means + rowSums(do.call(cbind, lapply(seq_along(lags), function(i){
        parms[[i]] %*% (Res[t - lags[i], ] - means)}))) + mvtnorm::rmvnorm(1, rep(0, Ni), residuals)
    }
  }
  return(as.data.frame(Res[-(1:burnin), ]))
}

################################################################################
trevSimGraph <- function(Nvar, sparsity = 0.5, parRange = c(0.5, 1), 
                         constant = 1.1, propPositive = 0.5){
  trueKappa <- matrix(0, Nvar, Nvar)
  totEdges <- sum(upper.tri(trueKappa))
  nEdges <- round((1 - sparsity) * totEdges)
  inclEdges <- sample(seq_len(totEdges), nEdges)
  trueKappa[upper.tri(trueKappa)][inclEdges] <- 1
  trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] * sample(
    c(-1, 1), sum(upper.tri(trueKappa)), TRUE, 
    prob = c(propPositive, 1 - propPositive)) * runif(
      sum(upper.tri(trueKappa)), min(parRange), max(parRange))
  trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]
  diag(trueKappa) <- constant * rowSums(abs(trueKappa))
  diag(trueKappa) <- ifelse(diag(trueKappa) == 0, 1, diag(trueKappa))
  trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
  trueKappa <- (trueKappa + t(trueKappa))/2
  return(as.matrix(qgraph::wi2net(trueKappa)))
}

simPcor <- function(x, sparsity = .5, maxiter = 100, ...){
  ow <- getOption("warn"); options(warn = 2)
  try <- 0
  while(try < maxiter){
    out <- tryCatch({
      cov2cor(solve(diag(x) - trevSimGraph(x, sparsity, ...)))}, 
      error = function(e){TRUE})
    if(!isTRUE(out)){break}
    try <- try + 1
  }
  options(warn = ow)
  if(try == maxiter){stop("Estimate failed to converge")}
  return(out)
}

################################################################################
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

################################################################################
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

################################################################################
runSim <- function(x = 1){
  assign("manual", TRUE, envir = .GlobalEnv)
  ff <- readLines("trevMLVARsim.R")
  f1.1 <- which(ff == "  if(is.null(thetaVar)){thetaVar <- rep(1, nNode)} ##### START")
  f1.2 <- which(ff == "  beta_fixed[c(mat) == 1] <- runif(sum(c(mat) == 1), 0, 1) ##### STOP")
  f2.1 <- which(ff == "    repeat { ##### START 2")
  f2.2 <- which(ff == "    } ##### STOP 2")
  f3.1 <- which(ff == "  if(is.matrix(parms)){parms <- list(parms)} ##### START 3")
  f3.2 <- which(ff == "  Res[1:maxLag, ] <- init ##### STOP 3")
  if(x == 1){
    slines("trevMLVARsim.R", c(f1.1:f1.2, f2.1:f2.2))
  } else if(x == 2){
    slines("trevMLVARsim.R", f3.1:f3.2)
  }
}

pcor2 <- function(x){
  x <- -corpcor::pseudoinverse(x)
  diag(x) <- -diag(x)
  x <- cov2cor(x) - diag(ncol(x))
  return((x + t(x))/2)
}

################################################################################
################################ CROSS-SECTIONAL ###############################
################################################################################
mgmSim <- function(N = 500, p = 10, prop = .2, bmin = .1, bmax = .2, 
                   propPos = .5, omega_mu = 1, muSD = 1, maxiter = 10){
  ff <- t(combn(p, 2))
  n <- nrow(ff)
  k <- floor(n * prop)
  edges <- ff[sample(1:n, k, replace = FALSE), ]
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
      mgm::mgmsampler(factors = factors, interactions = interactions, 
                      thresholds = intercepts, sds = sds, type = rep("g", p), 
                      level = rep(1, p), N = N, nIter = 100, pbar = TRUE)
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
  options(warn = oldw)
  return(out)
}

