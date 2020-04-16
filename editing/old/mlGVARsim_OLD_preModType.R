##### mlGVARsim: main workhorse for simulating VAR and mlGVAR data
mlGVARsim <- function(nTime = 50, nPerson = 10, nNode = 3, m = NULL, m2 = .25, m1 = .7, 
                      m0 = 1, lag = 1, thetaVar = NULL, mu_SD = NULL, init_beta_SD = NULL, 
                      fixedMuSD = 1, manual = FALSE, shrink_fixed = 0.9, propPos = .5,
                      m1SD = .1, m1_range = NULL, m2SD = .1, shrink_deviation = 0.9, 
                      m2_range = NULL, contemporaneous = "wishart", GGMsparsity = .5,
                      mcenter = TRUE, getM = FALSE, skew = FALSE, skewErr = FALSE,
                      ordinal = FALSE, nLevels = 5, ordWithin = TRUE, minOrd = 3, 
                      thresholds = NULL, mseed = NULL, keepNets = FALSE, modType = 'none'){
  modType <- match.arg(tolower(modType), c('none', 'full', 'partial', 'full2', 'partial2', 'zero'))
  if(identical(m, FALSE)){m <- NULL}
  if(minOrd < 2){minOrd <- 2}
  if(is.numeric(ordinal)){nLevels <- ordinal; ordinal <- TRUE}
  if(is.null(thetaVar)){thetaVar <- rep(1, nNode)} ##### START
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
  } else {
    if(contemporaneous == "randomGGM"){
      Theta <- lapply(1:nPerson, function(x) simCor(nNode, GGMsparsity))
      Theta <- do.call(abind::abind, c(Theta, along = 3))
      Theta_fixed <- apply(Theta, 1:2, mean)
    } else {
      Theta_fixed <- simCor(nNode, GGMsparsity)
      Theta <- lapply(1:nPerson, function(x) Theta_fixed)
      Theta <- do.call(abind::abind, c(Theta, along = 3))
    }
  }
  if(!is.null(m)){
    if(isTRUE(m)){m <- "random"}
    m <- match.arg(tolower(m), c(
      "fixed", "random", "mixed1", "mixed2", "ar", "binary", 
      "skewed", "mlm", "random0", "ordinal"), several.ok = TRUE)
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
    if("ordinal" %in% m & FALSE){ # Shut down, moved near end of function
      mm <- lapply(mm, function(z){
        ord <- c()
        while(length(unique(ord)) < minOrd){
          ord <- as.numeric(cut(z, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
        }
        return(ord)
      })
    }
    if("mlm" %in% m & FALSE){ # NOT WORKING YET
      #m10 <- rnorm(nNode, 0, mean(m1_range))
      #while(any(abs(m10) < min(m1_range)) | any(abs(m10) > max(m1_range))){
      #  m10 <- rnorm(nNode, 0, mean(m1_range))
      #}
      #m20 <- rnorm(nNode^2, 0, mean(m2_range))
      #while(any(abs(m20) < min(m2_range)) | any(abs(m20) > max(m2_range))){
      #  m20 <- rnorm(nNode^2, 0, mean(m2_range))
      #}
      #m10[order(abs(m10))][1:round(nNode * (1 - m1))] <- 0
      #m20[order(abs(m20))][1:round((nNode^2) * (1 - m2))] <- 0
      #Omega_m1 <- simCor(nNode, (1 - m1))
      #Omega_m2 <- simCor(nNode^2, (1 - m2))
      #mmb1 <- mvtnorm::rmvnorm(nPerson, m10, Omega_m1)
      #mmb2 <- mvtnorm::rmvnorm(nPerson, m20, Omega_m2) # NOPE......
      #mmb1 <- lapply(apply(mmb1, 1, list), unlist)
      #mmb2 <- lapply(lapply(apply(mmb2, 1, list), unlist), matrix, nNode, nNode)
    } else {
      if(isTRUE(m1) | is.null(m1)){m1 <- .7}
      if(m1 >= 0 & m1 < 1){
        m1 <- sample(c(0, 1), nNode, TRUE, prob = c(1 - m1, m1))
      } else if(m1 >= 1){
        m10 <- numeric(nNode)
        m10[sample(1:nNode, ifelse(m1 > nNode, nNode, round(m1)))] <- 1
        m1 <- m10
      }
      if(m2 >= 0 & m2 < 1){
        m20 <- 0
        while(all(m20 == 0)){
          m20 <- matrix(sample(c(0, 1), nNode^2, TRUE, prob = c(1 - m2, m2)), nNode, nNode)
        }
        m2 <- m20
      } else if(m2 >= 1){
        m20 <- numeric(nNode^2)
        m20[sample(1:(nNode^2), ifelse(m2 > nNode^2, nNode^2, round(m2)))] <- 1
        m2 <- m20
      }
      m1 <- m1 * sample(c(-1, 1), nNode, TRUE, prob = c(1 - propPos, propPos))
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
  if(!is.null(m) & modType != 'none'){
    cat('hooray')
  }
  if(lag > 0){
    repeat { ##### START 2
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
            if(keepNets){
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
        if(keepNets){return(DataList[[1]])}
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
