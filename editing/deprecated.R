################################################################################
################################### USELESS ###################################
################################################################################
##### datMat (DEPRECATED): Create a list of matrices for fitting lagged models
datMat <- function(data, type = "g", m = NULL, lags = 1, scale = FALSE, 
                   center = TRUE, name = TRUE, full = TRUE){
  if(length(type) == 1 & ncol(data) > 1){type <- rep(type, ncol(data))}
  if(!is.null(m)){if(m != 1){m <- NULL}}
  if(!all(c(center, scale) == TRUE)){
    data <- scale(data, center, scale)
    scale <- FALSE
  }
  dat <- list(); Y <- matrix(NA, ncol = ncol(data), nrow = nrow(data) - lags)
  for(i in 1:ncol(data)){
    dat[[i]] <- setup(data = data, type = type, y = i, lags = lags, scale = scale)
    Y[,i] <- dat[[i]][,1]
  }
  if(!is.null(colnames(data))){colnames(Y) <- paste0(colnames(data), ".y")}
  rownames(Y) <- 1:nrow(Y)
  X <- interactionMatrix(data = dat[[1]], y = 1, type = type, moderators = m, lags = lags)
  if(name == TRUE & is.null(m)){colnames(X) <- colnames(data)}
  mats <- list(Y = Y, X = X)
  if(full == TRUE){
    mats <- list(Y = Y, X = X, full = data.frame(do.call(cbind, mats)))
    p <- ncol(mats$Y)
    if(!is.null(m)){
      colnames(mats$full)[((p * 2) + 1):ncol(mats$full)] <- colnames(X)[(p + 1):ncol(X)]
    }
  }
  mats
}

##### makeEqs (DEPRECATED): Generates a set of linear equations based on lagged matrices
makeEqs <- function(dat = NULL, varMods = NULL, mod = "min", type = "g", m = 1, center = TRUE,
                    scale = FALSE, addLast = TRUE, names = TRUE, datBack = FALSE){
  if(!is.null(varMods)){
    if(attr(varMods, "criterion") != "CV"){mod <- "min"}
    mod <- match.arg(mod, choices = c("min", "1se"))
    mod <- ifelse(mod == "min" | attr(varMods, "method") == "regsubsets", 1, 2)
    if(class(dat) == "list"){
      data <- dat$X
      if(any(grepl(":", colnames(data)))){data <- data[,!grepl(":", colnames(data))]}
    }
    newDat <- datMat(data = data, type = type, lags = 1, m = m, center = center, scale = scale)
    if(names == TRUE){
      colnames(newDat$X) <- colnames(dat$X)
      colnames(newDat$full)[!colnames(newDat$full) %in% colnames(dat$Y)] <- colnames(dat$X)
    }
    Ys <- colnames(newDat$Y)
    Xs <- colnames(newDat$X)
    eqs <- list()
    for(i in 1:length(Ys)){
      mains <- varMods[[i]][[mod]][!grepl(":", varMods[[i]][[mod]])]
      ints <- varMods[[i]][[mod]][grepl(":", varMods[[i]][[mod]])]
      if(m == 1){ints <- ints[grepl("X1", ints)]}
      predictors <- c(mains, ints)
      if(length(predictors) == 0){predictors <- 1}
      eqs[[i]] <- as.formula(paste0(Ys[i], " ~ ", paste(predictors, collapse = " + ")))
    }
    if(addLast == TRUE){
      for(j in 1:length(dat)){
        newDat[[j]] <- rbind(newDat[[j]], rep(NA, ncol(newDat[[j]])))
        newDat[[j]][nrow(newDat[[j]]),] <- dat[[j]][nrow(dat[[j]]),]
        rownames(newDat[[j]]) <- NULL
      }
    }
    newDat2 <- list(Y = newDat$Y, X = newDat$X, full = newDat$full, eqs = eqs)
    if(datBack == TRUE){return(newDat2)} else {return(newDat2$eqs)}
  } else {
    if(!is.null(dat)){
      x <- dat$X
      y <- dat$Y
    }
    eqs <- list()
    for(i in 1:ncol(y)){
      eqs[[i]] <- as.formula(paste(colnames(y)[i], "~", paste(colnames(x), collapse = " + ")))
    }
    return(eqs)
  }
}

##### SURll_OLD: log-likelihood of SUR model with LRT to compare models ***FITOBJ***
SURll_OLD <- function(fit, s = "res"){
  if("SURfit" %in% names(fit)){fit <- fit$SURfit}
  s <- match.arg(tolower(s), choices = c("res", "dfres", "sigma"))
  resid <- residuals(fit)
  residCov <- getCoefs(fit = fit, mat = s)
  residCovInv <- solve(residCov)
  resid <- as.matrix(resid)
  nEq <- ncol(resid)
  if((nrow(resid) * ncol(resid)) <= 40000){
    ll <- sum((-nEq/2) * log(2 * pi) - .5 * log(det(residCov)) - .5 * diag(resid %*% residCovInv %*% t(resid)))
  } else {
    ll <- 0
    for(i in 1:nrow(resid)){
      ll <- ll - (nEq/2) * log(2 * pi) - .5 * log(det(residCov)) - .5 * resid[i, , drop = FALSE] %*% residCovInv %*% t(resid[i, , drop = FALSE])
    }
  }
  df <- fit$rank + (nEq * (nEq + 1))/2
  out <- c(LL = as.numeric(ll), df = df)
  out
}


##### SURnet_OLD (DEPRECATED): create temporal and contemporaneous network of SUR results
SURnet_OLD <- function(fit, dat, s = "sigma", m = NULL, type = "g"){
  s <- match.arg(tolower(s), choices = c("sigma", "res", "dfres"))
  if(class(dat) == "list"){data <- dat$X} else {data <- dat} #
  if(any(grepl(":", colnames(data)))){data <- data[,!grepl(":", colnames(data))]} #
  if(length(type) == 1 & ncol(data) > 1){type <- rep(type, ncol(data))} #
  if(is.null(m)){m <- ifelse(any(grepl(":", names(dat$full))), 1, 0)} #
  call <- list(type = type, moderators = m, lags = 1, residMatType = s) #
  beta <- sparsify(getCoefs(fit = fit, mat = "beta", data = dat))
  kappa <- sparsify(solve(getCoefs(fit = fit, mat = s)))
  mainEffects <- beta[,2:(ncol(kappa) + 1)] #
  interactions <- if(m == 1){beta[,(ncol(kappa) + 2):ncol(beta)]} else {NA} #
  PCC <- sparsify(getCoefs(fit = fit, mat = "pcor"))
  PDC <- sparsify(mainEffects/(sqrt(diag(solve(kappa)) %o% diag(kappa) + mainEffects^2))) #
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    if(!any(grepl("lag", colnames(adjMat)))){
      colnames(colMat) <- paste0(colnames(adjMat), ".lag1.")
    } else {
      colnames(colMat) <- colnames(adjMat)
    }
    rownames(colMat) <- rownames(adjMat)
    colMat
  }
  temporal <- list(adjMat = PDC, edgeColors = getEdgeColors(PDC), beta = beta)
  contemporaneous <- list(adjMat = PCC, edgeColors = getEdgeColors(PCC), kappa = kappa)
  colnames(temporal$adjMat) <- colnames(temporal$edgeColors)
  colnames(contemporaneous$adjMat) <- colnames(contemporaneous$edgeColors)
  mods <- lapply(1:nrow(beta), function(z) list(model = as.matrix(beta[z,], ncol = 1)))
  surNet <- list(call = call, temporal = temporal, contemporaneous = contemporaneous,
                 interactions = interactions, mods = mods, data = data)
  if(m == 0){surNet$interactions <- NULL}
  surNet
}


##### lmerVAR_OLD: fit lmerVAR models
lmerVAR_OLD <- function(data = NULL, k = NULL, temp = "correlated", contemp = "correlated", 
                    warnings = FALSE, int = NULL, X, Y, ids = NULL, fixed = NULL, 
                    rand = NULL, formulas = NULL){
  suppressMessages(invisible(c(require(lme4), require(lmerTest))))
  if(!warnings){oldw <- getOption("warn"); options(warn = -1)}
  if(class(data) == "list"){
    dat <- setupVAR(data = data, method = "lmer")
    if("vars" %in% names(data)){k <- length(data$vars)}
  } else {
    dat <- data
  }
  if(missing(X) | missing(Y)){
    if(is.null(dat)){stop("Need to supply data")}
    if(missing(Y)){
      if(is.null(k)){stop("Need to indicate the number of response variables 'k'")}
      Y <- dat[, 1:k]
    } else {
      k <- ncol(Y)
    }
    if(missing(X)){X <- dat[, (k + 1):((2 * k) + k)]}
  }
  if(is.null(k)){k <- ncol(Y)}
  if(is.null(ids)){
    if(is.null(dat)){
      stop("Must provide subject IDs")
    } else if(!any(grepl("ID", colnames(dat)))){
      stop("Must provide subject IDs")
    } else {
      ids <- dat[, grep("ID", colnames(dat))]
    }
  } else if(class(ids) == "character"){
    ids <- dat[, ids]
  }
  if(is.null(dat)){dat <- data.frame(Y, X, ID = ids)} else {dat <- as.data.frame(dat)}
  if(!is.null(int)){stopifnot(all(strsplit(int, ":")[[1]] %in% colnames(dat)))}
  if(!is.null(fixed)){
    if(class(fixed) != "list"){
      fixed <- NULL
    } else if(length(fixed) == 1){
      fixed <- lapply(1:k, function(z) fixed[[1]])
    }
  }
  if(!is.null(rand)){
    if(class(rand) != "list"){
      rand <- NULL
    } else if(length(rand) == 1){
      rand <- lapply(1:k, function(z) rand[[1]])
    }
  }
  tempNet <- match.arg(temp, c("correlated", "orthogonal", "fixed"))
  tempMods <- list()
  for(i in 1:k){
    if(tempNet == "orthogonal"){
      re1 <- paste0(paste0("(0 + ", colnames(X)[1:k], " | ID)"), collapse = " + ")
      if(!is.null(int)){re1 <- paste0(re1, " + (0 + ", int, " | ID)")}
      if(!is.null(rand)){re1 <- paste0(paste0("(0 + ", rand[[i]], " | ID)"), collapse = " + ")}
      re <- paste0("((1 | ID) + ", re1, ")")
    } else if(tempNet == "correlated"){
      re1 <- ifelse(!is.null(int), paste0(" + ", int), "")
      re <- paste0("(", paste0(colnames(X)[1:k], collapse = ' + '), re1, " | ID)")
      if(!is.null(rand)){re <- paste0("(", paste0(rand[[i]], collapse = " + "), " | ID)")}
    } else if(tempNet == "fixed"){
      re <- "(1 | ID)"
    }
    fe <- paste0(colnames(X)[-(i + k)], collapse = " + ")
    if(!is.null(int)){fe <- paste0(fe, " + ", int)}
    if(!is.null(fixed)){fe <- paste0(fixed[[i]], collapse = " + ")}
    tempMods[[i]] <- as.formula(paste0(colnames(Y)[i], " ~ ", fe, " + ", re))
  }
  if(!is.null(formulas)){if(class(formulas) == "list"){tempMods <- lapply(formulas, as.formula)}}
  mods <- lapply(tempMods, function(z) lmer(z, data = dat, REML = F))
  names(mods) <- colnames(Y)
  resDat <- do.call(cbind.data.frame, lapply(mods, resid))
  colnames(resDat) <- colnames(Y)
  contempNet <- match.arg(contemp, c("correlated", "orthogonal"))
  resMods <- list()
  for(i in 1:k){
    if(contempNet == "orthogonal"){
      re1 <- paste0(paste0("(0 + ", colnames(resDat)[-i], " | ID)"), collapse = " + ")
      re <- paste0("(", re1, ")")
    } else if(contempNet == "correlated"){
      re <- paste0("(0 + ", paste0(colnames(resDat)[-i], collapse = " + "), " | ID)")
    }
    fe <- paste0(colnames(resDat)[-i], collapse = " + ")
    resMods[[i]] <- as.formula(paste0(colnames(resDat)[i], " ~ 0 + ", fe, " + ", re))
  }
  resDat$ID <- ids
  contempMods <- lapply(resMods, function(z) lmer(z, data = resDat, REML = F))
  names(contempMods) <- colnames(Y)
  fit <- do.call(rbind, lapply(mods, function(z) c(AIC(z), BIC(z))))
  fit <- cbind.data.frame(var = colnames(Y), fit)
  colnames(fit)[2:3] <- c("aic", "bic")
  out <- list(tempMods = mods, contempMods = contempMods, fit = fit)
  attributes(out)$temporal <- paste0(tempNet, ifelse(!is.null(int), " (interaction)", ""))
  attributes(out)$contemporaneous <- contempNet
  if(!warnings){options(warn = oldw)}
  out
}



##### lmerNets_OLD (although some new edits were made): create networks from lmerVAR models
lmerNets_OLD <- function(model, data, threshold = FALSE, rule = "OR", 
                     ggm = "pcor", ggmAvg = TRUE){
  if(is.null(data) & "data" %in% names(model)){
    data <- model$data
  } else if(is(data, "list")){
    dat0 <- data
    idvar <- ifelse("idvar" %in% names(dat0), dat0[["idvar"]], "ID")
    data <- setupVAR(data = data, idvar = idvar, method = "lmer")
  }
  forcePositive <- function(x, avg = TRUE){
    if(avg){x <- (x + t(x))/2} else {x <- sign(x) * sqrt(x * t(x))}
    if(any(eigen(x)$values < 0)){x <- x - diag(nrow(x)) * min(eigen(x)$values) - 0.001}
    return(x)
  }
  if("inds" %in% names(model)){
    y <- model$inds$yvars
    x <- model$inds$mvars
    k <- length(y)
  } else {
    y <- model$fit$var
    k <- length(y)
    x <- colnames(data)[(2 * k + 1):(2 * k + k)]
  }
  y2 <- rep(list(y), 2)
  beta <- do.call(rbind, lapply(model$tempMods, function(z) fixef(z)[2:(k + 1)]))
  betaSE <- do.call(rbind, lapply(model$tempMods, function(z) arm::se.fixef(z)[2:(k + 1)]))
  beta_pvals <- (1 - pnorm(abs(beta/betaSE))) * 2
  dimnames(beta) <- dimnames(betaSE) <- dimnames(beta_pvals) <- y2
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
  dimnames(gammaTheta) <- dimnames(gammaThetaSE) <- dimnames(gammaTheta_pvals) <- y2
  D <- diag(1/sapply(model$contempMods, sigma)^2)
  inv <- tryCatch({forcePositive(D %*% (diag(k) - gammaTheta), avg = ggmAvg)}, error = function(e){
    warning("Non-convergent estimate of Theta"); return(diag(k))})
  Theta_prec <- inv
  Theta_cov <- corpcor::pseudoinverse(inv)
  Theta_cor <- cov2cor(Theta_cov)
  d <- 1/sqrt(diag(inv))
  Theta <- -t(d * inv) * d
  diag(Theta) <- 0
  PDC <- beta/(sqrt(diag(Theta_cov) %o% diag(Theta_prec) + beta^2))
  dimnames(Theta) <- dimnames(Theta_cov) <- dimnames(Theta_cor) <- dimnames(Theta_prec) <- y2
  gammaOmega <- lapply(model$tempMods, function(z){
    z2 <- fixef(z)
    return(z2[which(names(z2) %in% x)])
  })
  gammaOmegaSE <- lapply(model$tempMods, function(z){
    z2 <- arm::se.fixef(z)
    return(z2[which(names(z2) %in% x)])
  })
  go1 <- go2 <- matrix(NA, k, k)
  for(i in 1:k){
    go1[i, match(names(gammaOmega[[i]]), x)] <- gammaOmega[[i]]
    go2[i, match(names(gammaOmegaSE[[i]]), x)] <- gammaOmegaSE[[i]]
  }
  diag(go1) <- diag(go2) <- 0
  gammaOmega <- go1
  gammaOmegaSE <- go2
  gammaOmega_pvals <- (1 - pnorm(abs(gammaOmega/gammaOmegaSE))) * 2
  dimnames(gammaOmega) <- dimnames(gammaOmegaSE) <- dimnames(gammaOmega_pvals) <- y2
  mu_SD <- sapply(model$tempMods, function(z) attr(lme4::VarCorr(z)[[1]], "stddev")[1])
  D <- diag(1/mu_SD^2)
  inv <- tryCatch({forcePositive(D %*% (diag(k) - gammaOmega), avg = ggmAvg)}, error = function(e){
    warning("Non-convergent estimate of Omega"); return(diag(k))})
  Omega_prec <- inv
  Omega_cov <- corpcor::pseudoinverse(inv)
  Omega_cor <- cov2cor(Omega_cov)
  d <- 1/sqrt(diag(inv))
  Omega <- -t(d * inv) * d
  diag(Omega) <- 0
  dimnames(Omega) <- dimnames(Omega_cov) <- dimnames(Omega_cor) <- dimnames(Omega_prec) <- y2
  out <- list(Beta = list(mu = beta, SE = betaSE, Pvals = beta_pvals), 
              Theta = list(cor = Theta_cor, cov = Theta_cov, pcor = Theta, prec = Theta_prec), 
              Omega = list(cor = Omega_cor, cov = Omega_cov, pcor = Omega, prec = Omega_prec), 
              gammaTheta = list(mu = gammaTheta, SE = gammaThetaSE, Pvals = gammaTheta_pvals),
              gammaOmega = list(mu = gammaOmega, SE = gammaOmegaSE, Pvals = gammaOmega_pvals))
  out$model <- model
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    if(!any(grepl("lag", colnames(adjMat)))){
      colnames(colMat) <- paste0(colnames(adjMat), ".lag1.")
    } else {
      colnames(colMat) <- colnames(adjMat)
    }
    rownames(colMat) <- rownames(adjMat)
    colMat
  }
  if(threshold != FALSE){
    if(threshold == TRUE){threshold <- 0.05}
    thresh <- function(obj, alpha = 0.05, rule = "OR", ggm = "pcor"){
      rule <- tryCatch({match.arg(tolower(rule), c("or", "and"))}, 
                       error = function(e) return("or"))
      ggm <- tryCatch({match.arg(tolower(ggm), c("pcor", "cor", "cov", "prec"))}, 
                      error = function(e) return("pcor"))
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
  colnames(beta) <- colnames(PDC) <- paste0(colnames(beta), ".lag1.")
  out2 <- list(temporal = list(adjMat = beta, edgeColors = getEdgeColors(beta)),
               contemporaneous = list(adjMat = Theta, edgeColors = getEdgeColors(Theta)),
               between = list(adjMat = Omega, edgeColors = getEdgeColors(Omega)))
  out2$temporal <- append(out2$temporal, list(
    PDC = list(adjMat = PDC, edgeColors = getEdgeColors(PDC))))
  out2$contemporaneous$kappa <- Theta_prec
  out2$modelCoefs <- out
  out <- out2
  if(exists("dat0")){out$data <- dat0$data[, 1:k]} else {out$data <- data[, 1:k]}
  if(exists("idvar")){out$ID <- dat0$data[, idvar]} else {out$ID <- data[, "ID"]}
  colnames(out$contemporaneous$edgeColors) <- colnames(out$between$edgeColors) <- y
  attributes(out)$lmer <- TRUE
  out
}


################################################################################
################################ CAN REJUVENATE ################################
################################################################################
##### compareMats: compute estimated and true parameter matrices
compareMats <- function(B, S, fits, d = 5, ind = "cosine", similarity = TRUE,
                        labs = c("fitOLS", "fit0", "fit1se"), dat = NULL){
  if(all(sapply(fits, function(z) "SURfit" %in% names(z))) & 1 == 2){
    dats <- ifelse(length(unique(lapply(fits, '[[', "dat"))) == 1, TRUE, FALSE)
    mm <- length(unique(sapply(lapply(fits, '[[', "call"), '[[', "moderators"))) == 1
    cc <- length(unique(sapply(lapply(fits, '[[', "call"), '[[', "covariates"))) == 1
  } # unfinished...
  if(!is.null(dat)){
    if(is.logical(dat)){dat <- fits[[1]]$dat}
    if(fits[[1]]$call$exogenous & !is.null(fits[[1]]$call$moderators)){
      m <- fits[[1]]$call$moderators
      B <- B[-m, ]
      S <- S[-m, -m]
    }
  }
  if(all(sapply(fits, function(z) "SURfit" %in% names(z)))){
    fits <- lapply(fits, function(z) z$SURfit)}
  x <- matrix(NA, ncol = 5, nrow = 3); s <- similarity
  x[,1] <- sapply(fits, function(z){matrixDist(B, getCoefs(z, "beta", data = dat), ind, s)})
  x[,2] <- sapply(fits, function(z){matrixDist(S, getCoefs(z, "sigma"), ind, s)})
  x[,3] <- sapply(fits, function(z){matrixDist(S, getCoefs(z, "res"), ind, s)})
  x[,4] <- sapply(fits, function(z){matrixDist(S, getCoefs(z, "cor"), ind, s)})
  k <- tryCatch({-solve(S)}, error = function(e){-corpcor::pseudoinverse(S)})
  diag(k) <- -diag(k)
  K <- t((1/sqrt(diag(k))) * k) * (1/sqrt(diag(k)))
  x[,5] <- sapply(fits, function(z){matrixDist(K, getCoefs(z, "pcor"), ind, s)})
  if(is.null(names(fits))){
    if(all(labs == c("fitOLS", "fit0", "fit1se"))){
      rownames(x) <- labs[rev(order(sapply(fits, function(z) z$rank)))]
    } else {rownames(x) <- labs}
  }
  colnames(x) <- c("Beta", "Sigma", "ResidCov", "ResidCor", "ResidPcor")
  round(x, digits = d)
}

##### results: print or save SUR model output (true+estimated matrices for beta/sigma/cor)
results <- function(true, est, x = "beta", d = NULL, labels = NULL,
                    intercepts = TRUE, keep = FALSE){
  x <- match.arg(x, choices = c("beta", "sigma", "cor", "res", "dfres", "pcor"))
  if(dim(true)[1] == dim(true)[2] & all(diag(true) == 1) & x == "beta"){x <- "sigma"}
  if(dim(true)[1] != dim(true)[2] & x != "beta"){x <- "beta"}
  x2 <- ifelse(x == "beta", "BETA", ifelse(x %in% c("sigma", "res", "dfres"), "RESID.COV", 
                                           ifelse(x == "cor", "RESID.COR", "RESID.PCOR")))
  if(x == "pcor"){
    true <- tryCatch({-solve(true)}, error = function(e){
      -corpcor::pseudoinverse(true)})
    diag(true) <- -diag(true)
    delta <- (1/sqrt(diag(true)))
    true <- t(delta * true) * delta
    if(is.null(d)){d <- 2}
  }
  if(is.null(d)){d <- ifelse(x == "beta", 2, 1)}
  int <- ifelse(x == "beta", intercepts, TRUE)
  if(is.null(labels)){
    ilabs <- rownames(getCoefs(fit = est, mat = x))
    jlabs <- colnames(getCoefs(fit = est, mat = x))
  } else if(class(labels) == "systemfit"){
    ilabs <- rownames(getCoefs(labels, mat = x))
    jlabs <- colnames(getCoefs(labels, mat = x))
  } else {ilabs <- labels$ilabs; jlabs <- labels$jlabs}
  if(keep == FALSE){
    if(int == TRUE){
      cat(paste0("\nTRUE ", x2, "\n"))
      print(sparsify(mat = true, d = d, ilabs = ilabs, jlabs = jlabs))
      cat(paste0("\nESTIMATED ", x2, "\n"))
      print(sparsify(mat = getCoefs(fit = est, mat = x), d = d, 
                     ilabs = ilabs, jlabs = jlabs))
    } else {
      cat(paste0("\nTRUE ", x2, "\n"))
      print(sparsify(mat = true, d = d, ilabs = ilabs, jlabs = jlabs)[,-1])
      cat(paste0("\nESTIMATED ", x2, "\n"))
      print(sparsify(mat = getCoefs(fit = est, mat = x), d = d,
                     ilabs = ilabs, jlabs = jlabs)[,-1])
    }
  } else {
    m1 <- sparsify(mat = true, d = d, ilabs = ilabs, jlabs = jlabs)
    m2 <- sparsify(mat = getCoefs(fit = est, mat = x), d = d, 
                   ilabs = ilabs, jlabs = jlabs)
    if(int == FALSE){m1 <- m1[,-1]; m2 <- m2[,-1]}
    m <- list(m1, m2)
    names(m) <- c(paste0("true", Hmisc::capitalize(x)), paste0("est", Hmisc::capitalize(x)))
    m
  }
}


################################################################################
######################### CREATED FOR CONVENIENCE ONLY #########################
################################################################################
##### makeInts: create all interaction terms from a design matrix
makeInts <- function(x, mat = TRUE){
  k <- ncol(x)
  stopifnot(k > 1)
  k2 <- (k * (k - 1))/2
  x2 <- data.frame(matrix(NA, ncol = k2, nrow = nrow(x)))
  ints <- t(combn(x = c(1:k), m = 2))
  for(i in 1:k2){x2[, i] <- x[, ints[i, 1]] * x[, ints[i, 2]]}
  colnames(x2) <- paste0(colnames(x)[ints[,1]], ":", colnames(x)[ints[,2]])
  if(mat){x2 <- as.matrix(x2)}
  x2
}

##### fdr: False detection rates
fdr <- function(res, B, all = FALSE, ss = FALSE){
  if(ncol(B) == (2 * nrow(B))){B <- B[,-1]}
  B[B != 0] <- 1
  B <- t(B)
  m1 <- res$mod0$freqs * length(res$samps)
  m1 <- split(m1, rep(1:ncol(m1), each = nrow(m1)))
  for(i in 1:length(m1)){m1[[i]] <- cbind(m1[[i]], B[,i])}
  possiblePos <- as.vector(sapply(m1, function(z) sum(z[,2]))) * length(res$samps)
  possibleNeg <- as.vector(sapply(m1, function(z) sum(z[,2] == 0))) * length(res$samps)
  tp1 <- as.vector(sapply(m1, function(z) sum(z[,1] * z[,2])))
  fn1 <- possiblePos - tp1
  fp1 <- as.vector(sapply(m1, function(z) sum(z[z[,2] == 0, 1])))
  tn1 <- possibleNeg - fp1
  acc1 <- (tp1 + tn1)/(tp1 + fp1 + tn1 + fn1)
  prc1 <- tp1/(tp1 + fp1)
  fdr1 <- fp1/(tp1 + fp1)
  sens1 <- tp1/(tp1 + fn1)
  spec1 <- tn1/(tn1 + fp1)
  m2 <- res$mod1se$freqs * length(res$samps)
  m2 <- split(m2, rep(1:ncol(m2), each = nrow(m2)))
  for(i in 1:length(m2)){m2[[i]] <- cbind(m2[[i]], B[,i])}
  possiblePos <- as.vector(sapply(m2, function(z) sum(z[,2]))) * length(res$samps)
  possibleNeg <- as.vector(sapply(m2, function(z) sum(z[,2] == 0))) * length(res$samps)
  tp2 <- as.vector(sapply(m2, function(z) sum(z[,1] * z[,2])))
  fn2 <- possiblePos - tp2
  fp2 <- as.vector(sapply(m2, function(z) sum(z[z[,2] == 0, 1])))
  tn2 <- possibleNeg - fp2
  acc2 <- (tp2 + tn2)/(tp2 + fp2 + tn2 + fn2)
  prc2 <- tp2/(tp2 + fp2)
  fdr2 <- fp2/(tp2 + fp2)
  sens2 <- tp2/(tp2 + fn2)
  spec2 <- tn2/(tn2 + fp2)
  if(ss == TRUE){
    sss <- round(data.frame(Sens0 = sens1, Spec0 = spec1, Sens1se = sens2, Spec1se = spec2), 2)
    return(sss)
  }
  if(all == FALSE){
    round(data.frame(ACC0 = acc1, PRC0 = prc1, FDR0 = fdr1, Sens0 = sens1, 
                     Spec0 = spec1, ACC1se = acc2, PRC1se = prc2, FDR1se = fdr2, 
                     Sens1se = sens2, Spec1se = spec2), 2)
  } else {
    list(mod0 = data.frame(tp1, fp1, tn1, fn1, acc1, prc1, fdr1, sens1, spec1), 
         mod1se = data.frame(tp2, fp2, tn2, fn2, acc2, prc2, fdr2, sens2, spec2))
  }
}
makeFreqs <- function(res){
  niter <- length(res)
  allNames <- paste0("X", 1:4)
  allNames <- c(allNames, apply(combn(allNames, 2), 2, paste, collapse = ":"))[1:7]
  x <- list()
  for(i in 1:4){
    x[[i]] <- lapply(res, function(z) as.character(z[,i]))
    x[[i]] <- unname(unlist(x[[i]]))
    x[[i]] <- count(x[[i]][x[[i]] != "" & x[[i]] != "X2:X3" & x[[i]] != "X2:X4" & x[[i]] != "X3:X4"])
    x[[i]] <- x[[i]][match(allNames, x[[i]][,1]),]
  }
  m1 <- do.call(cbind, x)[,seq(2, 8, by = 2)]
  colnames(m1) <- paste0("X", 1:4, ".y")
  rownames(m1) <- allNames
  m1 <- as.matrix(m1)
  if(any(is.na(m1))){m1[is.na(m1)] <- 0}
  m1
}
getResids <- function(obj){
  r0 <- t(sapply(obj$fits, function(z) sapply(z$mod0$eq, function(zz) zz$df.residual)))
  r1 <- t(sapply(obj$fits, function(z) sapply(z$mod1se$eq, function(zz) zz$df.residual)))
  res0 <- lapply(obj$fits, function(z) residuals(z$mod0))
  res1 <- lapply(obj$fits, function(z) residuals(z$mod1se))
  colnames(r0) <- colnames(r1) <- colnames(res0[[1]])
  list(resid0 = res0, dfres0 = r0, resid1 = res1, dfres1 = r1)
}
mse <- function(x){
  rss0 <- lapply(x$residuals$resid0, function(z){
    apply(z, 2, function(zz) sum(zz^2))
  })
  rss0 <- do.call(rbind, rss0)
  mse0 <- matrix(0, ncol = ncol(x$residuals$dfres0), nrow = nrow(x$residuals$dfres0))
  for(i in 1:nrow(x$residuals$dfres0)){
    mse0[i,] <- rss0[i,]/x$residuals$dfres0[i,]
  }
  rss1 <- lapply(x$residuals$resid1, function(z){
    apply(z, 2, function(zz) sum(zz^2))
  })
  rss1 <- do.call(rbind, rss1)
  mse1 <- matrix(0, ncol = ncol(x$residuals$dfres1), nrow = nrow(x$residuals$dfres1))
  for(i in 1:nrow(x$residuals$dfres1)){
    mse1[i,] <- rss1[i,]/x$residuals$dfres1[i,]
  }
  list(mse0 = mse0, mse1 = mse1)
}
fdr_penalty <- function(p, alpha = 0.05, niter = 100){
  pmax(0.05, ((1 - log(ceiling(alpha * niter)/niter))/alpha) * p)
}

##### sparsify: convert any (i x j) matrix into a sparse matrix
# was removed from SURnet and getCoefs
sparsify <- function(mat, d = NULL, ilabs = NULL, jlabs = NULL){
  stopifnot(length(dim(mat)) == 2)
  p <- ncol(mat)
  n <- nrow(mat)
  if(!is.null(rownames(mat)) & is.null(ilabs)){ilabs <- rownames(mat)}
  if(!is.null(colnames(mat)) & is.null(jlabs)){jlabs <- colnames(mat)}
  if(class(mat) == "data.frame"){mat <- as.matrix(mat, ncol = p, nrow = n)}
  if(!is.null(d)){mat <- round(mat, digits = d)}
  dat <- data.frame(i = rep(1:n, p), j = rep(1:p, each = n), x = as.vector(mat))
  dat <- dat[dat$x != 0 & !is.na(dat$x), ]
  dat2 <- Matrix::sparseMatrix(i = dat$i, j = dat$j, x = dat$x, dims = c(n, p))
  if(!is.null(ilabs)){rownames(dat2) <- ilabs}
  if(!is.null(jlabs)){colnames(dat2) <- jlabs}
  dat2
}
