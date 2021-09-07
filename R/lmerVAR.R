#' Legit mixed-effects modeling for the GVAR in multilevel data
#'
#' Something or other
#'
#' @param data data.frame
#' @param m numeric
#' @param temporal character
#' @param contemp character
#' @param idvar character
#' @param intvars character vector
#' @param center logical
#' @param scale logical
#' @param centerWithin logical
#' @param scaleWithin logical
#' @param exogenous logical
#' @param covariates character or list
#' @param fix something
#' @param warnings logical
#' @param verbose logical
#' @param beepno something
#' @param dayno something
#' @param deleteMissing logical
#'
#' @return lmerVAR models
#' @export
#'
#' @examples
#' 1 + 1
lmerVAR <- function(data, m = NULL, temporal = "default", contemp = "default",
                    idvar = "ID", intvars = NULL, center = TRUE, scale = TRUE,
                    centerWithin = TRUE, scaleWithin = FALSE, exogenous = TRUE,
                    covariates = NULL, fix = NULL, warnings = FALSE, verbose = TRUE,
                    beepno = NULL, dayno = NULL, deleteMissing = TRUE){
  t1 <- Sys.time()
  #suppressMessages(invisible(c(require(lme4), require(lmerTest))))
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
    # CHECK THAT THIS SHOULD INDEED BE lmerTest!
    tm <- suppressMessages(lmerTest::lmer(tempForms[[z]], data = dat1, REML = FALSE))
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
    # CHECK THAT THIS SHOULD INDEED BE lmerTest!
    cm <- suppressMessages(lmerTest::lmer(contempForms[[z]], data = resDat, REML = FALSE))
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

#' Create networks out of lmerVAR models
#'
#' Gotta look into how this works
#'
#' @param model output from lmerVAR?
#' @param inds something
#' @param m mnames
#' @param threshold logical numeric
#' @param rule character
#' @param ggm character
#'
#' @return lmerNetworks
#' @export
#'
#' @examples
#' 1 + 1
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
    model$tempMods, function(z) lme4::fixef(z)[x]))
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
  gammaTheta <- lapply(model$contempMods, lme4::fixef)
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
    z2 <- lme4::fixef(z)
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
      lme4::fixef(z)[grepl(":", names(lme4::fixef(z)))]}))
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


#' Compare two to three lmerVAR models
#'
#' Gotta see how this works
#'
#' @param m1 model 1
#' @param m2 model 2
#' @param m3 model 3 (optional)
#' @param anova numeric?
#' @param type character
#'
#' @return Table of ANOVA
#' @export
#'
#' @examples
#' 1 + 1
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
