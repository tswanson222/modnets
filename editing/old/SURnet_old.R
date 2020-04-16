SURnet <- function(fit, dat, s = "sigma", m = NULL, threshold = FALSE, 
                   mval = NULL, medges = 1){
  s <- match.arg(tolower(s), choices = c("sigma", "res", "dfres"))
  beta <- getCoefs(fit = fit, mat = "beta", data = dat)
  kappa <- solve(getCoefs(fit = fit, mat = s))
  ynames <- gsub("[.]y$", "", rownames(beta))
  mainEffects <- beta[, ynames, drop = FALSE]
  intnames <- grepl(":", colnames(beta))
  interactions <- beta[, intnames, drop = FALSE]
  if(ncol(interactions) == 0){interactions <- NULL}
  PCC <- getCoefs(fit = fit, mat = "pcor")
  PDC <- mainEffects/(sqrt(diag(solve(kappa)) %o% diag(kappa) + mainEffects^2))
  dimnames(PCC) <- dimnames(kappa) <- dimnames(PDC)
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
    return(colMat)
  }
  fitobj <- lapply(fit$eq, function(z){
    attributes(z)$family <- "gaussian"
    return(z)
  })
  attr(fitobj, "rank") <- fit$rank
  pvals <- getCoefs(fit = fit, mat = "pvals", data = dat)
  if(!is.null(interactions)){
    stopifnot(!is.null(m))
    interactions <- list(beta = interactions, pvals = pvals[, intnames, drop = FALSE])
    modEdges <- ifelse(interactions$pvals <= ifelse(is.logical(threshold), .05, threshold), 1, 0)
    if(any(interactions$pvals == 0)){modEdges <- modEdges * ifelse(interactions$pvals == 0, 0, 1)}
    modEdges <- modEdges + 1
    if("covariates" %in% names(attributes(fit)) & length(m) == 1){
      m <- m - sum(attr(fit, "covariates") < m)
    }
    if(length(m) == 1 & !attr(fit, "exogenous")){
      modEdges0 <- matrix(1, length(ynames), length(ynames))
      modEdges0[, -m] <- modEdges
      modEdges0[, m] <- unname(apply(modEdges, 1, function(z) ifelse(any(z == 2), medges, 1)))
      rownames(modEdges0) <- rownames(modEdges)
      mcols0 <- character(length(ynames))
      mcols0[m] <- paste(rep(attr(fit, "moderators"), 2), collapse = ":")
      mcols0[-m] <- colnames(modEdges)
      colnames(modEdges0) <- mcols0
      modEdges <- modEdges0
    }
  } else {
    modEdges <- matrix(1, length(ynames), length(ynames))
  }
  if(threshold != FALSE){
    if(threshold == TRUE){threshold <- .05}
    mainEffects <- mainEffects * ifelse(pvals[, ynames, drop = FALSE] <= threshold, 1, 0)
    PDC <- PDC * ifelse(pvals[, ynames, drop = FALSE] <= threshold, 1, 0)
    if(!is.null(interactions)){
      interactions$beta <- interactions$beta * ifelse(interactions$pvals <= threshold, 1, 0)
    }
  }
  temporal <- list(adjMat = mainEffects, edgeColors = getEdgeColors(mainEffects), 
                   modEdges = modEdges, PDC = list(adjMat = PDC, edgeColors = getEdgeColors(PDC)),
                   coefs = list(beta = beta, pvals = pvals))
  if(is.null(interactions)){temporal$modEdges <- NULL}
  colnames(temporal$adjMat) <- colnames(temporal$PDC$adjMat) <- colnames(temporal$edgeColors)
  contemporaneous <- list(adjMat = PCC, edgeColors = getEdgeColors(PCC), kappa = kappa)
  colnames(contemporaneous$adjMat) <- colnames(contemporaneous$edgeColors)
  y <- dat$Y
  n <- nrow(y)
  yhat <- predict(fit)
  mods <- lapply(1:nrow(beta), function(z){
    model <- as.matrix(beta[z, ], ncol = 1)
    deviance <- sum((y[, z] - yhat[, z])^2)
    s <- sqrt(deviance/n)
    LL_model <- sum(dnorm(y[, z], mean = yhat[, z], sd = s, log = TRUE))
    k <- nrow(model) + 1
    aic <- (2 * k) - (2 * LL_model)
    bic <- (log(n) * k) - (2 * LL_model)
    out <- list(deviance = deviance, LL_model = LL_model, AIC = aic, BIC = bic, model = model)
    return(out)
  })
  mods0 <- lapply(mods, '[[', "model")
  names(mods) <- names(mods0) <- names(fitobj) <- sapply(fitobj, function(z) as.character(z$terms[[2]]))
  call <- list(type = rep("g", ncol(dat$Y)), moderators = m, mval = mval,
               lags = 1, residMat = s, threshold = threshold)
  surNet <- list(call = call, temporal = temporal, contemporaneous = contemporaneous)
  if(!is.null(interactions)){
    surNet$interactions <- interactions
    surNet$call$moderators <- mname <- attr(fit, "moderators")
    surNet$call$exogenous <- attr(fitobj, "exogenous") <- attr(fit, "exogenous")
    if(length(m) == 1){
      ints <- lapply(mods0, function(z) rownames(z)[z[, 1] != 0][grep(":", rownames(z)[z[, 1] != 0])])
      inds0 <- unlist(lapply(ints, function(z) gsub(paste0(mname, ":|:", mname), "", z)))
      inds1 <- data.frame(y = rep(names(mods), sapply(ints, length)), ints = unname(unlist(ints)))
      vars <- lapply(fitobj, vcov)
      surNet$interactions$coefvars <- vars1 <- data.frame(
        Y = inds1[, 1], X = unname(inds0), Z = mname, Int = inds1[, 2], 
        t(sapply(seq_len(nrow(inds1)), function(i){
          vb1 <- vars[[inds1[i, 1]]][inds0[i], inds0[i]]
          vb3 <- vars[[inds1[i, 1]]][inds1[i, 2], inds1[i, 2]]
          vb1b3 <- vars[[inds1[i, 1]]][inds0[i], inds1[i, 2]]
          return(c(varX = vb1, varInt = vb3, varCov = vb1b3))
        }))
      )
      if(attr(fit, "exogenous")){
        if(!is.null(mval)){
          margSE <- function(x, vars){
            if(length(vars) > 3 & !is.null(names(vars))){vars <- vars[grep("var", names(vars))]}
            as.numeric(sqrt(vars[1] + ((x^2) * vars[2]) + (2 * x * vars[3])))
          }
          ses <- getCoefs(fit = fit, mat = "ses", data = dat)
          for(i in 1:nrow(vars1)){
            B3Z <- mval * beta[vars1[i, "Y"], vars1[i, "Int"]]
            beta[vars1[i, "Y"], vars1[i, "X"]] <- beta[vars1[i, "Y"], vars1[i, "X"]] + B3Z
            ses[vars1[i, "Y"], vars1[i, "X"]] <- margSE(mval, vars1[i, ])
          }
          dfs <- matrix(rep(sapply(fitobj, '[[', "df.residual"), each = ncol(beta)), 
                        nrow(beta), ncol(beta), TRUE)
          pvals <- (2 * pt(abs(beta/ses), df = dfs, lower.tail = FALSE))
          if(any(is.na(pvals))){pvals[is.na(pvals)] <- 1}
        }
        mp <- length(ynames)
        madj <- modEdges0 <- matrix(0, mp + 1, mp + 1)
        madj[-m, -m] <- surNet$temporal$adjMat
        dimnames(madj) <- rep(list(1:(mp + 1)), 2)
        rownames(madj)[-m] <- rownames(surNet$temporal$adjMat)
        rownames(madj)[m] <- paste0(mname, ".y")
        colnames(madj)[-m] <- colnames(surNet$temporal$adjMat)
        colnames(madj)[m] <- paste0(mname, ".lag1.")
        if(threshold != FALSE){
          madj[-m, m] <- beta[, mname] * ifelse(pvals[, mname] <= threshold, 1, 0)
        } else {
          madj[-m, m] <- beta[, mname]
        }
        modEdges0[-m, -m] <- modEdges
        modEdges0[-m, m] <- unname(apply(modEdges, 1, function(z){
          ifelse(any(z == 2), medges, 1)}))
        dimnames(modEdges0) <- dimnames(madj)
        shape <- rep("circle", mp + 1)
        shape[m] <- "square"
        surNet$mnet <- list(adjMat = madj, edgeColors = getEdgeColors(madj),
                            modEdges = modEdges0, shape = shape)
      }
    } else if(length(m) > 1){
      surNet$temporal$modEdges <- NULL
    }
  }
  if(is.null(mval) | is.null(m) | length(m) > 1){surNet$call$mval <- NULL}
  surNet <- append(surNet, list(mods = mods, fitobj = fitobj, data = dat))
  names(surNet$fitobj) <- names(mods)
  return(surNet)
}
