#' Fit SUR model with or without constraints
#'
#' Actually not sure if this is a top-level function.
#'
#' @param data dataset
#' @param varMods output of varSelect
#' @param mod character
#' @param maxiter numeric
#' @param m moderator
#' @param type character
#' @param center logical
#' @param scale logical
#' @param exogenous logical
#' @param covs something
#' @param sur logical
#' @param consec something
#' @param ... Additional arguments.
#'
#' @return Models
#' @export
#'
#' @examples
#' 1 + 1
SURfit <- function(data, varMods = NULL, mod = "min", maxiter = 100, m = NULL,
                   type = "g", center = TRUE, scale = FALSE, exogenous = TRUE,
                   covs = NULL, sur = TRUE, consec = NULL, ...){
  if(!is(varMods, 'list')){type <- varMods; varMods <- NULL}
  eqs <- surEqs(data = data, varMods = varMods, mod = match.arg(mod, c("min", "1se")),
                m = m, exogenous = exogenous, covs = covs)
  dat <- lagMat(data = data, type = type, m = m, covariates = covs, center = center,
                scale = scale, exogenous = exogenous, consec = consec)
  fit <- systemfit::systemfit(formula = eqs, method = ifelse(sur, "SUR", "OLS"),
                              data = cbind.data.frame(dat$Y, dat$X),
                              maxiter = maxiter, ...) # FULLFIX
  for(i in 1:length(fit$eq)){attr(fit$eq[[i]], "family") <- "gaussian"}
  names(fit$eq) <- colnames(dat$Y)
  if(!is.null(m)){
    attr(fit, "exogenous") <- exogenous
    attr(fit, "moderators") <- colnames(data)[m]
  }
  if(!is.null(covs)){attr(fit, "covariates") <- covs}
  fit
}

#' Log-likelihood of SUR model with LRT to compare models
#'
#' Similar to modLL
#'
#' @param net0 network
#' @param net1 network
#' @param nodes logical
#' @param lrt logical
#' @param all logical
#' @param d numeric
#' @param alpha numeric
#' @param s character
#' @param orderBy character
#' @param decreasing logical
#' @param sysfits logical
#'
#' @return Results about SUR network
#' @export
#'
#' @examples
#' 1 + 1
SURll <- function(net0, net1 = NULL, nodes = FALSE, lrt = NULL, all = FALSE,
                  d = 4, alpha = .05, s = "res", orderBy = NULL,
                  decreasing = TRUE, sysfits = FALSE){
  sll <- function(fit, s = "res"){
    if("SURfit" %in% names(fit)){fit <- fit$SURfit}
    s <- match.arg(tolower(s), choices = c("res", "dfres", "sigma"))
    resid <- residuals(fit)
    residCov <- getCoefs(fit = fit, mat = s)
    residCovInv <- solve(residCov)
    resid <- as.matrix(resid)
    nEq <- ncol(resid)
    ll <- 0
    for(i in 1:nrow(resid)){
      ll <- ll - (nEq/2) * log(2 * pi) - .5 * log(det(residCov)) - .5 * resid[i, , drop = FALSE] %*% residCovInv %*% t(resid[i, , drop = FALSE])
    }
    df <- fit$rank + (nEq * (nEq + 1))/2
    out <- c(LL = as.numeric(ll), df = df)
    out
  }
  uni_sll <- function(fit){
    if("SURnet" %in% names(fit)){fit <- append(fit$SURnet, list(fitobj = fit$SURfit$eq))}
    getInd <- function(ind, x){
      ind <- c("deviance", "LL_model", "df.residual", "AIC", "BIC")[ind]
      X <- ifelse(ind == "df.residual", list(x$fitobj), list(x$mods))[[1]]
      return(unname(sapply(X, '[[', ind)))
    }
    out <- do.call(cbind.data.frame, lapply(1:5, getInd, x = fit))
    colnames(out) <- c("RSS", "LL", "df", "AIC", "BIC")
    rownames(out) <- names(fit$mods)
    out
  }
  omni_sll <- function(fit, s = "res"){
    k <- length(fit$SURnet$mods)
    n <- nrow(fit$SURnet$data$X) * k
    ll <- sll(fit = fit, s = s)
    aic <- (2 * ll[2]) - (2 * ll[1])
    bic <- (ll[2] * log(n)) - (2 * ll[1])
    out <- c(ll, AIC = unname(aic), BIC = unname(bic))
    out
  }
  sur_lrt <- function(object, d = 4, alpha = .05, N = NULL){
    if(is.list(object)){
      if(length(object) > 2){object <- object[1:2]}
      nn <- names(object)
      ll0 <- object[[1]]$LL; df0 <- object[[1]]$df
      ll1 <- object[[2]]$LL; df1 <- object[[2]]$df
      omnibus <- FALSE
    } else {
      if(nrow(object) > 2){object <- object[1:2, ]}
      nn <- rownames(object)
      ll0 <- object[1, 1]; df0 <- object[1, 2]
      ll1 <- object[2, 1]; df1 <- object[2, 2]
      omnibus <- TRUE
    }
    lldiff <- abs(ll0 - ll1) * 2
    dfdiff <- abs(df0 - df1)
    ps <- pchisq(q = lldiff, df = dfdiff, lower.tail = FALSE)
    decision <- c()
    for(i in seq_along(ps)){
      if(ps[i] <= alpha){
        decision[i] <- ifelse(ll0[i] > ll1[i], nn[1], nn[2])
      } else if(ps[i] == 1){
        decision[i] <- "- "
      } else if(ps[i] > alpha){
        if(omnibus){
          decision[i] <- ifelse(df0 < df1, nn[1], nn[2])
        } else {
          decision[i] <- ifelse(df0[i] > df1[i], nn[1], nn[2])
        }
      }
    }
    if(!omnibus){
      if(!is.null(d)){ps <- round(ps, d)}
      out <- data.frame(LL_diff2 = lldiff, Df_diff = dfdiff, pval = ps, decision = decision)
      rownames(out) <- rownames(object[[1]])
    } else {
      RMSEA <- function(X2, df, N){
        rmsea <- sqrt(max(c(((X2/N)/df) - (1/N), 0)))
        lower.l <- function(lambda){(pchisq(q = X2, df = df, ncp = lambda) - .95)}
        lambda.l <- tryCatch({uniroot(f = lower.l, lower = 0, upper = X2)$root}, error = function(e){0})
        rmsea.lower <- sqrt(lambda.l/(N * df))
        upper.l <- function(lambda){(pchisq(q = X2, df = df, ncp = lambda) - .05)}
        lambda.u <- tryCatch({uniroot(f = upper.l, lower = 0, upper = max(N, X2 * 4))$root}, error = function(e){1})
        rmsea.upper <- sqrt(lambda.u/(N * df))
        rmsea.pvalue <- 1 - pchisq(q = X2, df = df, ncp = (N * df * (.05^2)))
        return(c(lower = rmsea.lower, RMSEA = rmsea, upper = rmsea.upper, p.value = rmsea.pvalue))
      }
      rmsea <- RMSEA(lldiff, dfdiff, N)
      if(!is.null(d)){
        ps <- round(ps, d)
        lldiff <- round(lldiff, d)
        rmsea <- round(rmsea, d)
      }
      if(object[2, 2] < object[1, 2]){object <- object[order(object[, 2]), ]}
      out0 <- data.frame(LL_diff2 = c("", lldiff), Df_diff = c("", dfdiff),
                         pval = c("", ps), decision = c("", decision))
      out <- data.frame(object[, 1:2], out0, object[, 3:4])
      attr(out, "RMSEA") <- rmsea
    }
    return(out)
  }
  nn <- paste0("net", 0:1)
  if(length(net0) == 2 & ifelse(
    is.null(net1), TRUE, ifelse(is.logical(net1), net1, FALSE))){
    if(isTRUE(net1)){nodes <- net1}
    nn <- ifelse(!is.null(names(net0)), list(names(net0)), list(nn))[[1]]
    net1 <- net0[[2]]; net0 <- net0[[1]]
  }
  if(isTRUE(attr(net0, "mlGVAR"))){net0 <- net0$fixedNets}
  if("SURnet" %in% names(net0)){
    yLL <- isTRUE('SURll' %in% names(net0))
    omni0 <- switch(2 - yLL, net0$SURll$omnibus, omni_sll(fit = net0, s = s))
    uni0 <- switch(2 - yLL, net0$SURll$nodes, uni_sll(fit = net0))
    if(is.logical(net1)){nodes <- net1; net1 <- NULL}
    if(!is.null(net1)){
      if(isTRUE(attr(net1, "mlGVAR"))){net1 <- net1$fixedNets}
      if(is.null(lrt)){lrt <- TRUE}
      yLL <- isTRUE('SURll' %in% names(net1))
      omni1 <- switch(2 - yLL, net1$SURll$omnibus, omni_sll(fit = net1, s = s))
      uni1 <- switch(2 - yLL, net1$SURll$nodes, uni_sll(fit = net1))
      omni0 <- rbind(omni0, omni1)
      uni0 <- list(uni0, uni1)
      rownames(omni0) <- names(uni0) <- nn
    }
  } else {
    if(all(sapply(net0, function(z) isTRUE(attr(z, "mlGVAR"))))){
      net0 <- lapply(net0, '[[', "fixedNets")
    }
    stopifnot(all(sapply(net0, function(z) "SURnet" %in% names(z))))
    yLL <- all(sapply(net0, function(z) 'SURll' %in% names(z)))
    if(is.null(lrt)){lrt <- FALSE}
    if(is.logical(net1) & !sysfits){nodes <- net1}
    omni0 <- t(switch(2 - yLL, sapply(lapply(net0, '[[', 'SURll'), '[[', 'omnibus'), sapply(net0, omni_sll)))
    if(!is.null(orderBy)){
      orderBy <- switch(match.arg(tolower(as.character(
        orderBy)), c("true", "ll", "loglik", "aic", "bic", "df", "rank")),
        ll =, loglik = logLik, aic = AIC, bic = BIC, TRUE)
      if(is.function(orderBy)){
        orderBy <- c('LL', 'AIC', 'BIC')[which(sapply(
          list(logLik, AIC, BIC), function(z) identical(z, orderBy)))]
      }
      net0 <- net0[order(sapply(net0, function(z){
        if(isTRUE(orderBy)){
          return(sum(sapply(lapply(z$SURnet$mods, '[[', 'model'), nrow)))
        } else {
          return(omni0[, orderBy])
        }
      }), decreasing = decreasing)]
    }
    nn <- ifelse(!is.null(names(net0)), list(names(net0)),
                 list(paste0("net", 0:(length(net0) - 1))))[[1]]
    uni0 <- switch(2 - yLL, lapply(lapply(net0, '[[', 'SURll'), '[[', 'nodes'),
                   lapply(net0, uni_sll))
    rownames(omni0) <- names(uni0) <- nn
    sysgo <- all(sapply(net0, function(z) 'SURfit' %in% names(z)))
    if(sysfits & sysgo){
      net00 <- lapply(lapply(net0, '[[', "SURfit"), summary)
      ssr <- colSums(sapply(uni0, '[[', "RSS"))
      detrc <- sapply(net00, '[[', "detResidCov")
      olsr2 <- sapply(net00, '[[', "ols.r.squared")
      mcelroy <- sapply(net00, '[[', "mcelroy.r.squared")
      omni0 <- cbind(omni0, SSR = ssr, detSigma = detrc,
                     OLS.R2 = olsr2, McElroy.R2 = mcelroy)
      if(!is.null(d)){omni0 <- round(omni0, d)}
      return(omni0)
    }
  }
  out <- ifelse(all, list(list(nodes = uni0, omnibus = omni0)),
                ifelse(nodes, list(uni0), list(omni0)))[[1]]
  if(ifelse(!is.null(net1), lrt, FALSE)){
    lrt0 <- sur_lrt(object = omni0, d = d, alpha = alpha,
                    N = prod(dim(net1$SURnet$data$Y)))
    lrt1 <- sur_lrt(object = uni0, d = d, alpha = alpha)
    out <- list(nodes = uni0, LRT = lrt1, omnibus = lrt0)
    if(!all){out <- ifelse(nodes, list(lrt1), list(lrt0))[[1]]}
  }
  return(out)
}

##### SURtable: obtain all possible LRTs (with RMSEAs) comparing a list of models
#' Obtain all possible LRTs comparing a list of SUR models
#'
#' Includes RMSEAs
#'
#' @param fits list of networks
#' @param nodes logical
#' @param orderBy logical or character?
#' @param d numeric
#' @param alpha numeric
#' @param decreasing logical
#' @param names logical or character?
#' @param rmsea logical
#' @param s character
#'
#' @return Table of LRTs
#' @export
#'
#' @examples
#' 1 + 1
SURtable <- function(fits, nodes = FALSE, orderBy = TRUE, d = 4, alpha = .05,
                     decreasing = TRUE, names = NULL, rmsea = FALSE, s = "res"){
  n <- length(fits)
  stopifnot(is.list(fits) & n > 2)
  if(!is.null(names)){names(fits) <- names}
  if(is.null(names(fits))){names(fits) <- paste0("fit", 1:n)}
  if(all(sapply(fits, function(z) isTRUE(attr(z, "mlGVAR"))))){
    fits <- lapply(fits, '[[', "fixedNets")
  }
  stopifnot(all(sapply(fits, function(z) "SURnet" %in% names(z))))
  if(!is.null(orderBy)){
    oms <- SURll(fits)
    orderBy <- switch(match.arg(tolower(as.character(
      orderBy)), c("true", "ll", "loglik", "aic", "bic", "df", "rank")),
      ll =, loglik = logLik, aic = AIC, bic = BIC, TRUE)
    if(is.function(orderBy)){
      orderBy <- c('LL', 'AIC', 'BIC')[which(sapply(
        list(logLik, AIC, BIC), function(z) identical(z, orderBy)))]
    }
    fits <- fits[order(sapply(fits, function(z){
      if(isTRUE(orderBy)){
        return(sum(sapply(lapply(z$SURnet$mods, '[[', 'model'), nrow)))
      } else {
        return(oms[, orderBy])
      }
    }), decreasing = decreasing)]
  }
  tt <- combn(n, 2)
  lls0 <- SURll(fits, s = s)
  lrts0 <- lapply(seq_len(ncol(tt)), function(i){
    SURll(fits[tt[, i]], d = d, alpha = alpha, s = s)})
  out1 <- t(sapply(lrts0, rownames))
  out2 <- t(sapply(lrts0, function(z) as.numeric(z[2, 3:5])))
  out3 <- sapply(lrts0, function(z) z[2, 6])
  out4 <- cbind.data.frame(out1, out2, out3)
  colnames(out4) <- c("net0", "net1", "Chisq", "Df", "pval", "decision")
  rmsea0 <- data.frame(out4[, 1:2], do.call(rbind, lapply(lrts0, attr, "RMSEA")))
  select <- table(out4$decision)
  if(any(names(select) == "- ")){select <- select[names(select) != "- "]}
  lls0 <- cbind(lls0, LRT = numeric(nrow(lls0)))
  lls0[match(names(select), rownames(lls0)), "LRT"] <- unname(select)
  out <- list(LRT = out4, omnibus = lls0, RMSEA = rmsea0)
  if(!rmsea){out$RMSEA <- NULL}
  if(nodes != FALSE){
    lls1 <- SURll(fits, nodes = TRUE, s = s)
    stopifnot(length(unique(lapply(lls1, rownames))) == 1)
    nodenames <- unique(lapply(lls1, rownames))[[1]]
    lls2 <- lapply(nodenames, function(fit){
      nodemod <- matrix(unlist(sapply(lls1, function(node){
        node[fit, -1]})), nrow = length(lls1), ncol = 4, byrow = TRUE)
      dimnames(nodemod) <- list(names(fits), c("LL", "df", "AIC", "BIC"))
      return(nodemod)
    })
    names(lls2) <- nodenames
    lrts1 <- lapply(seq_len(ncol(tt)), function(i){
      SURll(fits[tt[, i]], nodes = TRUE, d = d, alpha = alpha, s = s)})
    netLRTs <- lapply(colnames(lrts1[[1]]), function(nn){
      n1 <- data.frame(out4[, 1:2], "|", do.call(rbind, lapply(lrts1, '[[', nn)))
      colnames(n1)[-c(1:2)] <- c("|", nodenames)
      return(n1)
    })
    names(netLRTs) <- colnames(lrts1[[1]])
    deci <- netLRTs$decision[, -c(1:3)]
    deci2 <- do.call(cbind.data.frame, lapply(1:ncol(deci), function(z){
      z <- table(deci[, z])
      if(any(names(z) == "- ")){z <- z[names(z) != "- "]}
      zz <- setNames(numeric(length(fits)), names(fits))
      zz[match(names(z), names(zz))] <- unname(z)
      return(zz)
    }))
    colnames(deci2) <- colnames(deci)
    out <- list(nodes = lls2, LRT = netLRTs, counts = deci2)
    if(is.character(nodes)){
      out$decision <- netLRTs$decision
      out$LRT <- NULL
    }
  }
  attr(out, "alpha") <- alpha
  return(out)
}

##### SURnet: create temporal and contemporaneous network of SUR results
#' Creates temporal and contemporaneous network of SUR results
#'
#' Used after something
#'
#' @param fit network
#' @param dat data
#' @param s character
#' @param m numeric
#' @param threshold logical
#' @param mval numeric
#' @param medges numeric
#' @param pcor character
#'
#' @return Temporal and contemporaneous networks
#' @export
#'
#' @examples
#' 1 + 1
SURnet <- function(fit, dat, s = "sigma", m = NULL, threshold = FALSE,
                   mval = NULL, medges = 1, pcor = "none"){
  y <- dat$Y
  p <- ncol(y)
  fitobj <- fit$eq
  yhat <- predict(fit)
  attr(fitobj, "rank") <- fit$rank
  ynames <- gsub("[.]y$", "", colnames(y))
  beta <- getCoefs(fit = fit, mat = "beta", data = dat)
  pvals <- getCoefs(fit = fit, mat = "pvals", data = dat)
  mods <- lapply(1:nrow(beta), function(z){
    model <- as.matrix(beta[z, ], ncol = 1)
    deviance <- sum((y[, z] - yhat[, z])^2)
    s <- sqrt(deviance/nrow(y))
    LL_model <- sum(dnorm(y[, z], mean = yhat[, z], sd = s, log = TRUE))
    k <- nrow(model) + 1
    aic <- (2 * k) - (2 * LL_model)
    bic <- (log(nrow(y)) * k) - (2 * LL_model)
    out <- list(deviance = deviance, LL_model = LL_model, AIC = aic, BIC = bic, model = model)
    return(out)
  })
  mods0 <- lapply(mods, '[[', "model")
  names(mods) <- names(mods0) <- names(fitobj) <- colnames(y)
  s <- match.arg(tolower(s), choices = c("sigma", "res", "dfres"))
  pcor <- ifelse(is.logical(pcor), "none", pcor)
  call <- list(type = rep("g", ncol(dat$Y)), moderators = m, mval = mval,
               lags = 1, residMat = s, threshold = threshold, pcor = pcor)
  call$mval <- mval <- ifelse(length(m) != 1, list(NULL), ifelse(
    attr(fit, "exogenous"), list(mval), list(NULL)))[[1]]
  if(!is.null(m)){
    mname <- call$moderators <- attr(fit, "moderators")
    exogenous <- call$exogenous <- attr(fitobj, "exogenous") <- attr(fit, "exogenous")
    intnames <- colnames(beta)[grep(":", colnames(beta))]
    beta2 <- beta[, intnames, drop = FALSE]
    pvals2 <- pvals[, intnames, drop = FALSE]
    interactions <- list(beta = beta2, pvals = pvals2)
    modEdges <- ifelse(pvals2 <= ifelse(!is.numeric(threshold), .05, threshold), 1, 0)
    #if(any(pvals2 == 0)){modEdges <- modEdges * ifelse(pvals2 == 0, 0, 1)}
    modEdges <- modEdges + 1
    if(length(m) == 1){
      if(is.character(m)){
        m <- which(colnames(dat$X) %in% m)
      } else if("covariates" %in% names(attributes(fit))){
        m <- m - sum(attr(fit, "covariates") < m)
      }
      ints <- lapply(mods0, function(z) rownames(z)[z[, 1] != 0][grep(":", rownames(z)[z[, 1] != 0])])
      inds0 <- unlist(lapply(ints, function(z) gsub(paste0(mname, ":|:", mname), "", z)))
      inds1 <- cbind(y = rep(names(mods), sapply(ints, length)), ints = unname(unlist(ints)))
      vars <- lapply(fitobj, vcov)
      interactions$coefvars <- data.frame(
        Y = inds1[, 1], X = unname(inds0), Z = mname, Int = inds1[, 2],
        t(sapply(seq_len(nrow(inds1)), function(i){
          vb1 <- vars[[inds1[i, 1]]][inds0[i], inds0[i]]
          vb3 <- vars[[inds1[i, 1]]][inds1[i, 2], inds1[i, 2]]
          vb1b3 <- vars[[inds1[i, 1]]][inds0[i], inds1[i, 2]]
          return(c(varX = vb1, varInt = vb3, varCov = vb1b3))
        }))
      )
      vars0 <- as.matrix(interactions$coefvars[, 1:4])
      vars1 <- interactions$coefvars[, -c(1:4)]
      if(!exogenous){
        modEdges0 <- matrix(1, p, p)
        modEdges0[, -m] <- modEdges
        modEdges0[, m] <- unname(apply(modEdges, 1, function(z) ifelse(any(z == 2), medges, 1)))
        rownames(modEdges0) <- rownames(modEdges)
        mcols0 <- character(p)
        mcols0[m] <- paste(rep(attr(fit, "moderators"), 2), collapse = ":")
        mcols0[-m] <- colnames(modEdges)
        colnames(modEdges0) <- mcols0
        modEdges <- modEdges0
      } else if(!is.null(mval)){
        margSE <- function(x, vars){
          as.numeric(sqrt(vars[1] + ((x^2) * vars[2]) + (2 * x * vars[3])))}
        ses <- getCoefs(fit = fit, mat = "ses", data = dat)
        for(i in 1:nrow(vars1)){
          B3Z <- mval * beta[vars0[i, "Y"], vars0[i, "Int"]]
          beta[vars0[i, "Y"], vars0[i, "X"]] <- beta[vars0[i, "Y"], vars0[i, "X"]] + B3Z
          ses[vars0[i, "Y"], vars0[i, "X"]] <- margSE(mval, vars1[i, ])
        }
        dfs <- matrix(rep(sapply(fitobj, '[[', "df.residual"), each = ncol(beta)),
                      nrow = nrow(beta), ncol = ncol(beta), byrow = TRUE)
        pvals <- (2 * pt(abs(beta/ses), df = dfs, lower.tail = FALSE))
        if(any(is.na(pvals))){pvals[is.na(pvals)] <- 1}
      }
    }
  } else {
    modEdges <- matrix(1, p, p)
  }
  b <- beta[, ynames, drop = FALSE]
  kappa <- solve(getCoefs(fit = fit, mat = s))
  PCC <- getCoefs(fit = fit, mat = "pcor")
  PDC <- b/(sqrt(diag(solve(kappa)) %o% diag(kappa) + b^2))
  pvals3 <- psych::corr.p(r = PCC, n = nrow(y), adjust = pcor)[[4]]
  pvals4 <- matrix(0, ncol(PCC), ncol(PCC))
  pvals4[upper.tri(pvals4)] <- pvals3[upper.tri(pvals3)]
  pvals4 <- as.matrix(Matrix::forceSymmetric(pvals4))
  dimnames(PCC) <- dimnames(kappa) <- dimnames(pvals4) <- dimnames(PDC)
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    colnames(colMat) <- paste0(colnames(adjMat), ifelse(
      !any(grepl("lag", colnames(adjMat))), ".lag1.", ""))
    rownames(colMat) <- rownames(adjMat)
    return(colMat)
  }
  if(threshold != FALSE){
    if(!is.character(threshold)){
      if(isTRUE(threshold)){threshold <- .05}
      b <- b * ifelse(pvals[, ynames, drop = FALSE] <= threshold, 1, 0)
      PDC <- PDC * ifelse(pvals[, ynames, drop = FALSE] <= threshold, 1, 0)
      if(!is.null(m)){interactions$beta <- beta2 * ifelse(pvals2 <= threshold, 1, 0)}
    }
    kdiag <- diag(kappa)
    kappa <- kappa * ifelse(pvals4 <= ifelse(!is.numeric(threshold), .05, threshold), 1, 0)
    PCC <- PCC * ifelse(pvals4 <= ifelse(!is.numeric(threshold), .05, threshold), 1, 0)
    diag(kappa) <- kdiag
  }
  temporal <- list(adjMat = b, edgeColors = getEdgeColors(b), modEdges = modEdges,
                   PDC = list(adjMat = PDC, edgeColors = getEdgeColors(PDC)),
                   coefs = list(beta = beta, pvals = pvals))
  if(length(m) != 1){temporal$modEdges <- NULL}
  colnames(temporal$adjMat) <- colnames(temporal$PDC$adjMat) <- colnames(temporal$edgeColors)
  contemporaneous <- list(adjMat = PCC, edgeColors = getEdgeColors(PCC), pvals = pvals4, kappa = kappa)
  colnames(contemporaneous$adjMat) <- colnames(contemporaneous$edgeColors)
  surNet <- list(call = call, temporal = temporal, contemporaneous = contemporaneous)
  if(!is.null(m)){surNet$interactions <- interactions}
  if(length(m) == 1 & ifelse(is.null(m), FALSE, exogenous)){
    madj <- modEdges0 <- matrix(0, p + 1, p + 1)
    madj[-m, -m] <- surNet$temporal$adjMat
    dimnames(madj) <- rep(list(1:(p + 1)), 2)
    rownames(madj)[-m] <- rownames(surNet$temporal$adjMat)
    rownames(madj)[m] <- paste0(mname, ".y")
    colnames(madj)[-m] <- colnames(surNet$temporal$adjMat)
    colnames(madj)[m] <- paste0(mname, ".lag1.")
    if(threshold != FALSE & !is.character(threshold)){
      madj[-m, m] <- beta[, mname] * ifelse(pvals[, mname] <= threshold, 1, 0)
    } else {
      madj[-m, m] <- beta[, mname]
    }
    modEdges0[-m, -m] <- modEdges
    modEdges0[-m, m] <- unname(apply(
      modEdges, 1, function(z){ifelse(any(z == 2), medges, 1)}))
    dimnames(modEdges0) <- dimnames(madj)
    shape <- rep("circle", p + 1)
    shape[m] <- "square"
    surNet$mnet <- list(adjMat = madj, edgeColors = getEdgeColors(madj),
                        modEdges = modEdges0, shape = shape)
  }
  surNet <- append(surNet, list(mods = mods, data = dat))
  return(surNet)
}

##### surEqs: create regression equations for fitting SUR models
surEqs <- function(data, varMods = NULL, mod = "min", m = NULL,
                   exogenous = TRUE, covs = NULL){
  if(!is.null(varMods)){
    if(attr(varMods, "criterion") != "CV"){mod <- "min"}
    mod <- match.arg(mod, c("min", "1se"))
    mod <- ifelse(mod == "min" | attr(varMods, "method") == "regsubsets", "mod0", "mod1se")
    x <- lapply(varMods, '[[', mod)
    y <- names(x)
    if("covs" %in% names(attributes(varMods))){covs <- attr(varMods, "covs")}
    if(!is.null(covs) & !is(data, 'list')){
      covs <- colnames(data)[covs]
      x <- lapply(x, function(z) c(z[!grepl(":", z)], covs, z[grepl(":", z)]))
    }
  } else {
    if(!is(data, 'list')){
      data <- lagMat(data = data, m = m, exogenous = exogenous,
                     covariates = covs, checkType = TRUE)
    }
    y <- colnames(data$Y)
    x <- rep(list(colnames(data$X)), length(y))
  }
  eqs <- lapply(seq_along(y), function(z){
    as.formula(paste(y[z], "~", paste(x[[z]], collapse = " + ")))})
  return(eqs)
}

##### getCoefs: extract beta matrix and/or sigma matrix from systemfit model
getCoefs <- function(fit, mat = "beta", data = NULL){
  if("SURfit" %in% names(fit)){fit <- fit$SURfit}
  if(class(fit) == "systemfit"){
    mat <- match.arg(arg = mat, c(
      "beta", "pvals", "ses", "sigma", "res", "dfres", "cor", "pcor"))
    if(mat %in% c("beta", "pvals", "ses")){
      ynames <- c(); n <- list()
      for(i in 1:length(fit$eq)){
        ynames[i] <- as.character(fit$eq[[i]]$terms[[2]])
        n[[i]] <- names(coef(fit$eq[[i]]))
      }
      if(!is.null(data)){
        if(class(data) == "list"){data <- data$X}
        if(is.null(colnames(data))){colnames(data) <- paste0("X", 1:ncol(data))}
        bnames <- c("(Intercept)", colnames(data))
        N <- length(bnames)
      } else {
        N <- ifelse(any(grep(":", names(coef(fit)))), length(fit$eq) * 2, length(fit$eq) + 1)
        bnames <- n[[which.max(sapply(n, length))]]
      }
      b <- t(matrix(ifelse(mat == "beta", 0, 1), ncol = length(fit$eq), nrow = N))
      rownames(b) <- ynames
      colnames(b) <- bnames
      for(i in 1:nrow(b)){
        bb <- ifelse(
          mat == "beta", list(coef(fit$eq[[i]])), ifelse(
            mat == "pvals", list(summary(fit$eq[[i]])$coefficients[, 4]),
            list(summary(fit$eq[[i]])$coefficients[, 2])))[[1]]
        for(j in 1:length(n[[i]])){b[i, names(bb)[j]] <- bb[j]}
      }
    } else {
      if(mat %in% c("res", "dfres")){
        e <- as.matrix(residuals(fit))
        N <- ifelse(mat == "res", nrow(e), nrow(e) - 1)
        b <- (t(e) %*% e)/N
      } else if(mat == "pcor"){
        b <- tryCatch({-solve(cor(residuals(fit)))}, error = function(e){
          -corpcor::pseudoinverse(cor(residuals(fit)))})
        diag(b) <- -diag(b)
        delta <- (1/sqrt(diag(b)))
        b <- t(delta * b) * delta
        diag(b) <- 0
      } else {
        b <- fit$residCov
        if(mat == "cor"){b <- cov2cor(b)}
      }
    }
  } else {
    b <- do.call(cbind, sapply(fit$mods, '[', "model"))
    colnames(b) <- rownames(fit$adjMat)
  }
  return(b)
}

##### SURsampler: Sample data for fitting system of lagged SURs
SURsampler <- function(B = NULL, S, n, seed = NULL, beta, beta2 = NULL,
                       full = TRUE, mu = 0, cholesky = FALSE,
                       time = FALSE, allDat = TRUE){
  p <- ncol(S)
  if(!is.null(B)){
    if(is(B, 'list')){B <- do.call(rbind, B)}
    if(missing(beta)){
      if(ncol(B) == (p + 1)){
        beta <- B
      } else {
        beta <- B[, 1:(p + 1)]
        beta2 <- B[, -c(1:(p + 1)), drop = FALSE]
      }
    }
  }
  if(is(beta, 'list')){beta <- do.call(rbind, beta)}
  if(dim(beta)[2] != dim(beta)[1] + 1){
    if(dim(beta)[2] == dim(beta)[1]){
      beta <- cbind(0, beta)
    } else {stop("Invalid dimensions for beta matrix")}
  }
  if(!is.null(beta2)){
    if(is(beta2, 'list')){beta2 <- do.call(rbind, beta2)}
    stopifnot(dim(beta2)[2] == dim(beta2)[1] - 1)
  }
  if(length(mu) == 1){mu <- rep(mu, p)}
  if(!is.null(seed)){set.seed(seed)}
  tx <- Sys.time()
  if(!cholesky){
    R <- mvtnorm::rmvnorm(n = n, mean = mu, sigma = S)
    #R <- MASS::mvrnorm(n = n, mu = mu, Sigma = S)
  } else {
    R <- matrix(rnorm(p * n), ncol = p) %*% chol(S)
  }
  X <- matrix(NA, ncol = p, nrow = n + 1)
  if(!is.null(seed)){set.seed(seed)}
  X[1, ] <- rnorm(p)
  it <- 0
  if(is.null(beta2)){
    repeat{
      it <- it + 1
      for(i in 1:p){
        X[it + 1, i] <- beta[i, 1] + sum(X[it, ] * beta[i, -1]) + R[it, i]
      }
      if(it == n){break}
    }
    Y <- X[-1, ]
  } else {
    Y <- X
    X2 <- X[, 1] * X[, -1]
    X <- cbind(X, X2)
    beta <- cbind(beta, beta2)
    repeat{
      it <- it + 1
      for(i in 1:p){
        Y[it + 1, i] <- beta[i, 1] + sum(X[it, ] * beta[i, -1]) + R[it, i]
        X <- cbind(Y, (Y[, 1] * Y[, -1]))
      }
      if(it == n){break}
    }
    Y <- Y[-1, ]
  }
  tx <- Sys.time() - tx
  if(time){cat("Completed in", round(tx, 3), attr(tx, "units"), "\n")}
  colnames(Y) <- paste0("X", 1:ncol(Y), ".y")
  X <- X[-nrow(X), ]
  colnames(X) <- paste0("X", 1:ncol(X))
  if(!is.null(beta2)){
    ints <- (ncol(beta) - ncol(beta2)):ncol(X)
    for(i in ints){colnames(X)[i] <- paste0("X1:X", which(ints == i) + 1)}
  }
  dat <- list(Y = Y, X = X)
  if(full){
    dat <- list(Y = Y, X = X, full = data.frame(do.call(cbind, dat)))
    if(!is.null(beta2)){colnames(dat$full)[ncol(dat$Y) + ints] <- colnames(dat$X)[ints]}
  }
  if(time){attributes(dat)$time <- tx}
  if(!allDat){
    X <- dat$X[, 1:p]
    X <- rbind(X, dat$Y[nrow(dat$Y), ])
    dat <- data.frame(X)
  }
  dat
}
