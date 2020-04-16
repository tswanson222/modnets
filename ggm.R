### ======================================================================== ###
### ========================= PRIMARY GGM FUNCTIONS ======================== ###
### ======================================================================== ###
##### nodewise: nodewise regression, with option to simulate/add moderators
nodewise <- function(data, mods = NULL, varMods = NULL, lambda = "min", center = TRUE, 
                     scale = FALSE, covariates = NULL, exogenous = TRUE, getEqs = FALSE){
  data <- dat <- dat0 <- data.frame(data)
  vs <- colnames(data)
  if(!is.null(varMods)){
    if(class(varMods) == "list"){
      if(all(sapply(varMods, class) != "list")){
        varMods <- lapply(varMods, list)
        lambda <- "min"
      }
      if(is.null(names(varMods))){
        if((is.null(mods) | (!is.null(mods) & class(mods) == "list")) & 
           (is.null(covariates) | (!is.null(covariates) & class(covariates) == "list"))){
          names(varMods) <- colnames(data)
        } else {
          if(!is.null(mods) & class(mods) %in% c("numeric", "integer")){mc <- mods} else {mc <- NULL}
          if(!is.null(covariates) & class(covariates) %in% c("numeric", "integer")){mc <- c(mc, covariates)}
          names(varMods) <- colnames(data)[-mc]
        }
      }
      if(!is.null(mods)){if(length(mods) > 1){mods <- NULL}}
      if(length(varMods) == ncol(data)){exogenous <- FALSE}
      if(is.null(unlist(lapply(varMods, '[[', "mod1se")))){lambda <- "min"}
      lambda <- ifelse(lambda == "min", 1, 2)
      type <- unname(sapply(varMods, attr, "family"))
      if(all(type == "c")){center <- FALSE}
    } else {
      type <- varMods
      if(is.null(mods)){
        varMods <- NULL
      } else if(length(mods) == 1){
        varMods <- NULL
      } else {
        varMods <- lapply(1:ncol(data), function(z){
          mains <- vs[-z]
          ints <- apply(combn(mains, 2), 2, paste0, collapse = ":")
          if(!z %in% mods){ints <- ints[apply(combn(mains, 2), 2, function(zz) any(vs[mods] %in% zz))]}
          return(list(mod0 = c(mains, ints)))
        })
        names(varMods) <- vs
        exogenous <- FALSE
        mods <- NULL
        lambda <- 1
      }
    }
  } else if(!is.null(mods)){
    if(length(mods) > 1){
      varMods <- lapply(1:ncol(data), function(z){
        mains <- vs[-z]
        ints <- apply(combn(mains, 2), 2, paste0, collapse = ":")
        if(!z %in% mods){ints <- ints[apply(combn(mains, 2), 2, function(zz) any(vs[mods] %in% zz))]}
        return(list(mod0 = c(mains, ints)))
      })
      names(varMods) <- vs
      type <- unname(ifelse(apply(data, 2, function(z) dim(table(z))) <= 2, "c", "g"))
      exogenous <- FALSE
      mods <- NULL
      lambda <- 1
    }
  }
  intMatrix <- function(data, mods = NULL, covariates = NULL){
    if(class(mods) == "list"){stopifnot(!is.null(names(mods)))}
    data <- as.data.frame(data)
    vars <- colnames(data)
    eqs <- list()
    if(!is.null(mods)){
      for(i in 1:length(vars)){
        modTerms <- list()
        for(j in 1:length(mods)){modTerms[[j]] <- paste0(names(mods[j]), " + ", paste(paste0(vars[-i], ":", names(mods[j])), collapse = " + "))}
        modTerms <- paste(modTerms, collapse = " + ")
        eqs[[i]] <- paste0(vars[i], " ~ ", paste(vars[-i], collapse = " + "), " + ", modTerms)
      }
    } else {
      for(i in 1:length(vars)){eqs[[i]] <- paste0(vars[i], " ~ ", paste(vars[-i], collapse = " + "))}
    }
    if(!is.null(covariates)){
      for(i in 1:length(vars)){eqs[[i]] <- paste0(eqs[[i]], " + ", paste(names(covariates), collapse = " + "))}
    }
    return(lapply(eqs, as.formula))
  }
  if(!is.null(covariates)){
    if(class(covariates) %in% c("numeric", "integer")){
      if(length(covariates) > 1){
        covs <- as.list(data[, covariates])
      } else {
        covs <- list(data[, covariates])
        names(covs) <- colnames(data)[covariates]
      }
      data <- dat <- dat0 <- data[, -covariates]
      covariates <- covs
    }
  }
  if(!is.null(mods)){
    if(class(mods) %in% c("numeric", "integer")){
      mod <- list(data[, mods])
      names(mod) <- colnames(data)[mods]
      data <- data[, -mods]
      if(length(type) > ncol(data)){type <- c(type[-mods], type[mods])}
      mods <- mod
    }
    if(length(mods) != 1){stop("Cannot specify more than on exogenous moderator")}
    dat <- dat0 <- data.frame(data, mods)
  }
  if(center != FALSE){
    binary <- unname(which(apply(dat, 2, function(z) dim(table(z)) <= 2)))
    if(length(binary) == ncol(dat) | ifelse(!is.null(mods), length(binary) == ncol(dat) - 1, FALSE)){
      type <- "binomial"
    }
    if(!is.null(mods) & dim(table(mods[[1]])) <= 2 | center != TRUE){
      if(length(binary) > 0){
        dat[, -union(binary, ncol(dat))] <- apply(dat[, -union(binary, ncol(dat))], 2, scale, TRUE, scale)
      } else {
        dat[, -ncol(dat)] <- apply(dat[, -ncol(dat)], 2, scale, TRUE, scale)
      }
    } else {
      if(length(binary) > 0){
        dat[, -binary] <- apply(dat[, -binary], 2, scale, TRUE, scale)
      } else {
        dat <- apply(dat, 2, scale, TRUE, scale)
      }
    }
    if(!is.null(covariates)){
      covariates <- lapply(covariates, function(z){
        ifelse(dim(table(z)) <= 2 | center != TRUE, list(z), list(scale(z, TRUE, scale)))[[1]]
      })
    }
    dat <- dat0 <- data.frame(dat)
  }
  if(!is.null(covariates)){dat <- data.frame(dat, covariates)}
  if(!is.null(varMods)){
    ints <- as.list(paste(names(varMods), "~", lapply(varMods, function(z) paste0(z[[lambda]], collapse = " + "))))
    if(!is.null(covariates) & ifelse("covariates" %in% names(attributes(varMods)), FALSE, TRUE)){
      ints <- lapply(ints, paste0, " + ", paste(names(covariates), collapse = " + "))
    }
  } else {
    ints <- intMatrix(data, mods, covariates)
    if(!is.null(mods) & !exogenous){ints <- append(ints, list(as.formula(paste0(names(mods), " ~ .^2"))))}
  }
  if(getEqs){return(lapply(ints, as.formula))}
  if(exists("type", inherits = FALSE)){
    if(length(type) == 1){
      type <- rep(match.arg(type, c("g", "c", "gaussian", "binomial")), length(ints))
    }
    if(any(type %in% c("g", "c"))){
      type <- unname(sapply(type, switch, "g" = "gaussian", "c" = "binomial"))}
    m <- suppressWarnings(lapply(1:length(ints), function(z){
      #if(type[z] == "gaussian"){mglm <- lm(ints[[z]], dat)}
      #if(type[z] == "binomial"){mglm <- glm(ints[[z]], data = dat, family = type[z])}
      mglm <- glm(ints[[z]], data = dat, family = type[z])
      attr(mglm, "family") <- type[z]
      return(mglm)
    }))
  } else {
    m <- lapply(ints, function(z) lm(z, dat))
  }
  if(!is.null(mods) & !exogenous){
    data0 <- data
    data <- dat0
  }
  mm <- lapply(lapply(m, coef), function(z) z[which(names(z) %in% colnames(data))])
  ps1 <- lapply(m, function(z){
    z1 <- is.na(coef(z))
    z2 <- t(data.frame(summary(z)$coefficients))
    z0 <- ifelse(any(z1), list(names(z1[!z1])), list(names(coef(z))))[[1]]
    z3 <- z2[4, which(z0 %in% colnames(data)), drop = FALSE]
    rownames(z3) <- NULL
    return(z3[1, ])
  })
  if(!is.null(mods)){
    mname <- names(mods)
    m2 <- lapply(lapply(m, coef), function(z) z[grep(paste0(":", mname, "$"), names(z))])
    ps2 <- lapply(m, function(z){
      z2 <- t(data.frame(summary(z)$coefficients))
      z3 <- z2[4, grep(paste0(":", mname, "$"), colnames(z2)), drop = FALSE]
      rownames(z3) <- NULL
      return(z3[1, ])
    })
    vs <- ifelse(exogenous, list(colnames(data)), list(colnames(data0)))[[1]]
    mx <- paste0(vs, ":", mname)
    psx <- bx <- suppressWarnings(diag(mx))
    rownames(psx) <- rownames(bx) <- vs
    colnames(psx) <- colnames(bx) <- mx
    diag(psx) <- diag(bx) <- 1
    m3 <- lapply(m, function(z) summary(z)$coefficients)
    names(m3) <- vs
    bm <- do.call(rbind, lapply(m3, function(z) z[rownames(z) == mname, ]))
    if(!exogenous){
      m2 <- m2[-ncol(data)]
      ps2 <- ps2[-ncol(data)]
      m3 <- m3[-ncol(data)]
    }
  }
  if(!is.null(covariates)){
    m4 <- lapply(m, function(z) summary(z)$coefficients)
    names(m4) <- ifelse(!is.null(mods), list(vs), list(colnames(data)))[[1]]
    cnames <- names(covariates)
    bcovs <- list()
    for(i in 1:length(cnames)){
      bcovs[[i]] <- do.call(rbind, lapply(m4, function(z) z[rownames(z) == cnames[i], ]))
    }
    names(bcovs) <- cnames
  }
  b <- ps <- matrix(NA, nrow = ncol(data), ncol = ncol(data))
  for(i in 1:ncol(data)){
    b[i, match(names(mm[[i]]), colnames(data))] <- mm[[i]]
    ps[i, match(names(ps1[[i]]), colnames(data))] <- ps1[[i]]
    if(!is.null(mods)){
      if(exogenous | (!exogenous & i < ncol(data))){
        bx[i, match(names(m2[[i]]), mx)] <- m2[[i]]
        psx[i, match(names(ps2[[i]]), mx)] <- ps2[[i]]
      }
    }
  }
  diag(b) <- diag(ps) <- 1
  if(any(is.na(b))){b[is.na(b)] <- 0}
  if(any(is.na(ps))){ps[is.na(ps)] <- 0}
  out <- list(models = m, B = list(b = b, ps = ps))
  if(is.null(mods)){
    attributes(out$models)$noMods <- TRUE
  } else {
    attributes(out$models)$exogenous <- exogenous
    if(nrow(bm) >= 1 & exogenous){out$Bm <- bm}
    out$Bx <- list(bx = bx, px = psx)
  }
  if(!is.null(varMods)){
    if(!any(grepl(":", unlist(sapply(out$models, function(z) names(coef(z))))))){
      attributes(out$models)$noMods <- TRUE
      out$Bx <- NULL
    }
    attributes(out$models)$varMods <- c("min", "1se")[lambda]
  }
  if(!is.null(covariates)){
    out$dat <- dat0
    out$covariates <- list(covs = do.call(cbind.data.frame, covariates), Bcovs = bcovs)
  } else {
    out$dat <- dat
  }
  if(exists("type", inherits = FALSE)){attr(out, "type") <- type}
  out
}

##### modNet: create moderated network from nodewise regression models
modNet <- function(models, data = NULL, threshold = FALSE, rule = "AND", mval = NULL, 
                   pcor = FALSE, useCIs = FALSE, nsims = 5000, mlty = 2, binarize = FALSE){
  if("models" %in% names(models)){
    mods0 <- models
    models <- models$models
    if(is.null(data)){
      if("noMods" %in% names(attributes(models)) & length(models) == ncol(mods0$dat)){
        data <- mods0$dat
      } else if(attr(models, "exogenous") == TRUE){
        data <- mods0$dat[, -ncol(mods0$dat)]
      } else {
        data <- mods0$dat
      }
    }
  }
  p <- ncol(data)
  n <- nrow(data)
  vs <- colnames(data)
  mods <- lapply(models, function(z){
    z2 <- matrix(coef(z), ncol = 1)
    rownames(z2) <- names(coef(z))
    return(z2)
  })
  mods2 <- lapply(mods, function(z) z[which(rownames(z) %in% vs), ])
  pvals <- lapply(models, function(z){
    z1 <- is.na(coef(z))
    z2 <- t(data.frame(summary(z)$coefficients))
    z0 <- ifelse(any(z1), list(names(z1[!z1])), list(names(coef(z))))[[1]]
    z3 <- z2[4, which(z0 %in% vs), drop = FALSE]
    rownames(z3) <- NULL
    return(z3[1, ])
  })
  ses <- lapply(models, function(z){
    z1 <- is.na(coef(z))
    z2 <- t(data.frame(summary(z)$coefficients))
    z0 <- ifelse(any(z1), list(names(z1[!z1])), list(names(coef(z))))[[1]]
    z3 <- z2[2, which(z0 %in% vs), drop = FALSE]
    rownames(z3) <- NULL
    return(z3[1, ])
  })
  b <- matrix(0, p, p)
  pvals2 <- ses2 <- matrix(1, p, p)
  for(i in 1:p){
    b[i, match(names(mods2[[i]]), vs)] <- mods2[[i]]
    ses2[i, match(names(ses[[i]]), vs)] <- ses[[i]]
    pvals2[i, match(names(pvals[[i]]), vs)] <- pvals[[i]]
  }
  results <- lapply(1:p, function(z){
    notype <- !"type" %in% names(attributes(mods0))
    if(notype | ifelse(!notype, ifelse(attr(mods0, "type")[z] == "gaussian", TRUE, FALSE), TRUE)){
      yhat <- predict(models[[z]])
      deviance <- sum((data[, z] - yhat)^2)
      s <- sqrt(deviance/n)
      LL_model <- sum(dnorm(data[, z], mean = yhat, sd = s, log = TRUE))
      k <- nrow(mods[[z]]) + 1
      aic <- (2 * k) - (2 * LL_model)
      bic <- (log(n) * k) - (2 * LL_model)
    } else {
      deviance <- deviance(models[[z]])
      LL_model <- as.numeric(logLik(models[[z]]))
      aic <- AIC(models[[z]])
      bic <- BIC(models[[z]])
    }
    return(list(deviance = deviance, LL_model = LL_model, AIC = aic, BIC = bic, model = mods[[z]]))
  })
  if(!"noMods" %in% names(attributes(models))){pcor <- FALSE}
  if("noMods" %in% names(attributes(models))){
    inds <- ints <- mval <- vars1 <- intMats <- NULL
  } else if("varMods" %in% names(attributes(models))){
    nmods <- ifelse(attr(models, "exogenous") == TRUE, length(mods), length(mods) - 1)
    inds0 <- lapply(mods, function(z) rownames(z)[grepl(":", rownames(z))])[1:nmods]
    inds1 <- unlist(lapply(inds0, function(z) gsub(":.*", "", z)))
    inds <- data.frame(outcome = rep(1:nmods, sapply(inds0, length)), interaction = unlist(inds0))
    vars <- lapply(models, vcov)[1:nmods]
    vars1 <- ints <- list()
    for(i in 1:nrow(inds)){
      ints[[i]] <- mods[[inds[i, 1]]][inds[i, 2], 1]
      vars1[[i]] <- c(vars[[inds[i, 1]]][inds1[i], inds1[i]], vars[[inds[i, 1]]][inds[i, 2], inds[i, 2]], vars[[inds[i, 1]]][inds1[i], inds[i, 2]])
      names(vars1[[i]]) <- c(inds1[i], inds[i, 2], "cov")
    }
  } else {
    nmods <- ifelse(attr(models, "exogenous") == TRUE, length(mods), length(mods) - 1)
    inds <- t(combn(1:nmods, 2))
    ints <- vector("list", nrow(inds))
    for(i in 1:nrow(inds)){
      ints[[i]][1] <- mods[[inds[i, 1]]][nmods + inds[i, 2], ]
      ints[[i]][2] <- mods[[inds[i, 2]]][nmods + inds[i, 1] + 1, ]
    }
    vars <- lapply(models, vcov)[1:nmods]
    vars1 <- list()
    for(i in 1:nmods){
      vars1[[i]] <- vector("list", nmods - 1)
      for(j in 1:(nmods - 1)){
        vars1[[i]][[j]] <- c(vars[[i]][j + 1, j + 1], vars[[i]][j + nmods + 1, j + nmods + 1], vars[[i]][j + 1, j + nmods + 1])
        names(vars1[[i]][[j]]) <- c(colnames(vars[[i]])[j + 1], colnames(vars[[i]])[j + nmods + 1], "cov")
      }
      names(vars1[[i]]) <- colnames(vars[[i]])[2:nmods]
    }
    names(vars1) <- colnames(data)[1:nmods]
  }
  b2ggm <- function(b, rule = "AND", pcor = FALSE, threshold = FALSE, n = NULL){
    rule <- match.arg(tolower(rule), c("and", "or"))
    if(rule == "and" & pcor == FALSE){
      bb <- cbind(b[upper.tri(b)], t(b)[upper.tri(t(b))])
      notBoth <- !apply(bb, 1, function(z) (z[1] == 0 & z[2] == 0) | (z[1] != 0 & z[2] != 0))
      if(any(notBoth)){
        bb[notBoth, ] <- 0
        b[upper.tri(b)] <- bb[, 1]
        b <- t(b)
        b[upper.tri(b)] <- bb[, 2]
        b <- t(b)
      }
    }
    if(pcor != FALSE){
      bb <- sign(b) * sqrt(b * t(b))
      if(pcor == 'cor'){
        diag(bb) <- 1
        bb <- corpcor::pcor2cor(bb)
      } else if(grepl('cor_auto', pcor)){
        bb <- qgraph::cor_auto(data, npn.SKEPTIC = TRUE)
        if(pcor == 'cor_auto2'){bb <- corpcor::cor2pcor(bb)}
      }
      diag(bb) <- 0
      if(threshold != FALSE){
        pcor <- ifelse(isTRUE(pcor) | grepl('cor', pcor), 'none', pcor)
        if(is.character(threshold)){pcor <- threshold}
        if(!is.numeric(threshold)){threshold <- .05}
        dimnames(bb) <- rep(list(paste0("X", 1:ncol(bb))), 2)
        sigMat <- ifelse(psych::corr.p(bb, n, adjust = pcor)[[4]] <= threshold, 1, 0)
        sigMat0 <- matrix(0, ncol(bb), ncol(bb))
        sigMat0[upper.tri(sigMat0)] <- sigMat[upper.tri(sigMat)]
        sigMat0 <- as.matrix(Matrix::forceSymmetric(sigMat0))
        bb <- bb * sigMat0
        bb <- unname(bb)
      }
      return(bb)
    } else {
      return((b + t(b))/2)
    }
  }
  if(!is.null(mval)){
    margSE <- function(x, vars){sqrt(vars[1] + ((x^2) * vars[2]) + (2 * x * vars[3]))}
    if("varMods" %in% names(attributes(models))){
      inds1.1 <- match(inds1, vs)
      for(i in 1:length(ints)){
        b[inds[i, 1], inds1.1[i]] <- b[inds[i, 1], inds1.1[i]] + (mval * ints[[i]])
        ses[[inds[i, 1]]][inds1[i]] <- margSE(mval, vars1[[i]])
      }
      for(i in 1:p){ses2[i, match(names(ses[[i]]), vs)] <- ses[[i]]}
    } else {
      for(i in 1:length(ints)){
        b[inds[i,1], inds[i,2]] <- b[inds[i,1], inds[i,2]] + (mval * ints[[i]][1])
        b[inds[i,2], inds[i,1]] <- b[inds[i,2], inds[i,1]] + (mval * ints[[i]][2])
      }
      inds2 <- rbind(inds, cbind(inds[, 2], inds[, 1]))
      inds2 <- inds2[order(inds2[, 1]), ]
      inds3 <- cbind(inds2[, 1], rep(c(1:(p - 1)), p))
      for(i in 1:nrow(inds3)){ses[[inds3[i, 1]]][inds3[i, 2]] <- margSE(mval, vars1[[inds3[i, 1]]][[inds3[i, 2]]])}
      for(i in 1:p){ses2[i, match(names(ses[[i]]), vs)] <- ses[[i]]}
    }
    bmval <- b
    dimnames(bmval) <- rep(list(colnames(data)), 2)
    dfs <- matrix(sapply(models, function(z) z$df.residual), ncol = p, nrow = p)
    pvals2 <- (2 * pt(abs(b/ses2), df = dfs, lower.tail = FALSE))
    if(any(is.na(pvals2))){pvals2[is.na(pvals2)] <- 1}
  }
  if(threshold != FALSE & pcor == FALSE){
    if(threshold == TRUE){threshold <- .05}
    b <- b * ifelse(pvals2 <= threshold, 1, 0)
  }
  bb <- b2ggm(b, rule = rule, pcor = pcor, threshold = threshold, n = n)
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    colnames(colMat) <- colnames(adjMat)
    rownames(colMat) <- rownames(adjMat)
    colMat
  }
  if(!"noMods" %in% names(attributes(models))){
    if(useCIs & nsims > 0 & attr(models, "exogenous") == TRUE){
      cis <- margCIs(mods = mods0, alpha = ifelse(threshold == FALSE, .05, threshold), nsims = nsims)
      modEdgesNW <- getInts(x = cis, allInts = TRUE)
    } else {
      modEdgesNW <- ifelse(mods0$Bx$px <= ifelse(threshold == FALSE, .05, threshold), 1, 0)
      if(any(mods0$Bx$px == 0)){modEdgesNW <- modEdgesNW * ifelse(mods0$Bx$px == 0, 0, 1)}
      colnames(modEdgesNW) <- rownames(modEdgesNW)
    }
    modEdges <- t(modEdgesNW) * modEdgesNW
    modEdgesNW <- (modEdgesNW * (mlty - 1)) + 1
    modEdges <- (modEdges * (mlty - 1)) + 1
    rule <- match.arg(tolower(rule), c("and", "or"))
    if(rule == "or"){
      modEdges <- modEdgesNW + t(modEdgesNW)
      modEdges[modEdges == 2] <- 1
      modEdges[modEdges > 2] <- mlty
    }
    intMat1 <- mods0$Bx$bx
    intPs <- mods0$Bx$px
    if(any(intPs == 0)){intPs[intPs == 0] <- 1}
    if(threshold != FALSE){intMat1 <- intMat1 * ifelse(intPs <= threshold, 1, 0)}
    intMat2 <- b2ggm(intMat1, rule = rule, pcor = FALSE)
    diag(intMat1) <- diag(intMat2) <- 0
    intMats <- list(avgInts = t(intMat2), nwInts = list(adj2NW = t(intMat1), pvals2NW = t(intPs)))
  } else {
    modEdges <- modEdgesNW <- matrix(1, p, p)
  }
  diag(modEdges) <- diag(modEdgesNW) <- 0
  dimnames(b) <- dimnames(bb) <- dimnames(pvals2) <- rep(list(colnames(data)), 2)
  names(mods) <- names(models) <- names(results) <- colnames(data)
  out <- list(adjMat = bb, edgeColors = getEdgeColors(bb), modEdges = t(modEdges),
              nodewise = list(adjNW = t(b), edgeColsNW = getEdgeColors(t(b)), 
                              pvalsNW = t(pvals2), modEdgesNW = t(modEdgesNW)),
              interactions = list(intMats = intMats, inds = inds, ints = ints, coefvars = vars1))
  if(binarize){
    out$adjMat[out$adjMat != 0] <- 1
    out$edgeColors[out$edgeColors != 'darkgrey'] <- 'darkgrey'
  }
  if("noMods" %in% names(attributes(models))){
    out[c("interactions", "modEdges")] <- NULL
    out[["nodewise"]]["modEdgesNW"] <- NULL
  } else if(attr(models, "exogenous") == FALSE){
    out[["modEdges"]] <- rbind(cbind(out[["modEdges"]], 1), 1)
    out[["nodewise"]][["modEdgesNW"]] <- rbind(cbind(out[["nodewise"]][["modEdgesNW"]], 1), 1)
    colnames(out$modEdges)[p] <- rownames(out$modEdges)[p] <- colnames(mods0$dat)[p]
    colnames(out$nodewise$modEdgesNW)[p] <- rownames(out$nodewise$modEdgesNW)[p] <- colnames(mods0$dat)[p]
    attributes(out)$moderator <- colnames(mods0$dat)[ncol(mods0$dat)]
  } else {
    if(useCIs){out$interactions$cis <- cis}
    attributes(out)$moderator <- colnames(mods0$dat)[ncol(mods0$dat)]
    if("Bm" %in% names(mods0)){
      mp <- ncol(bb)
      madj <- medges <- matrix(0, mp + 1, mp + 1)
      md <- matrix(FALSE, mp + 1, mp + 1)
      madj[1:mp, 1:mp] <- bb
      dimnames(madj) <- dimnames(medges) <- dimnames(md) <- rep(list(colnames(mods0$dat)), 2)
      if(nrow(mods0$Bm) != p){mpp <- which(colnames(data) %in% rownames(mods0$Bm))} else {mpp <- 1:mp}
      if(threshold != FALSE){
        madj[(mp + 1), mpp] <- ifelse(mods0$Bm[, 4] <= threshold, mods0$Bm[, 1], 0)
      } else {
        madj[(mp + 1), mpp] <- mods0$Bm[, 1]
      }
      medges[1:mp, 1:mp] <- t(modEdges)
      medges[(mp + 1), 1:mp] <- 1
      md[(mp + 1), 1:mp] <- TRUE
      ints <- out$interactions
      out$interactions <- NULL
      out$mnet <- list(adjMat = madj, edgeColors = getEdgeColors(madj), modEdges = medges, d = md)
      out$interactions <- ints
    }
  }
  out$mods <- results
  out$fitobj <- models
  out$data <- data
  if(!is.null(mval)){
    attributes(out)$mval <- mval
    if(threshold != FALSE){out$nodewise[[paste0("adjMval", mval)]] <- t(bmval)}
  }
  attributes(out)$ggm <- TRUE
  out
}

##### modLL: log-likelihood of whole network, & LRT comparing two networks
modLL <- function(net0, net1 = NULL, nodes = FALSE, lrt = NULL, all = FALSE, 
                  d = 4, alpha = .05, orderBy = NULL, decreasing = TRUE){
  mll <- function(object){
    k <- length(object$fitobj)
    res <- as.matrix(do.call(cbind, lapply(object$fitobj, resid)))
    n <- nrow(res)
    sigma <- (t(res) %*% res)/n
    inv <- solve(sigma)
    ll <- sum(((-k/2) * log(2 * pi)) - (.5 * log(det(sigma))) - (.5 * diag(res %*% inv %*% t(res))))
    df <- sum(sapply(lapply(object$mods, '[[', "model"), nrow)) + 1
    return(c(LL = ll, df = df))
  }
  uni_mll <- function(object){
    getInd <- function(ind, x){
      ind <- c("deviance", "LL_model", "df.residual", "AIC", "BIC")[ind]
      X <- ifelse(ind == "df.residual", list(x$fitobj), list(x$mods))[[1]]
      return(unname(sapply(X, '[[', ind)))
    }
    out <- do.call(cbind.data.frame, lapply(1:5, getInd, x = object))
    colnames(out) <- c("RSS", "LL", "df", "AIC", "BIC")
    rownames(out) <- names(object$mods)
    out
  }
  omni_mll <- function(object){
    k <- length(object$fitobj)
    n <- nrow(object$data) * k
    ll <- mll(object = object)
    aic <- (2 * ll[2]) - (2 * ll[1])
    bic <- (ll[2] * log(n)) - (2 * ll[1])
    out <- c(ll, AIC = unname(aic), BIC = unname(bic))
    out
  }
  mod_lrt <- function(object, d = 4, alpha = .05, N = NULL){
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
    if(is.logical(net1)){nodes <- net1}
    nn <- ifelse(!is.null(names(net0)), list(names(net0)), list(nn))[[1]]
    net1 <- net0[[2]]; net0 <- net0[[1]]
  }
  if(isTRUE(attr(net0, "mlGVAR"))){net0 <- net0$betweenNet}
  if("ggm" %in% names(attributes(net0))){
    omni0 <- switch(2 - isTRUE('modLL' %in% names(net0)), net0$modLL$omnibus, omni_mll(net0))
    uni0 <- switch(2 - isTRUE('modLL' %in% names(net0)), net0$modLL$nodes, uni_mll(net0))
    if(is.logical(net1)){nodes <- net1; net1 <- NULL}
    if(!is.null(net1)){
      if(isTRUE(attr(net1, "mlGVAR"))){net1 <- net1$betweenNet}
      if(is.null(lrt)){lrt <- TRUE}
      omni1 <- switch(2 - isTRUE('modLL' %in% names(net1)), net1$modLL$omnibus, omni_mll(net1))
      uni1 <- switch(2 - isTRUE('modLL' %in% names(net1)), net1$modLL$nodes, uni_mll(net1))
      omni0 <- rbind(omni0, omni1)
      uni0 <- list(uni0, uni1)
      rownames(omni0) <- names(uni0) <- nn
    }
  } else {
    if(all(sapply(net0, function(z) isTRUE(attr(z, "mlGVAR"))))){
      net0 <- lapply(net0, '[[', "betweenNet")
    }
    stopifnot(all(sapply(net0, function(z) isTRUE(attr(z, "ggm")))))
    yLL <- all(sapply(net0, function(z) 'modLL' %in% names(z)))
    if(is.null(lrt)){lrt <- FALSE}
    if(is.logical(net1)){nodes <- net1}
    omni0 <- t(switch(2 - yLL, sapply(lapply(net0, '[[', 'modLL'), '[[', 'omnibus'), sapply(net0, omni_mll)))
    if(!is.null(orderBy)){
      orderBy <- switch(match.arg(tolower(as.character(
        orderBy)), c("true", "ll", "loglik", "aic", "bic", "df", "rank")), 
        ll =, loglik = "LL", aic = "AIC", bic = "BIC", "df")
      net0 <- net0[order(omni0[, orderBy], decreasing = decreasing)]
      #net0 <- net0[order(sapply(net0, function(z){
      #  omni_mll(z)[orderBy]}), decreasing = decreasing)]
    }
    nn <- ifelse(!is.null(names(net0)), list(names(net0)), 
                 list(paste0("net", 0:(length(net0) - 1))))[[1]]
    uni0 <- switch(2 - yLL, lapply(lapply(net0, '[[', 'modLL'), '[[', 'nodes'), lapply(net0, uni_mll))
    #omni0 <- do.call(rbind, lapply(net0, omni_mll))
    #uni0 <- lapply(net0, uni_mll)
    rownames(omni0) <- names(uni0) <- nn
  }
  out <- ifelse(all, list(list(nodes = uni0, omnibus = omni0)), 
                ifelse(nodes, list(uni0), list(omni0)))[[1]]
  if(ifelse(!is.null(net1), lrt, FALSE)){
    lrt0 <- mod_lrt(object = omni0, d = d, alpha = alpha, N = prod(dim(net1$dat$Y)))
    lrt1 <- mod_lrt(object = uni0, d = d, alpha = alpha)
    out <- list(nodes = uni0, LRT = lrt1, omnibus = lrt0)
    if(!all){out <- ifelse(nodes, list(lrt1), list(lrt0))[[1]]}
  }
  return(out)
}

##### modTable: obtain all possible LRTs (with RMSEAs) comparing a list of models
modTable <- function(fits, nodes = FALSE, orderBy = TRUE, d = 4, alpha = .05, 
                     decreasing = TRUE, names = NULL, rmsea = FALSE){
  n <- length(fits)
  stopifnot(is.list(fits) & n > 2)
  if(!is.null(names)){names(fits) <- names}
  if(is.null(names(fits))){names(fits) <- paste0("fit", 1:n)}
  if(all(sapply(fits, function(z) isTRUE(attr(z, "mlGVAR"))))){
    fits <- lapply(fits, '[[', "betweenNet")
  }
  stopifnot(all(sapply(fits, function(z) isTRUE(attr(z, "ggm")))))
  if(!is.null(orderBy)){
    orderBy <- switch(match.arg(tolower(as.character(
      orderBy)), c("true", "ll", "loglik", "aic", "bic", "df", "rank")), 
      ll =, loglik = "LL", aic = "AIC", bic = "BIC", "df")
    fits <- fits[order(sapply(fits, function(z){
      modLL(z)[orderBy]}), decreasing = decreasing)]
  }
  tt <- combn(n, 2)
  lls0 <- modLL(fits)
  lrts0 <- lapply(seq_len(ncol(tt)), function(i){
    modLL(fits[tt[, i]], d = d, alpha = alpha)})
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
    lls1 <- modLL(fits, nodes = TRUE)
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
      modLL(fits[tt[, i]], nodes = TRUE, d = d, alpha = alpha)})
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

##### binomialPL: pseudo-likelihood for binomial network
binomialPL <- function(object, m = NULL){
  dat <- as.matrix(object$data)
  k <- ncol(dat)
  n <- nrow(dat)
  if(is.null(m)){
    LLs <- sapply(1:k, function(z){
      Y <- dat[, z]
      X <- cbind(1, dat[, -z])
      b0 <- object$mods[[z]]$model[[1]][1, 1]
      B <- c(b0, object$adjMat[z, -z])
      yhat <- (X %*% B)[, 1]
      LL <- sum(Y * yhat - log(1 + exp(yhat)))
      return(LL)
    })
    return(sum(LLs))
  }
}

##### covNet: create network for plotting that includes covariates as nodes
covNet <- function(object, mnet = TRUE, threshold = .05){
  if(threshold == TRUE){threshold <- .05}
  if(class(mnet) == "character"){
    object$mods0$covariates$Bcovs <- list(object$mods0$Bm)
    names(object$mods0$covariates$Bcovs) <- mnet
    data <- object$mods0$dat
    mnet <- FALSE
  } else {
    data <- if(mnet){object$mods0$dat} else {object$data}
    data <- data.frame(data, object$mods0$covariates$covs)
  }
  if(mnet){if(!"mnet" %in% names(object)){mnet <- FALSE}}
  cs <- length(object$mods0$covariates$Bcovs)
  bb <- if(mnet){object$mnet$adjMat} else {object$adjMat}
  cp <- ncol(bb)
  cadj <- cedges <- matrix(0, cp + cs, cp + cs)
  cd <- matrix(FALSE, cp + cs, cp + cs)
  cadj[1:cp, 1:cp] <- bb
  dimnames(cadj) <- dimnames(cedges) <- dimnames(cd) <- rep(list(colnames(data)), 2)
  if(mnet){cd[1:cp, 1:cp] <- object$mnet$d}
  for(i in 1:cs){
    np <- nrow(object$mods0$covariates$Bcovs[[i]])
    if(threshold != FALSE){
      cadj[cp + i, 1:np] <- ifelse(object$mods0$covariates$Bcovs[[i]][, 4] <= threshold, 
                                   object$mods0$covariates$Bcovs[[i]][, 1], 0)
    } else {
      cadj[cp + i, 1:np] <- object$mods0$covariates$Bcovs[[i]][, 1]
    }
  }
  p <- ifelse(mnet, cp - 1, cp)
  if("modEdges" %in% names(object)){
    if(mnet){cedges[1:cp, 1:cp] <- object$mnet$modEdges} else {cedges[1:cp, 1:cp] <- object$modEdges}
    cedges[-c(1:cp), 1:p] <- 1
  }
  cd[-c(1:cp), 1:p] <- TRUE
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    colnames(colMat) <- colnames(adjMat)
    rownames(colMat) <- rownames(adjMat)
    colMat
  }
  out <- list(adjMat = cadj, edgeColors = getEdgeColors(cadj), modEdges = cedges, d = cd, data = data)
  if(all(cedges == 0)){out$modEdges <- NULL}
  attributes(out)$mnet <- mnet
  attributes(out)$threshold <- threshold
  attributes(out)$covs <- ifelse(mnet, cs + 1, cs)
  out
}

##### show which variables were selected in each model
selected <- function(object, threshold = FALSE, mod = 1){
  ints <- NULL
  if(threshold != FALSE & !is.numeric(threshold)){threshold <- .05}
  if(isTRUE(attr(object, "mlGVAR"))){object <- object[[mod + 1]]}
  if(isTRUE(attr(object, "SURnet"))){
    object <- object$SURnet
    mods0 <- object$temporal$coefs[[ifelse(threshold == FALSE, 1, 2)]]
    mods <- lapply(seq_len(nrow(mods0)), function(z){
      if(threshold == FALSE){
        z <- colnames(mods0)[mods0[z, ] != 0]
      } else {
        z <- colnames(mods0)[mods0[z, ] <= threshold]
      }
      z <- replace(z, length(z) == 0, "")
      z <- replace(z, is.null(z), "")
      return(z[z != "(Intercept)"])
    })
    names(mods) <- rownames(mods0)
  } else {
    if(threshold != FALSE & "fitobj" %in% names(object)){
      mods <- lapply(object$fitobj, function(z){
        z <- summary(z)$coefficients[, 4]
        z <- ifelse(z <= threshold, 1, 0)
        z <- names(z)[which(z == 1)]
        z <- replace(z, length(z) == 0, "")
        z <- replace(z, is.null(z), "")
        return(z[z != "(Intercept)"])
      })
    } else {
      mods <- lapply(object$mods, function(z){
        z <- rownames(z$model)[-1]
        return(replace(z, length(z) == 0, ""))
      })
    }
  }
  if(any(grepl(":", mods))){
    ints <- lapply(mods, function(z){
      z <- z[grepl(":", z)]
      return(replace(z, length(z) == 0, ""))
    })
    ix <- max(sapply(ints, length))
    ints <- do.call(cbind.data.frame, lapply(ints, function(z){
      if(length(z) < ix){z <- c(z, rep("", ix - length(z)))}
      return(z)
    }))
    mods <- lapply(mods, function(z){
      z <- z[!grepl(":", z)]
      return(replace(z, length(z) == 0, ""))
    })
  }
  mx <- max(sapply(mods, length))
  mods <- do.call(cbind.data.frame, lapply(mods, function(z){
    if(length(z) < mx){z <- c(z, rep("", mx - length(z)))}
    return(z)
  }))
  out <- mods
  if(!is.null(ints)){out <- list(mods = mods, ints = ints)}
  return(out)
}

##### posteriorProbs
posteriorProbs <- function(bics, plot = FALSE, mods = NULL){
  if(is(bics, 'list')){
    if(all(sapply(bics, function(z) isTRUE(attr(z, 'ggm'))))){
      bics <- modLL(bics)[, 'BIC']
    } else if(all(sapply(bics, function(z) isTRUE(attr(z, 'SURnet'))))){
      bics <- SURll(bics)[, 'BIC']
    }
  } else if(is(bics, 'matrix')){
    if('BIC' %in% colnames(bics)){bics <- bics[, 'BIC']}
  }
  BICs <- Rmpfr::mpfr(bics, 100)
  BICtrans <- exp(-0.5 * BICs)
  sumBICtrans <- sum(BICtrans)
  modelProbability <- as.numeric(BICtrans/sumBICtrans)
  if(identical(plot, FALSE)){
    return(modelProbability)
  } else {
    if(is.null(mods) & is.null(names(bics))){
      mods <- paste0('Model', 1:length(bics))
    } else {
      mods <- switch(2 - !is.null(mods), mods, names(bics))
    }
    dat <- cbind.data.frame(model = mods, modelProbability = modelProbability)
    invisible(suppressMessages(require(ggplot2)))
    g <- ggplot(dat, aes(x = model, y = modelProbability)) + geom_bar(stat = 'identity') + 
      xlab('') + ylab('Posterior model probability') + theme_bw() + ylim(0, 1)
    return(g)
  }
}


### ======================================================================== ###
### ============================= BOOTSTRAPPING ============================ ###
### ======================================================================== ###
##### bootNet: essentially 'bootnet', but with fitNetwork for all estimation
bootNet <- function(data, m = NULL, nboots = 10, lags = NULL, caseDrop = FALSE, rule = 'OR',
                    ci = .95, caseMin = .05, caseMax = .75, caseN = 10, threshold = FALSE,
                    fits = NULL, type = 'g', saveMods = TRUE, verbose = TRUE, fitCoefs = FALSE, 
                    size = NULL, nCores = 1, cluster = 'mclapply', directedDiag = FALSE, ...){
  if(identical(m, 0)){m <- NULL}
  args <- tryCatch({list(...)}, error = function(e){list()})
  call <- as.list(match.call())
  if(is.numeric(threshold)){threshold <- ifelse(threshold == 0, 'fit0', ifelse(threshold == 1, 'fits', threshold))}
  sampThresh <- switch(2 - (threshold == 'fit0'), TRUE, ifelse(threshold == 'fits', FALSE, threshold))
  threshold <- ifelse(threshold == 'fits', TRUE, ifelse(threshold == 'fit0', FALSE, threshold))
  if(any(sapply(list(data, fits), function(z) isTRUE(attr(z, 'resample'))))){
    dat <- switch(2 - isTRUE(attr(data, 'resample')), data, replace(fits, 'data', list(data = data)))
    stopifnot('data' %in% names(dat))
    caseDrop <- FALSE
    inds <- do.call(cbind.data.frame, lapply(dat$samples$iters, '[[', 'samp_inds'))
    fits <- lapply(dat$samples$iters, '[[', 'fit')
    attr(fits, 'resample') <- TRUE
    lags <- ifelse('lags' %in% names(dat$call), dat$call$lags, FALSE)
    nboots <- dat$call$niter
    m <- call$m <- dat$call$moderators
    attributes(fits)[c('inds', 'lags', 'm')] <- list(inds, lags, m)
    if('rule' %in% names(dat$call)){rule <- dat$call$rule}
    args0 <- dat$call[intersect(names(dat$call), formalArgs(fitNetwork))]
    args0$threshold <- sampThresh
    data <- args0$data <- dat$data
    if(!'fit0' %in% names(args)){
      fit0 <- do.call(fitNetwork, append(args0, list(saveMods = FALSE)))
    } else {
      fit0 <- switch(2 - isTRUE(args$fit0), modSelect(dat, data, TRUE), args$fit0)
      args$fit0 <- NULL
    }
  }
  ci <- (1 - ci)/2
  n <- n0 <- nrow(data)
  lags <- switch(2 - !is.null(lags), ifelse(all(lags != 0), 1, 0), 0)
  if(is.null(fits)){
    if(!caseDrop){
      n <- ifelse(is.null(size), n, size)
      inds <- data.frame(replicate(nboots, sample(1:(n0 - lags), n, replace = TRUE)))
    } else {
      mm <- as.numeric(!is.null(m))
      p <- ncol(data) - mm
      if(lags){
        mm <- as.numeric(isTRUE(!is.null(m)))
        p <- ncol(data) - mm
        fixMax <- round((1 - caseMax) * n) < ((p + 1) * (mm + 2))
        if(fixMax){
          caseMax <- 1 - (((p + 1) * (mm + 2))/n)
          message(paste0('caseMax too large; setting to max value: ', round(caseMax, 2)))
        }
      }
      subCases <- round((1 - seq(caseMin, caseMax, length = caseN)) * n)
      subNs <- sample(subCases, nboots, replace = TRUE)
      inds <- lapply(subNs, function(z) sort(sample(seq_len(n - lags), z)))
    }
    if(verbose & identical(nCores, 1)){pb <- txtProgressBar(min = 0, max = nboots + 1, style = 3)}
    args0 <- list(data = data, moderators = m, type = type, lags = lags, rule = rule, 
                  threshold = sampThresh, verbose = FALSE, 
                  saveMods = FALSE, fitCoefs = fitCoefs)
    args <- args[setdiff(names(args), names(args0))]
    args0 <- append(args0, args[intersect(names(args), formalArgs(fitNetwork))])
    fit0 <- do.call(fitNetwork, args0)
    if('fit0' %in% names(args)){fit0 <- args$fit0}
    if(verbose & identical(nCores, 1)){setTxtProgressBar(pb, 1)}
    args0$threshold <- threshold
    if(nCores > 1 | isTRUE(nCores)){
      if(isTRUE(nCores)){nCores <- parallel::detectCores()}
      if(grepl('Windows', sessionInfo()$running)){cluster <- 'SOCK'}
      if(tolower(cluster) != 'mclapply'){
        cluster <- match.arg(toupper(cluster), c('SOCK', 'FORK'))
        cl <- parallel::makeCluster(nCores, type = cluster)
        if(cluster == 'SOCK'){
          obj1 <- switch(2 - !as.logical(lags), c('nodewise', 'modNet', 'modLL'), 
                         c('lagMat', 'SURfit', 'SURnet', 'SURll', 'surEqs', 'getCoefs', 'systemfit'))
          obj1 <- c(obj1, 'lags', 'inds', 'data', 'args0', 'm', 'caseDrop')
          parallel::clusterExport(cl, c('fitNetwork', 'Matrix', 'net', 'netInts', obj1), envir = environment())
        }
      } else {
        cl <- nCores
      }
      if(verbose){
        pbapply::pboptions(type = 'timer', char = '-')
        fits <- suppressWarnings(structure(pbapply::pblapply(seq_len(nboots), function(z){
          newargs <- args0
          newargs$data <- switch(lags + 1, data[inds[[z]], ], structure(data, samp_ind = inds[[z]]))
          fit <- tryCatch({do.call(fitNetwork, newargs)}, error = function(e){
            stop(paste0("Sampling error", ifelse(caseDrop, ": Try reducing 'caseMax'", "")))})
          return(fit)
        }, cl = cl), class = c('list', 'bootFits'), inds = inds, lags = lags, m = m))
      } else if(tolower(cluster) == 'mclapply'){
        fits <- suppressWarnings(structure(parallel::mclapply(seq_len(nboots), function(z){
          args0$data <- switch(lags + 1, data[inds[[z]], ], structure(data, samp_ind = inds[[z]]))
          fit <- tryCatch({do.call(fitNetwork, args0)}, error = function(e){
            stop(paste0("Sampling error", ifelse(caseDrop, ": Try reducing 'caseMax'", "")))})
          return(fit)
        }, mc.cores = nCores), class = c('list', 'bootFits'), inds = inds, lags = lags, m = m))
      } else {
        fits <- suppressWarnings(structure(parallel::parLapply(cl, seq_len(nboots), function(z){
          args0$data <- switch(lags + 1, data[inds[[z]], ], structure(data, samp_ind = inds[[z]]))
          fit <- tryCatch({do.call(fitNetwork, args0)}, error = function(e){
            stop(paste0("Sampling error", ifelse(caseDrop, ": Try reducing 'caseMax'", "")))})
          return(fit)
        }), class = c('list', 'bootFits'), inds = inds, lags = lags, m = m))
      }
      if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
      rm(cl)
    } else {
      fits <- structure(lapply(seq_len(nboots), function(z){
        args0$data <- switch(lags + 1, data[inds[[z]], ], structure(data, samp_ind = inds[[z]]))
        fit <- tryCatch({do.call(fitNetwork, args0)}, error = function(e){
          stop(paste0("Sampling error", ifelse(caseDrop, ": Try reducing 'caseMax'", "")))})
        if(verbose){setTxtProgressBar(pb, z + 1)}
        return(fit)
      }), class = c('list', 'bootFits'), inds = inds, lags = lags, m = m)
    }
    if(any(sapply(fits, class) == 'try-error')){
      err <- which(sapply(fits, class) == 'try-error')
      inds <- inds[-err]
      nboots <- nboots - length(err)
      fits <- structure(fits[-err], class = c('list', 'bootFits'), inds = inds, lags = lags, m = m)
      if(verbose){message(paste0(length(err), ' iterations failed'))}
    }
    net2 <- net3 <- FALSE
  } else {
    if(is(fits, 'bootNet')){
      if(!'fit0' %in% names(args)){fit0 <- fits$fit0}
      if('bootFits' %in% names(fits)){fits <- fits$bootFits} else {stop('No fits')}
    }
    if('lags' %in% names(attributes(fits))){lags <- attr(fits, 'lags')}
    if('inds' %in% names(attributes(fits))){inds <- attr(fits, 'inds')}
    if('m' %in% names(attributes(fits))){m <- attr(fits, 'm')}
    net2 <- ifelse('net2' %in% names(args), args$net2, FALSE)
    net3 <- ifelse('net3' %in% names(args), args$net3, FALSE)
    nboots <- length(fits)
    if(length(args) > 0){
      if('fit0' %in% names(args)){fit0 <- args$fit0}
      if('inds' %in% names(args)){inds <- args$inds}
      if('subNs' %in% names(args)){subNs <- args$subNs}
    }
  }
  if(directedDiag){
    strFUN <- switch(2 - net3, function(x){rowSums(x)}, function(x){colSums(x)})
  } else {
    strFUN <- switch(2 - net3, function(x){diag(x) <- 0; rowSums(x)}, function(x){diag(x) <- 0; colSums(x)})
  }
  nodes <- switch(lags + 1, names(fit0$mods), gsub('[.]y$', '', names(fit0$SURnet$mods)))
  v <- length(nodes)
  node01 <- node1 <- rep(nodes[-v], (v - 1):1)
  node02 <- node2 <- matrix(nodes, v, v)[lower.tri(matrix(nodes, v, v))]
  if(lags & net2){
    node1 <- rep(nodes, each = v)
    node2 <- rep(nodes, v)
  } else if(lags){
    m <- NULL
  }
  edges <- paste0(node1, ifelse(lags & net2, '-->', '--'), node2)
  e <- length(edges)
  edge1 <- rep(edges[-e], (e - 1):1)
  edge2 <- matrix(edges, e, e)[lower.tri(matrix(edges, e, e))]
  adj0 <- switch(factorial(lags + net2), 
                 net(fit0, 'PCC', sampThresh, rule)[lower.tri(net(fit0, 'PCC', sampThresh, rule))], 
                 c(net(fit0, threshold = sampThresh, rule = rule)))
  which.net <- ifelse(lags & net2, 'beta', 'PCC')
  str0 <- strFUN(abs(net(fit0, which.net, sampThresh, rule)))
  ei0 <- strFUN(net(fit0, which.net, sampThresh, rule))
  adj1 <- do.call(cbind, lapply(fits, function(z){
    switch(factorial(lags + net2), 
           net(z, 'PCC', threshold, rule)[lower.tri(net(z, 'PCC', threshold, rule))], 
           c(net(z, threshold = threshold, rule = rule)))
  }))
  colnames(adj1) <- names(inds) <- paste0('boot', 1:nboots)
  if(!caseDrop){
    adj2 <- as.data.frame(t(apply(adj1, 1, function(z) quantile(z, c(ci, 1 - ci)))))
    colnames(adj2) <- c('boot_lower', 'boot_upper')
    adjmeans <- data.frame(boot_mean = rowMeans(adj1), boot_sd = apply(adj1, 1, sd))
    adj2 <- data.frame(edge = edges, node1 = node1, node2 = node2, adjmeans, adj2,
                       sample = adj0, sample_lower = adj0 + (adjmeans$boot_sd * qnorm(ci)),
                       sample_upper = adj0 + (adjmeans$boot_sd * qnorm(1 - ci)))
  } else {
    uniqueNs <- sort(unique(subNs))
    uN <- length(uniqueNs)
    if(any(adj0 == 0)){
      groupEdges <- lapply(split(data.frame(t(adj1))[, -which(adj0 == 0)], subNs), 
                           function(z){apply(z, 1, function(zz) cor(zz, adj0[adj0 != 0]))})
    } else {
      groupEdges <- lapply(split(data.frame(t(adj1)), subNs), 
                           function(z){apply(z, 1, function(zz) cor(zz, adj0))})
    }
    adjmeans <- lapply(groupEdges, function(z) c(mean = mean(z), sd = sd(z)))
    qs <- c(.01, .025, .05, .25, .5, .75, .95, .975, .99)
    adjqs <- data.frame(t(sapply(groupEdges, quantile, probs = qs)))
    colnames(adjqs) <- paste0('q', c(1, 2.5, 5, 25, 50, 75, 95, 97.5, 99))
    adjmeans <- data.frame(do.call(rbind, adjmeans))
    adj2 <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2), 
                       subN = uniqueNs, n = sapply(groupEdges, length), 
                       adjmeans, adjqs)
    adj2 <- adj2[order(adj2$subN, decreasing = FALSE), ]
    rownames(adj2) <- 1:nrow(adj2)
  }
  rownames(adj1) <- edges
  if(!is.null(m)){
    mname <- ifelse(length(fit0$call$moderators) == 1, fit0$call$moderators, 'MOD')
    inodes <- paste0(nodes, ':', mname); iedges <- paste0(edges, '|', mname)
    inode1 <- paste0(node1, ':', mname); inode2 <- paste0(node2, ':', mname)
    iedge1 <- paste0(edge1, '|', mname); iedge2 <- paste0(edge2, '|', mname)
    if(isTRUE(attr(fit0, 'ggm'))){
      ints0 <- netInts(fit0, threshold = threshold, rule = rule, avg = TRUE)
      istr0 <- colSums(abs(ints0))
      iei0 <- colSums(ints0)
      ints0 <- ints0[lower.tri(ints0)]
      ints1 <- do.call(cbind, lapply(1:nboots, function(z){
        i1 <- netInts(fits[[z]], threshold = threshold, rule = rule, avg = TRUE)
        if(length(i1) == 0){i1 <- diag(0, v)}
        return(i1[lower.tri(i1)])
      }))
      #ints0 <- fit0$interactions$intMats$avgInts[lower.tri(fit0$interactions$intMats$avgInts)]
      #istr0 <- colSums(abs(fit0$interactions$intMats$avgInts))
      #iei0 <- colSums(fit0$interactions$intMats$avgInts)
      #ints1 <- do.call(cbind, lapply(1:nboots, function(z){
      #  fits[[z]]$interactions$intMats$avgInts[lower.tri(fits[[z]]$interactions$intMats$avgInts)]
      #}))
    } else {
      ints0 <- c(netInts(fit0, threshold = sampThresh))
      istr0 <- strFUN(abs(netInts(fit0, threshold = sampThresh)))
      iei0 <- strFUN(netInts(fit0, threshold = sampThresh))
      ints1 <- do.call(cbind, lapply(1:nboots, function(z){
        z1 <- netInts(fits[[z]], threshold = threshold)
        if(length(z1) == 0){z1 <- diag(0, v)}
        return(c(z1))
      }))
    }
    colnames(ints1) <- paste0('boot', 1:nboots)
    if(!caseDrop){
      ints2 <- as.data.frame(t(apply(ints1, 1, function(z) quantile(z, c(ci, 1 - ci)))))
      colnames(ints2) <- c('boot_lower', 'boot_upper')
      intsmeans <- data.frame(boot_mean = rowMeans(ints1), boot_sd = apply(ints1, 1, sd))
      ints2 <- data.frame(edge = iedges, node1 = inode1, node2 = inode2, 
                          intsmeans, ints2, sample = ints0, 
                          sample_lower = ints0 + (intsmeans$boot_sd * qnorm(ci)),
                          sample_upper = ints0 + (intsmeans$boot_sd * qnorm(1 - ci)))
    } else {
      if(any(ints0 == 0)){
        groupInts <- lapply(split(data.frame(t(ints1))[, -which(ints0 == 0)], subNs), 
                            function(z){apply(z, 1, function(zz) cor(zz, ints0[ints0 != 0]))})
      } else {
        groupInts <- lapply(split(data.frame(t(ints1)), subNs), 
                            function(z){apply(z, 1, function(zz) cor(zz, ints0))})
      }
      intsmeans <- lapply(groupInts, function(z) c(mean = mean(z), sd = sd(z)))
      intsqs <- data.frame(t(sapply(groupInts, quantile, probs = qs)))
      colnames(intsqs) <- paste0('q', c(1, 2.5, 5, 25, 50, 75, 95, 97.5, 99))
      intsmeans <- data.frame(do.call(rbind, intsmeans))
      ints2 <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2), 
                          subN = uniqueNs, n = sapply(groupInts, length), 
                          intsmeans, intsqs)
      ints2 <- ints2[order(ints2$subN, decreasing = FALSE), ]
      rownames(ints2) <- 1:nrow(ints2)
    }
    rownames(ints1) <- iedges
  }
  strengths <- intStr <- EIs <- intEIs <- list()
  if(lags & net2){
    fnets <- lapply(fits, net, which.net, threshold, rule)
    strengths <- data.frame(do.call(rbind, lapply(fnets, function(z) strFUN(abs(z)))))
    EIs <- data.frame(do.call(rbind, lapply(fnets, function(z) strFUN(z))))
    if(!is.null(m)){
      fnets2 <- lapply(fits, netInts, threshold = threshold)
      intStr <- data.frame(do.call(rbind, lapply(fnets2, function(z) strFUN(abs(z)))))
      intEIs <- data.frame(do.call(rbind, lapply(fnets2, function(z) strFUN(z))))
    }
    node1 <- node01
    node2 <- node02
  } else {
    for(i in 1:length(nodes)){
      strengths[[i]] <- adj1[which(node1 == nodes[i] | node2 == nodes[i]), ]
      EIs[[i]] <- colSums(strengths[[i]])
      strengths[[i]] <- colSums(abs(strengths[[i]]))
      if(!is.null(m)){
        intStr[[i]] <- ints1[which(node1 == nodes[i] | node2 == nodes[i]), ]
        intEIs[[i]] <- colSums(intStr[[i]])
        intStr[[i]] <- colSums(abs(intStr[[i]]))
      }
    }
    strengths <- do.call(cbind.data.frame, strengths)
    EIs <- do.call(cbind.data.frame, EIs)
  }
  colnames(strengths) <- colnames(EIs) <- nodes
  if(!is.null(m)){
    if(!(lags & net2)){
      intStr <- do.call(cbind.data.frame, intStr)
      intEIs <- do.call(cbind.data.frame, intEIs)
    }
    colnames(intStr) <- colnames(intEIs) <- inodes
  }
  if(!caseDrop){
    estDiffs <- function(ests, ci = .025){
      comps <- t(combn(ncol(ests), 2))
      diffOut <- as.data.frame(t(sapply(1:nrow(comps), function(z){
        ints <- quantile(ests[, comps[z, 1]] - ests[, comps[z, 2]], c(ci, 1 - ci))
        return(c(ints, ifelse(ints[1] <= 0 & ints[2] >= 0, FALSE, TRUE)))
      })))
      colnames(diffOut) <- c('lower', 'upper', 'sig')
      diffOut$sig <- as.logical(diffOut$sig)
      return(diffOut)
    }
    sdiffs <- data.frame(node1 = node1, node2 = node2, estDiffs(strengths, ci))
    eidiffs <- data.frame(node1 = node1, node2 = node2, estDiffs(EIs, ci))
    smeans <- data.frame(node = nodes, boot_mean = colMeans(strengths), boot_sd = apply(strengths, 2, sd))
    squants <- as.data.frame(t(apply(strengths, 2, function(z) quantile(z, c(ci, 1 - ci)))))
    eimeans <- data.frame(node = nodes, boot_mean = colMeans(EIs), boot_sd = apply(EIs, 2, sd))
    eiquants <- as.data.frame(t(apply(EIs, 2, function(z) quantile(z, c(ci, 1 - ci)))))
    colnames(squants) <- colnames(eiquants) <- c('boot_lower', 'boot_upper')
    smeans <- data.frame(smeans, squants, sample = str0, sample_lower = str0 + (smeans$boot_sd * qnorm(ci)),
                         sample_upper = str0 + (smeans$boot_sd * qnorm(1 - ci)))
    eimeans <- data.frame(eimeans, eiquants, sample = ei0, sample_lower = ei0 + (eimeans$boot_sd * qnorm(ci)),
                          sample_upper = ei0 + (eimeans$boot_sd * qnorm(1 - ci)))
    rownames(smeans) <- rownames(eimeans) <- 1:v
    edge_diffs <- data.frame(edge1 = edge1, edge2 = edge2, estDiffs(t(adj1), ci))
  } else {
    pairStr <- pairEI <- multiStr <- multiEI <- vector('list', length = uN)
    psmeans <- pemeans <- msmeans <- memeans <- vector('list', length = uN)
    for(i in 1:uN){
      colN <- which(subNs == uniqueNs[i])
      pairStr[[i]] <- apply(t(strengths)[, colN, drop = FALSE], 2, function(z) cor(z, str0))
      pairEI[[i]] <- apply(t(EIs)[, colN, drop = FALSE], 2, function(z) cor(z, ei0))
      psmeans[[i]] <- c(mean = mean(pairStr[[i]]), sd = sd(pairStr[[i]]))
      pemeans[[i]] <- c(mean = mean(pairEI[[i]]), sd = sd(pairEI[[i]]))
      if(i == uN){
        psmeans <- data.frame(do.call(rbind, psmeans))
        pemeans <- data.frame(do.call(rbind, pemeans))
        psqs <- data.frame(t(sapply(pairStr, quantile, probs = qs)))
        peqs <- data.frame(t(sapply(pairEI, quantile, probs = qs)))
        colnames(psqs) <- colnames(peqs) <- paste0('q', c(1, 2.5, 5, 25, 50, 75, 95, 97.5, 99))
        pairStrengths <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2),
                                    subN = uniqueNs, n = sapply(pairStr, length),
                                    psmeans, psqs)
        pairExInf <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2),
                                subN = uniqueNs, n = sapply(pairEI, length),
                                pemeans, peqs)
        pairStrengths <- pairStrengths[order(pairStrengths$subN, decreasing = FALSE), ]
        pairExInf <- pairExInf[order(pairExInf$subN, decreasing = FALSE), ]
        rownames(pairExInf) <- rownames(pairStrengths) <- 1:nrow(pairStrengths)
      }
      if(!is.null(m)){
        multiStr[[i]] <- apply(t(intStr)[, colN, drop = FALSE], 2, function(z) cor(z, istr0))
        multiEI[[i]] <- apply(t(intEIs)[, colN, drop = FALSE], 2, function(z) cor(z, iei0))
        msmeans[[i]] <- c(mean = mean(multiStr[[i]]), sd = sd(multiStr[[i]]))
        memeans[[i]] <- c(mean = mean(multiEI[[i]]), sd = sd(multiEI[[i]]))
        if(i == uN){
          msmeans <- data.frame(do.call(rbind, msmeans))
          memeans <- data.frame(do.call(rbind, memeans))
          msqs <- data.frame(t(sapply(multiStr, quantile, probs = qs)))
          meqs <- data.frame(t(sapply(multiEI, quantile, probs = qs)))
          colnames(msqs) <- colnames(meqs) <- paste0('q', c(1, 2.5, 5, 25, 50, 75, 95, 97.5, 99))
          intStrengths <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2),
                                     subN = uniqueNs, n = sapply(multiStr, length),
                                     msmeans, msqs)
          intExInf <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2),
                                 subN = uniqueNs, n = sapply(multiEI, length),
                                 memeans, meqs)
          intStrengths <- intStrengths[order(intStrengths$subN, decreasing = FALSE), ]
          intExInf <- intExInf[order(intExInf$subN, decreasing = FALSE), ]
          rownames(intExInf) <- rownames(intStrengths) <- 1:nrow(intStrengths)
        }
      }
    }
  }
  if(!is.null(m) & !caseDrop){
    inode1 <- paste0(node1, ':', mname); inode2 <- paste0(node2, ':', mname)
    iedge1 <- paste0(edge1, '|', mname); iedge2 <- paste0(edge2, '|', mname)
    int_squants <- as.data.frame(t(apply(intStr, 2, function(z) quantile(z, c(ci, 1 - ci)))))
    int_eiquants <- as.data.frame(t(apply(intEIs, 2, function(z) quantile(z, c(ci, 1 - ci)))))
    colnames(int_squants) <- colnames(int_eiquants) <- c('boot_lower', 'boot_upper')
    int_edge_diffs <- data.frame(edge1 = iedge1, edge2 = iedge2, estDiffs(t(ints1), ci))
    int_sdiffs <- data.frame(node1 = inode1, node2 = inode2, estDiffs(intStr, ci))
    int_eidiffs <- data.frame(node1 = inode1, node2 = inode2, estDiffs(intEIs, ci))
    int_smeans <- data.frame(node = inodes, boot_mean = colMeans(intStr), boot_sd = apply(intStr, 2, sd))
    int_smeans <- data.frame(int_smeans, int_squants, sample = istr0, 
                             sample_lower = istr0 + (int_smeans$boot_sd * qnorm(ci)),
                             sample_upper = istr0 + (int_smeans$boot_sd * qnorm(1 - ci)))
    int_eimeans <- data.frame(node = inodes, boot_mean = colMeans(intEIs), boot_sd = apply(intEIs, 2, sd))
    int_eimeans <- data.frame(int_eimeans, int_eiquants, sample = iei0, 
                              sample_lower = iei0 + (int_eimeans$boot_sd * qnorm(ci)),
                              sample_upper = iei0 + (int_eimeans$boot_sd * qnorm(1 - ci)))
    rownames(int_smeans) <- rownames(int_eimeans) <- 1:v
    out <- list(pairwise = list(
      edges = adj2, strength = smeans, EI = eimeans, diffs = list(
        edge_diffs = edge_diffs, str_diffs = sdiffs, EI_diffs = eidiffs)),
      interactions = list(edges = ints2, strength = int_smeans, EI = int_eimeans, diffs = list(
        int_edge_diffs = int_edge_diffs, int_str_diffs = int_sdiffs, int_EI_diffs = int_eidiffs)),
      boots = list(boot_edges = adj1, boot_strengths = strengths, boot_EIs = EIs, boot_int_edges = ints1, 
                   boot_int_strengths = intStr, boot_int_EI = intEIs, boot_inds = inds))
  } else if(is.null(m) & !caseDrop){
    out <- list(edges = adj2, strength = smeans, EI = eimeans,
                diffs = list(edge_diffs = edge_diffs, str_diffs = sdiffs, EI_diffs = eidiffs),
                boots = list(boot_edges = adj1, boot_strengths = strengths, 
                             boot_EIs = EIs, boot_inds = inds))
  } else {
    pairwise <- list(edges = adj2, strength = pairStrengths, EI = pairExInf)
    boots <- list(boot_edges = adj1, boot_strengths = strengths, boot_EIs = EIs)
    out <- append(pairwise, list(boots = boots))
    if(!is.null(m)){
      interactions <- list(edges = ints2, strength = intStrengths, EI = intExInf)
      boots2 <- list(boot_int_edges = ints1, boot_int_strengths = intStr, boot_int_EI = intEIs)
      out <- list(pairwise = pairwise, interactions = interactions, boots = append(boots, boots2))
    }
    out$boots$boot_inds <- inds
    out$subN <- subNs
  }
  if(lags & !net2){
    if(isTRUE(attr(fits, 'resample'))){
      attr(fits, 'resample') <- NULL
      if(!'nboots' %in% names(call)){call$nboots <- nboots}
      caseDrop <- FALSE
    }
    args1 <- list(data = data, fit0 = fit0, inds = inds, net2 = TRUE,
                  fits = fits, lags = lags, saveMods = FALSE)
    args2 <- append(args1, args[setdiff(names(args), names(args1))])
    args2 <- append(args2, call[-1][setdiff(names(call[-1]), names(args2))])
    if(caseDrop){args2$subNs <- subNs} else {args2$caseDrop <- FALSE}
    out2 <- tryCatch({do.call(match.fun(call[[1]]), args2)}, error = function(e){list()})
    out3 <- tryCatch({do.call(match.fun(call[[1]]), append(args2, list(net3 = TRUE)))}, error = function(e){list()})
    attributes(out)[c('ci', 'caseDrop', 'class')] <- attributes(out2)[c('ci', 'caseDrop', 'class')]
    attr(out, 'net') <- 'contemporaneous'; attr(out2, 'net') <- 'temporal'
    makeTemporal <- function(out2, out3, caseDrop, subNs = NULL){
      makeTemp <- function(out2, out3, caseDrop, subNs = NULL){
        temp1 <- append(list(edges = out2$edges), list(
          strength = list(outStrength = out2$strength, inStrength = out3$strength), 
          EI = list(outEI = out2$EI, inEI = out3$EI)
        ))
        boot1 <- append(list(boot_edges = out2$boots$boot_edges), list(
          boot_strengths = list(boot_outStrengths = out2$boots$boot_strengths, boot_inStrengths = out3$boots$boot_strengths),
          boot_EIs = list(boot_outEIs = out2$boots$boot_EIs, boot_inEIs = out3$boots$boot_EIs),
          boot_inds = out2$boots$boot_inds
        ))
        if(!caseDrop){
          diffs1 <- append(list(edge_diffs = out2$diffs$edge_diffs), list(
            str_diffs = list(outStr_diffs = out2$diffs$str_diffs, inStr_diffs = out3$diffs$str_diffs),
            EI_diffs = list(outEI_diffs = out2$diffs$EI_diffs, inEI_diffs = out3$diffs$EI_diffs)
          ))
          out2 <- append(temp1, list(diffs = diffs1, boots = boot1))
        } else {
          out2 <- append(temp1, list(boots = boot1, subN = subNs))
        }
        return(out2)
      }
      if(!'interactions' %in% names(out2)){
        out <- makeTemp(out2, out3, caseDrop, subNs)
      } else {
        bootNames <- names(out2$boots)
        pairs2 <- append(out2$pairwise, list(boots = out2$boots[!grepl('int', names(out2$boots))]))
        pairs3 <- append(out3$pairwise, list(boots = out3$boots[!grepl('int', names(out3$boots))]))
        ints2 <- append(out2$interactions, list(boots = out2$boots[grepl('int|inds', names(out2$boots))]))
        ints3 <- append(out3$interactions, list(boots = out3$boots[grepl('int|inds', names(out3$boots))]))
        names(ints2$boots) <- names(ints3$boots) <- names(pairs2$boots)
        names(ints2$diffs) <- names(ints3$diffs) <- names(pairs2$diffs)
        pairwise <- makeTemp(pairs2, pairs3, caseDrop, subNs)
        interactions <- makeTemp(ints2, ints3, caseDrop, subNs)
        if(caseDrop){pairwise$subN <- interactions$subN <- NULL}
        out <- list(pairwise = pairwise, interactions = interactions)
      }
      return(out)
    }
    atts1 <- attributes(out2)[-1]
    out2 <- makeTemporal(out2, out3, caseDrop, subNs)
    attributes(out2)[names(atts1)] <- atts1
    out <- list(temporal = out2, contemporaneous = out, fit0 = fit0)
    if(caseDrop){out$subN <- subNs}
  }
  names(fits) <- paste0('boot', 1:nboots)
  if(!lags){out$fit0 <- fit0}
  if(saveMods){out$bootFits <- fits}
  attr(out, 'n') <- n
  attr(out, 'ci') <- 1 - ci * 2
  attr(out, 'caseDrop') <- caseDrop
  class(out) <- c('list', 'bootNet')
  return(out)
}

##### plotBoot: plot results of bootNet; currently only available for caseDrop
plotBoot <- function(obj, type = 'edges', net = 'temporal', plot = 'all', cor = .7,
                     order = 'mean', ci = .95, pairwise = TRUE, interactions = TRUE, 
                     labels = NULL, title = NULL, cis = 'quantile', true = NULL,
                     errbars = FALSE, vline = FALSE, threshold = FALSE,
                     difference = FALSE, color = FALSE, text = FALSE,
                     textPos = 'value', multi = NULL, directedDiag = FALSE, ...){
  if(!is.null(multi)){
    call <- replace(as.list(match.call())[-1], 'multi', NULL)
    if(isTRUE(multi)){multi <- list(plot = c('pairwise', 'interactions'))}
    n <- names(multi)
    p1 <- invisible(do.call(plotBoot, replace(call, n, setNames(list(multi[[1]][1]), n))))
    legend <- switch(2 - ('legend' %in% names(call)), eval(call$legend), g_legend(p1))
    p2 <- invisible(do.call(plotBoot, replace(call, n, setNames(list(multi[[1]][2]), n))))
    return(gridExtra::grid.arrange(legend, gridExtra::arrangeGrob(
      p1 + theme(legend.position = 'none'), 
      p2 + theme(legend.position = 'none'), nrow = 1), 
      nrow = 2, heights = c(1, 10)))
  }
  args <- tryCatch({list(...)}, error = function(e){list()})
  fit0 <- switch(2 - ('fit0' %in% names(obj)), obj$fit0, list())
  if(isTRUE(attr(obj, 'resample'))){
    if(!'data' %in% names(obj)){return(plotCoefs(fit = obj, ...))}
    if(length(fit0) == 0){
      fit0 <- switch(2 - ('fit0' %in% names(args)), args$fit0, TRUE)
    }
    obj <- do.call(bootNet, append(list(data = obj, fit0 = fit0, directedDiag = directedDiag), 
                                   replace(args, 'fit0', NULL)))
  }
  if(!identical(threshold, FALSE) & 'bootFits' %in% names(obj)){
    dat <- paste0('obj$fit0$', ifelse('temporal' %in% names(obj), 'SURnet$data', 'data'))
    if(isTRUE(attr(obj$bootFits, 'resample'))){attributes(obj$bootFits)$resample <- NULL}
    obj <- bootNet(data = eval(parse(text = dat)), fits = obj$bootFits, 
                   threshold = threshold, fit0 = obj$fit0, 
                   directedDiag = directedDiag)
  }
  runonce <- TRUE
  plots <- c('none', 'pairwise', 'interactions', 'all', 'both')
  nets <- c('ggm', 'temporal', 'contemporaneous', 'between')
  types <- c('edges', 'strength', 'outstrength', 'instrength', 'ei', 'outei', 'inei')
  if(length(type) > 1){
    type <- match.arg(tolower(type), types, several.ok = TRUE)[1:2]
    stopifnot(all(type %in% c('outstrength', 'instrength')) | all(type %in% c('outei', 'inei')))
    runonce <- type[2]
    type <- type[1]
  }
  nty <- list(type, net, plot)
  plot <- which(!sapply(nty, is.character))
  plot <- ifelse(length(plot) == 0, unlist(nty)[nty %in% plots], nty[[plot]])
  nty <- unlist(nty)[-match(plot, unlist(nty))]
  nty <- match.arg(tolower(nty), c(nets, types), several.ok = TRUE)
  type <- nty[nty %in% types][1]
  net <- nty[nty %in% nets][1]
  net <- switch(net, between = 'ggm', 'NA' = 'temporal', net)
  type <- switch(type, outstrength = 'outStrength', instrength = 'inStrength', 
                 outei = 'outEI', inei = 'inEI', ei = 'EI', 'NA' = 'edges', type)
  cis <- match.arg(tolower(cis), c('quantile', 'se'))
  invisible(suppressMessages(require(ggplot2)))
  if(isTRUE(attr(obj, 'mlGVAR')) & 'varMods' %in% names(obj)){
    mlnet <- switch(net, ggm = 'between|means', 'fixed')
    if(is.null(title) & startsWith(mlnet, 'b')){title <- 'Between-subjects network\t'}
    fits <- obj$varMods[[grep(mlnet, names(obj$varMods))]]
    fit0 <- switch(2 - ('fit0' %in% names(args)), args$fit0, obj[[grep(mlnet, names(obj))]])
    data <- obj$netData[[grep(mlnet, names(obj$netData))]]
    args <- append(list(data = data, fits = fits, fit0 = fit0), replace(args, 'fit0', NULL))
    obj <- do.call(bootNet, args)
  }
  lags <- any(c('temporal', 'contemporaneous') %in% names(obj))
  if(lags){
    net <- switch(net, ggm = 'temporal', net)
    obj <- obj[[net]]
    if(net == 'contemporaneous'){boots <- obj$boots}
    title <- switch(2 - is.null(title), paste(Hmisc::capitalize(net), 'Network\t'), title)
    if(identical(title, FALSE)){title <- NULL}
    if(type == 'strength' & net == 'temporal'){
      type <- 'outStrength'
    } else if(!type %in% c('edges', 'strength', 'EI') & net == 'contemporaneous'){
      type <- 'edges'
    }
    if(type == 'EI' & net == 'temporal'){
      type <- 'outEI'
    }
  } else {
    if(!type %in% c('edges', 'strength', 'EI')){type <- 'edges'}
    net <- 'ggm'
    boots <- obj$boots
  }
  lags2 <- lags & (net == 'temporal')
  if(is.logical(plot)){plot <- ifelse(plot, 'all', 'none')}
  if(is.numeric(plot)){plot <- c('none', 'pairwise', 'interactions', 'all')[plot + 1]}
  pp <- match.arg(tolower(plot), c('all', 'both', 'pairwise', 'interactions', 'none'))
  pp <- switch(pp, both = 'all', none = FALSE, pp)
  if(pp == 'all' & difference | lags & !lags2 & difference){pp <- 'pairwise'}
  if(pp == 'interactions'){interactions <- TRUE; pairwise <- FALSE}
  if(pp == 'pairwise'){interactions <- FALSE; pairwise <- TRUE}
  type0 <- ifelse(type %in% c('outStrength', 'inStrength'), paste0('strength$', type), type)
  type0 <- ifelse(type %in% c('outEI', 'inEI'), paste0('EI$', type), type0)
  if(!"interactions" %in% names(obj)){interactions <- FALSE}
  obj1 <- list()
  if(interactions | !pairwise){
    obj1 <- obj$interactions
    if(lags){boots <- obj1$boots}
    p1 <- dat <- eval(parse(text = paste0('obj$interactions$', type0)))
    if(pairwise){
      text <- FALSE
      obj1 <- obj$pairwise
      obj2 <- obj$interactions
      if('diffs' %in% names(obj1) & 'diffs' %in% names(obj2)){
        obj1$diffs <- setNames(lapply(seq_along(obj1$diffs), function(z){
          rbind(obj1$diffs[[z]], obj2$diffs[[z]])
        }), names(obj1$diffs))
      }
      p2 <- eval(parse(text = paste0('obj$pairwise$', type0)))
      dat <- rbind.data.frame(p2, p1)
      dat$type <- rep(c("Pairwise", "Interactions"), each = nrow(dat)/2)
    } else {
      dat <- data.frame(dat)
      dat$type <- rep("Interactions", nrow(dat))
    }
  } else {
    boots <- obj$boots
    if("pairwise" %in% names(obj)){
      obj1 <- obj$pairwise
      if(lags){boots <- obj1$boots}
      type0 <- paste0('pairwise$', type0)
    }
    dat <- eval(parse(text = paste0('obj$', type0)))
    dat$type <- rep("Pairwise", nrow(dat))
  }
  if('diffs' %in% names(obj1)){
    obj1 <- obj1$diffs
  } else if(lags & difference){
    obj1 <- obj$diffs
  }
  if(!identical(text, FALSE) & type == 'edges'){ # had !lags
    boots <- boots[[switch(2 - interactions, ifelse(
      !lags, 'boot_int_edges', 'boot_edges'), 'boot_edges')]]
  }
  if(grepl('EI', type)){type <- gsub('I', 'Influence', gsub('E', 'Expected', type))}
  if(isTRUE(attr(obj, "caseDrop"))){
    dat <- dat2 <- dat[order(dat$subN), ]
    ci <- paste0("q", c((1 - ci)/2, 1 - (1 - ci)/2) * 100)
    if(length(unique(dat$type)) == 1 & FALSE){ # MAY DELETE THIS WHOLE PART
      css <- dat[dat[, ci[1]] > cor, 'drop']
      css <- attr(dat2, 'CS') <- ifelse(length(css) == 0, 0, max(css))
      if(pp != FALSE){cat(paste0('CS: ', css, ' (cor = ', cor, ')'))}
    } else if(FALSE){ # MAY DELETE
      cssPair <- dat[dat$type == 'Pairwise', 'drop'][dat[dat$type == 'Pairwise', ci[1]] > cor]
      cssPair <- attr(dat2, 'CS_Pair') <- ifelse(length(cssPair) == 0, 0, max(cssPair))
      cssInt <- dat[dat$type == 'Interactions', 'drop'][dat[dat$type == 'Interactions', ci[1]] > cor]
      cssInt <- attr(dat2, 'CS_Int') <- ifelse(length(cssInt) == 0, 0, max(cssInt))
      if(pp != FALSE){cat(paste0('CS_Pair: ', cssPair, ' (cor = ', cor, ')'), '\n')}
      if(pp != FALSE){cat(paste0('CS_Int: ', cssInt, ' (cor = ', cor, ')'))}
    }
    legLab <- Hmisc::capitalize(type)
    N <- dat$N[1]
    p <- ggplot(dat, aes(x = subN, y = mean, group = type, colour = type, fill = type))
    p <- p + geom_line(lwd = 1) + geom_point()
    p <- p + geom_ribbon(colour = NA, alpha = .1, aes_string(ymin = ci[1], ymax = ci[2]))
    p <- p + scale_x_reverse(breaks = seq(.9, .1, by = -.1) * N, 
                             labels = paste0(seq(90, 10, by = -10), "%"))
    p <- p + ylim(-1, 1) + xlab("Sampled cases") + ylab("Average correlation with original sample")
    p <- p + labs(fill = legLab, colour = legLab, group = legLab) + theme_bw()
    p <- p + geom_hline(yintercept = 0) + geom_hline(yintercept = cor, linetype = 2, alpha = .3)
    if(!is.null(title)){p <- p + ggtitle(title)}
  } else if(difference){
    if(order %in% c(6, 'true')){order <- TRUE}
    if(is.logical(order)){order <- ifelse(order, 'mean', 'id')}
    if(is.numeric(order)){order <- switch(order, 'id', 'sample', 'mean', 'v1', 'v2', 'true')}
    order <- match.arg(order, c('id', 'sample', 'mean', 'v1', 'v2', 'true'))
    if(!order %in% c('id', 'sample', 'mean')){order <- 'id'}
    order <- switch(order, mean = 'boot_mean', id = ifelse(type == 'edges', 'edge', 'node'), order)
    dat0 <- dat
    if(all(startsWith(names(obj1), 'int_'))){names(obj1) <- gsub('^int_', '', names(obj1))}
    if(startsWith(type, 'out') | startsWith(type, 'in')){
      tt <- ifelse(gsub('out|in', '', type) == 'Strength', 's', 'E')
      dat <- obj1[[which(sapply(names(obj1), startsWith, tt))]]
      dat <- dat[[which(sapply(c('^out', '^in'), grepl, type))]]
    } else {
      dat <- obj1[[which(sapply(names(obj1), startsWith, strsplit(type, '')[[1]][1]))]]
    }
    dat <- rbind(dat, cbind(setNames(dat[, 2:1], names(dat)[1:2]), dat[, -(1:2)]))
    dat$sig <- ifelse(dat$sig, 'sig', 'nonsig')
    d2 <- setNames(data.frame(unique(dat[, 1]), unique(dat[, 1]), 0, 0, 'same'), names(dat))
    dat <- dat2 <- rbind(dat, d2)
    colorValues <- c(same = 'white', nonsig = 'lightgray', sig = 'black')
    colnames(dat)[1:2] <- paste0('id', 1:2)
    id0 <- ifelse(type == 'edges', 'edge', 'node')
    dat$id1 <- factor(dat$id1, levels = dat0[[id0]][order(dat0[[order]])])
    dat$id2 <- factor(dat$id2, levels = dat0[[id0]][order(dat0[[order]])])
    dat$type <- factor(paste0(ifelse(pp == 'pairwise', 'Pairwise (', 'Interactions ('), Hmisc::capitalize(type), ')'))
    if(color & type == 'edges'){
      if(net == 'ggm'){net <- 'temporal'}
      nn <- switch(2 - pairwise, net(fit0, n = net), netInts(fit0))
      if(lags2 & !pairwise){nn <- t(nn)}
      nn <- rownames(nn)
      obj0 <- switch(2 - pairwise, fit0, t(netInts(fit0, threshold = threshold, avg = !lags2)))
      graph <- plotNet(obj0, which.net = net, threshold = threshold, plot = FALSE)
      edgelist <- data.frame(from = nn[graph$Edgelist$from], to = nn[graph$Edgelist$to], 
                             col = graph$graphAttributes$Edges$color, stringsAsFactors = FALSE)
      if(all(grepl(':', edgelist$from))){
        n1 <- gsub(':.*', '', edgelist$from)
        n2 <- gsub(':.*', '', edgelist$to)
        mm <- unique(gsub('.*:', '', edgelist$from))
        edgelist$id <- paste0(paste0(n1, ifelse(lags2, '-->', '--'), n2), '|', mm)
      } else {
        if(lags){
          edgelist$from <- gsub('[.]y$', '', edgelist$from)
          edgelist$to <- gsub('[.]y$', '', edgelist$to)
        }
        edgelist$id <- paste0(edgelist$from, ifelse(lags2, '-->', '--'), edgelist$to)
      }
      if(any(is.na(edgelist$col))){edgelist$col[is.na(edgelist$col)] <- 'white'}
      dat$sig <- ifelse(dat$sig == 'same', as.character(dat$id1), dat$sig)
      colorValues <- c(setNames(edgelist$col, edgelist$id), colorValues)
    }
    p <- ggplot(dat, aes(x = id1, y = id2, fill = sig)) + 
      geom_tile(colour = 'white') + xlab('') + ylab('') + 
      scale_fill_manual(values = colorValues) + theme(legend.position = 'none') + facet_grid(~ type)
    p <- p + theme_grey(base_size = 9) + 
      theme(legend.position = 'none', axis.text.x = element_text(
        size = 7.2, angle = 270, hjust = 0, colour = 'grey50'))
    if(text){
      n3 <- as.character(dat[!dat$sig %in% c('sig', 'nonsig'), 'id1'])
      dat3 <- data.frame(id1 = as.character(dat[!dat$sig %in% c('sig', 'nonsig'), 'id1']), value = NA)
      dat3[match(dat0[[id0]], dat3$id1), 'value'] <- format(signif(dat0$sample, 2), scientific = FALSE)
      dat3$id1 <- dat3$id2 <- factor(dat3$id1)
      dat3$sig <- dat[!dat$sig %in% c('sig', 'nonsig'), 'sig']
      p <- p + geom_text(data = dat3, aes(label = value))
    }
  } else {
    nets <- c('mean', 'sample')
    if(!is.null(true)){
      getType <- switch(type, edges = switch(net, temporal = c, function(x){x[lower.tri(x)]}), 
                        strength = , outStrength = function(x){diag(x) <- 0; colSums(abs(x))}, 
                        inStrength = function(x){diag(x) <- 0; rowSums(abs(x))}, ExpectedInfluence = , 
                        outExpectedInfluence = function(x){diag(x) <- 0; colSums(x)}, 
                        inExpectedInfluence = function(x){diag(x) <- 0; rowSums(x)})
      dat$true <- numeric(nrow(dat))
      dat[dat$type == 'Pairwise', 'true'] <- getType(net(true, n = switch(net, ggm = 'between', net)))
      if('Interactions' %in% dat$type){dat[dat$type == 'Interactions', 'true'] <- getType(netInts(true))}
      nets <- c(nets, 'true')
    }
    if(is.null(true) & order %in% c('true', 6)){order <- TRUE}
    if(is.logical(order)){order <- ifelse(order, 'mean', 'id')}
    if(is.numeric(order)){order <- switch(order, 'id', 'sample', 'mean', 'v1', 'v2', 'true')}
    order <- match.arg(order, c('id', 'sample', 'mean', 'v1', 'v2', 'true'))
    order <- which(order == c('id', 'sample', 'mean', 'v1', 'v2', 'true'))
    id <- ifelse(type == 'edges', 'edge', 'node')
    lci <- ifelse(cis == 'quantile', 'boot_lower', 'sample_lower')
    uci <- ifelse(cis == 'quantile', 'boot_upper', 'sample_upper')
    if(ci != .95 & cis == 'quantile' & type == 'edges'){
      dat$boot_lower <- unname(apply(obj$boots$boot_edges, 1, quantile, probs = (1 - ci)/2))
      dat$boot_upper <- unname(apply(obj$boots$boot_edges, 1, quantile, probs = 1 - (1 - ci)/2))
    }
    v2 <- all(unique(dat$type) == c('Pairwise', 'Interactions'))
    if(!v2 & order > 3 & order != 6){order <- 1}
    if(order > 3 & order != 6){
      dat <- dat[order(dat$type, decreasing = TRUE), ]
      ord <- c('sample', 'boot_mean')[order - 3]
      o1 <- order(dat[dat$type == 'Pairwise', ord])
      o2 <- order(dat[dat$type == 'Interactions', ord]) + max(o1)
      dat <- dat[c(o1, o2), ]
      order <- 1
    } else if(order == 6){
      order <- 4
    }
    v2 <- TRUE
    type2 <- switch(2 - v2, paste0(rep(dat$type, length(nets)), ' (', Hmisc::capitalize(type), ')'), type)
    dat2 <- cbind.data.frame(
      id = rep(dat[[id]], length(nets)), value = unname(unlist(dat[, gsub('mean', 'boot_mean', nets)])), 
      lower = rep(dat[[lci]], length(nets)), upper = rep(dat[[uci]], length(nets)), 
      net = rep(nets, each = nrow(dat)), type = factor(Hmisc::capitalize(type2))
    )
    colnames(dat2)[1] <- id
    dat2[[id]] <- factor(dat2[[id]], levels = dat[[id]][switch(
      order, 1:nrow(dat), order(dat$sample), order(dat$boot_mean), order(dat$true))])
    colnames(dat2)[1] <- 'id'
    if(is.null(labels)){labels <- isTRUE(nrow(dat) <= 50)}
    scaleSize <- c('mean' = .7, 'sample' = 1, 'true' = .9)[seq_along(nets)]
    scaleAlpha <- c('mean' = .5, 'sample' = 1, 'true' = .5)[seq_along(nets)]
    legLabs <- c('Bootstrap mean', 'Sample', 'True network')[seq_along(nets)]
    legCols <- c('black', 'darkred', 'blue')[seq_along(nets)]
    if(!isTRUE(runonce) & !(pairwise & interactions)){
      call <- replace(as.list(match.call())[-1], c('type', 'pairwise', 'interactions', 'plot'), 
                      list(type = runonce, pairwise = pairwise, interactions = interactions, plot = FALSE))
      dat2 <- rbind.data.frame(dat2, do.call(plotBoot, call))
    }
    p <- ggplot(dat2, aes(x = id, y = value, group = net, colour = net, fill = net))
    p <- p + geom_point(aes(alpha = net), show.legend = c(alpha = FALSE))
    if(!errbars){
      p <- p + geom_line(aes(size = net, alpha = net), show.legend = c(colour = FALSE, alpha = FALSE, size = FALSE))
      p <- p + geom_ribbon(aes(ymin = lower, ymax = upper, colour = NULL, fill = NULL), alpha = .1)
    } else {
      p <- p + geom_errorbar(aes(ymin = lower, ymax = upper, colour = NULL, fill = NULL), width = .2)
    }
    p <- p + facet_grid(~ type) + coord_flip() + theme_bw() + xlab('') + ylab('')
    p <- p + guides(colour = guide_legend(title = title), fill = 'none') + theme(legend.position = "top")
    p <- p + scale_size_manual(values = scaleSize) + scale_alpha_manual(values = scaleAlpha)
    p <- p + scale_color_manual('', values = legCols, labels = legLabs)
    if(!labels){p <- p + theme(axis.text.y = element_blank())}
    if(!identical(vline, FALSE)){
      p <- p + geom_hline(yintercept = ifelse(!is.numeric(vline), 0, vline), linetype = 2, alpha = .35)
    }
    if(!identical(text, FALSE) & type == 'edges'){
      if(isTRUE(text)){text <- 1.2}
      rnames <- rownames(boots)
      dat2$prop0 <- rep(NA, nrow(dat2))
      dat2 <- dat2[dat2$net == 'mean', ]
      dat2[match(rnames, as.character(dat2$id)), 'prop0'] <- rowMeans(data.frame((boots != 0) * 1))
      if(!identical(pp, FALSE) & textPos != 'value'){dat2$value <- textPos}
      p <- p + geom_label(aes(y = value, label = format(round(prop0, 2), nsmall = 2), fill = NULL), 
                          data = dat2, cex = text * 2, label.padding = unit(0.1, 'lines'), 
                          label.size = 0.1, alpha = .8, colour = 'black')
    }
  }
  if(pp != FALSE){return(p)} else {return(dat2)}
}

##### cscoef: get correlation stability coefficients from caseDrop procedure
cscoef <- function(obj, cor = .7, ci = .95, first = FALSE, verbose = TRUE){
  stopifnot(isTRUE(attr(obj, 'caseDrop')))
  ci0 <- ci * 100
  ci <- paste0("q", c((1 - ci)/2, 1 - (1 - ci)/2) * 100)[1]
  lags <- any(c('temporal', 'contemporaneous') %in% names(obj))
  inds <- c('edges', switch(2 - !lags, c('strength', 'EI'), c('outStrength', 'inStrength', 'outEI', 'inEI')))
  if(lags){
    obj0 <- obj
    obj <- obj$temporal
  }
  if('interactions' %in% names(obj)){
    inds2 <- c('pairwise', 'interactions')
    out <- do.call(rbind, setNames(lapply(inds2, function(z){
      obj2 <- obj[[z]]
      setNames(sapply(inds, function(i){
        if(grepl('^out|^in', i)){
          ii <- tolower(gsub('^out|^in', '', i))
          if(ii == 'ei'){ii <- 'EI'}
          dat <- obj2[[ii]][[i]]
          dat <- dat[order(dat$subN), ]
        } else {
          dat <- obj2[[i]][order(obj2[[i]]$subN), ]
        }
        if(!first){
          css <- dat[dat[, ci] > cor, 'drop']
          css <- ifelse(length(css) == 0, 0, max(css))
        } else {
          css <- which(!dat[, ci] > cor)
          css <- dat[ifelse(length(css) == 0, 1, max(css) + 1), 'drop']
        }
        return(css)
      }), inds)
    }), inds2))
  } else {
    out <- setNames(sapply(inds, function(i){
      if(grepl('^out|^in', i)){
        ii <- tolower(gsub('^out|^in', '', i))
        if(ii == 'ei'){ii <- 'EI'}
        dat <- obj[[ii]][[i]]
        dat <- dat[order(dat$subN), ]
      } else {
        dat <- obj[[i]][order(obj[[i]]$subN), ]
      }
      if(!first){
        css <- dat[dat[, ci] > cor, 'drop']
        css <- ifelse(length(css) == 0, 0, max(css))
      } else {
        css <- which(!dat[, ci] > cor)
        css <- dat[ifelse(length(css) == 0, 1, max(css) + 1), 'drop']
      }
      return(css)
    }), inds)
  }
  if(lags){
    call <- replace(as.list(match.call())[-1], c('obj', 'verbose'), 
                    list(obj = obj0$contemporaneous, verbose = FALSE))
    out <- list(temporal = out, contemporaneous = do.call(cscoef, call))
  }
  if(verbose){cat(paste0('CS = ', cor, '(', ci0, '%)'), '\n')}
  return(out)
}


##### splitNets: conduct a sample-split test using either NCT or FGL
splitNets <- function(data, niter = 100, fun = 'NCT', gamma = .5, rule = 'and', 
                      m = 'M', alpha = .05, true = NULL, mval = 'median', 
                      alloutput = TRUE, verbose = TRUE, ...){
  args0 <- tryCatch({list(...)}, error = function(e){list()})
  rule <- match.arg(tolower(rule), c('and', 'or'))
  FUN <- fun <- match.arg(tolower(fun), c('nct', 'estimategroupnetwork', 'groupnetwork', 'group', 'fgl'))
  FUN <- switch(FUN, nct = NetworkComparisonTest::NCT, EstimateGroupNetwork::EstimateGroupNetwork)
  if(is(data, 'list')){
    if('data' %in% names(data)){
      true <- netInts(data)
      data <- data$data
    }
  }
  m0 <- ifelse(is.character(m), which(colnames(data) == m), m)
  m <- data[, m]
  if(dim(table(m)) == 2){
    mm <- unique(m)
    dat1 <- data[m == mm[1], -m0]
    dat2 <- data[m == mm[2], -m0]
  } else if(mval == 'median'){
    dat1 <- data[m < median(m), -m0]
    dat2 <- data[m >= median(m), -m0]
  } else {
    dat1 <- data[m < mval, -m0]
    dat2 <- data[m > mval, -m0]
  }
  p <- ncol(dat1)
  edges <- t(combn(p, 2))
  edges <- split(edges, 1:nrow(edges))
  args <- list(data1 = dat1, data2 = dat2)
  if(fun == 'nct'){
    args <- append(args, list(it = niter, test.edges = TRUE, edges = edges, 
                              AND = (rule == 'and'), gamma = gamma, 
                              verbose = verbose))
    if(!verbose){args$progressbar <- FALSE}
  } else {
    args <- append(list(X = args, gamma = gamma), args0)
    args$nlambda1 <- args$nlambda2 <- niter
    args$simplifyOutput <- !alloutput
  }
  if(length(args0) > 0){args <- append(args, args0[setdiff(names(args0), names(args))])}
  out <- out0 <- suppressMessages(suppressWarnings(do.call(FUN, args)))
  ints <- diag(0, p)
  if(fun == 'nct'){
    if(any(as.numeric(out$einv.pvals[, 3]) <= alpha)){
      ps <- which(as.numeric(out$einv.pvals[, 3]) <= alpha)
      for(i in seq_along(ps)){ints[edges[[ps[i]]][1], edges[[ps[i]]][2]] <- 1}
      ints <- Matrix::forceSymmetric(ints)
    }
    out <- list(ints = as.matrix(ints), NCT = out)
    if(!is.null(true)){
      out$performance <- performance(out$ints, true, inds = 'between')
    }
  } else {
    if(alloutput){out <- out0$network}
    diffs <- out$data2 - out$data1
    ints[out$data1 != out$data2] <- diffs[out$data1 != out$data2]
    if(!is.null(true)){
      ints <- (sign(true) + sign(ints) != 0) * 1
      out$performance <- performance(ints, sign(true), inds = 'between')
    }
    out$ints <- ints
    if(alloutput){out$output <- out0[setdiff(names(out0), 'network')]}
  }
  class(out) <- c('list', 'splitNets')
  return(out)
}


### ======================================================================== ###
### =========================== PLOTTING FUNCTIONS ========================= ###
### ======================================================================== ###
##### condPlot: plot conditional effects
condPlot <- function(out, to, from, swap = FALSE, avg = FALSE, compare = NULL, 
                     hist = FALSE, xlab = NULL, mods = NULL, nsims = 500, 
                     xn = NULL, getCIs = FALSE, discrete = FALSE, 
                     ylab = NULL, main = NULL, midline = TRUE){
  require(ggplot2)
  if("adjMat" %in% names(out)){out <- out$mods0}
  if(any(c("models", "SURnet") %in% names(out))){
    out <- condEffects(out, xn = xn, x = compare)}
  if(is.null(xlab)){xlab <- attributes(out)$moderator}
  if(length(to) == 2){from <- to[2]; to <- to[1]}
  if(class(to) == "numeric"){to <- names(out)[to]}
  if(class(from) == "numeric"){from <- names(out)[from]}
  if(swap){
    toX <- to
    to <- from
    from <- toX
  }
  data <- out[[to]][[from]]
  stopifnot(!is.null(data))
  stopifnot(nrow(data) > 1)
  if(avg & !"SURnet" %in% names(attributes(out))){
    stopifnot(!is.null(out[[from]][[to]]))
    newdat <- data.frame(matrix(NA, ncol = 5, nrow = nrow(data)))
    newdat[, 1] <- data$x
    newdat[, 2] <- rowMeans(cbind(data$y, out[[from]][[to]]$y))
    newdat[, 3] <- rowMeans(cbind(data$se, out[[from]][[to]]$se))
    newdat[, 4] <- rowMeans(cbind(data$lower, out[[from]][[to]]$lower))
    newdat[, 5] <- rowMeans(cbind(data$upper, out[[from]][[to]]$upper))
    colnames(newdat) <- c("x", "y", "se", "lower", "upper")
    data <- newdat
  }
  to2 <- Hmisc::capitalize(to)
  fr2 <- Hmisc::capitalize(from)
  xlab <- Hmisc::capitalize(xlab)
  if(!avg){
    if(is.null(ylab)){ylab <- paste0(fr2, " ---> ", to2)}
    if(is.null(main)){mainlab <- paste0("Estimated Effect of ", fr2, " on ", to2, " across levels of ", xlab)}
  } else if(!"SURnet" %in% names(attributes(out))){
    if(is.null(ylab)){ylab <- paste0(fr2, " ~ ", to2)}
    if(is.null(main)){mainlab <- paste0("Mean Relation between ", fr2, " and ", to2, " across levels of ", xlab)}
  }
  if(!is.null(main)){mainlab <- main}
  alpha <- attributes(out)$alpha
  if("mods" %in% names(attributes(out))){mods <- attributes(out)$mods}
  if(!is.null(mods) & !"SURnet" %in% names(attributes(out))){
    ci_diff <- margCIs(mods = mods, alpha = alpha, nsims = nsims, compare = compare)
    if(avg){ci_diff2 <- ci_diff[[from]][to, ]}
    ci_diff <- ci_diff[[to]][from, ]
  } else {
    ci_diff <- NULL
  }
  if(nrow(data) != 2 & discrete == FALSE){
    if(!hist){
      pp <- ggplot(data = data, aes(x = x, y = y)) + 
        geom_line(color = "red") +
        geom_ribbon(alpha = .2, aes(ymin = lower, ymax = upper))
      if(min(data$lower) < 0 & max(data$upper) > 0){pp <- pp + geom_hline(yintercept = 0, linetype = 2)}
      pp <- pp + xlab(xlab) + ylab(ylab) + ggtitle(mainlab) + theme_bw()
    } else {
      var2_dt <- ifelse(
        "SURnet" %in% names(attributes(out)), 
        list(mods), mods$dat[, ncol(mods$dat)])[[1]]
      yrange <- c(data$upper, data$lower)
      maxdiff <- (max(yrange) - min(yrange))
      breaks_var2 <- nrow(data)
      if(hist == "unique"){breaks_var2 <- length(unique(var2_dt))}
      hist.out <- hist(var2_dt, breaks = breaks_var2, plot = F)
      n.hist <- length(hist.out$mids)
      dist <- hist.out$mids[2] - hist.out$mids[1]
      hist.max <- max(hist.out$counts)
      histX <- data.frame(ymin = rep(min(yrange) - maxdiff/5, n.hist), 
                          ymax = hist.out$counts/hist.max * maxdiff/5 + min(yrange) - maxdiff/5,
                          xmin = hist.out$mids - dist/2, xmax = hist.out$mids + dist/2)
      pp <- ggplot()
      pp <- pp + geom_rect(data = histX, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                           colour = "gray50", alpha = 0, size = .5)
      pp <- pp + geom_line(data = data, aes(x = x, y = y), color = "red")
      pp <- pp + geom_ribbon(data = data, aes(x = x, ymin = lower, ymax = upper), alpha = .2)
      if(min(data$lower) < 0 & max(data$upper) > 0){pp <- pp + geom_hline(yintercept = 0, linetype = 2)}
      pp <- pp + xlab(xlab) + ylab(ylab) + ggtitle(mainlab) + theme_bw()
    }
  } else {
    pp <- ggplot(data = data, aes(x = x, y = y)) + 
      geom_point() + geom_errorbar(aes(ymin = lower, ymax = upper), width = .05)
    if(midline){pp <- pp + geom_hline(yintercept = mean(data$y), linetype = 3, alpha = .6, color = "red")}
    if(min(data$lower) < 0 & max(data$upper) > 0){pp <- pp + geom_hline(yintercept = 0, linetype = 2, alpha = .6)}
    pp <- pp + scale_x_discrete(limits = data$x)
    pp <- pp + xlab(xlab) + ylab(ylab) + ggtitle(mainlab) + theme_bw()
  }
  if(!is.null(ci_diff)){
    if(!is.null(compare)){
      compare <- round(compare, 2)
      citxt <- paste0(100 - (alpha * 100), "% CI(", compare[2], " - ", compare[1], "): [")
    } else {
      citxt <- paste0(100 - (alpha * 100), "% CI(Max - Min): [")
    }
    if(avg){
      a1 <- paste0(fr2, " ---> ", to2, ":  ")
      b1 <- paste0(to2, " ---> ", fr2, ":  ")
      a2 <- paste0(round(ci_diff[2], 3), ", ", round(ci_diff[3], 3), "]")
      b2 <- paste0(round(ci_diff2[2], 3), ", ", round(ci_diff2[3], 3), "]")
      pp <- pp + labs(subtitle = paste0(a1, citxt, a2, "\n", b1, citxt, b2))
      if(getCIs){
        ci_diff0 <- list(ci_diff, ci_diff2)
        names(ci_diff0) <- c(paste0(fr2, ":", to2), paste0(to2, ":", fr2))
        return(ci_diff0)
      }
    } else {
      pp <- pp + labs(subtitle = paste0(citxt, round(ci_diff[2], 3), ", ", round(ci_diff[3], 3), "]"))
    }
  }
  pp
}

##### intsPlot: plot CI(max - min) for each predictor with interactions against Y
intsPlot <- function(out, y = 'all', moderator = NULL, nsims = 500, alpha = .05){
  require(ggplot2)
  if("adjMat" %in% names(out)){out <- out$mods0}
  if("models" %in% names(out)){out <- margCIs(out, modname = moderator, nsims = nsims, alpha = alpha)}
  if(is.null(moderator)){moderator <- attributes(out)$moderator}
  if(class(y) == "character"){
    if(y == "all"){
      y0 <- as.vector(sapply(seq_along(out), function(z) 
        paste0(names(out)[z], ":", rownames(out[[z]]))))
      dd <- suppressWarnings(data.frame(do.call(rbind, out)))
      if(any(grepl(":$", y0))){y0 <- y0[-grep(":$", y0)]}
      if(class(y0) == "list"){y0 <- unlist(y0)}
      rownames(dd) <- y0
    } else if(y == "sig"){
      out2 <- lapply(seq_along(out), function(z){
        z2 <- data.frame(out[[z]], sig = as.numeric(!(out[[z]][, 2] < 0 & out[[z]][, 3] > 0)), 
                         check.names = FALSE)
        return(z2[z2$sig == 1, -4])
      })
      names(out2) <- names(out)
      if(any(sapply(out2, length) == 0)){
        out2 <- out2[-which(sapply(out2, length) == 0)]
      }
      y0 <- unlist(as.vector(sapply(seq_along(out2), function(z)
        paste0(names(out2)[z], ":", rownames(out2[[z]])))))
      dd <- suppressWarnings(data.frame(do.call(rbind, out2)))
      if(any(grepl(":$", y0))){y0 <- y0[-grep(":$", y0)]}
      if(class(y0) == "list"){y0 <- unlist(y0)}
      rownames(dd) <- y0
    } else {
      y <- which(names(out) == y)
      dd <- data.frame(out[[y]])
      Y <- Hmisc::capitalize(names(out)[y])
    }
  } else {
    dd <- data.frame(out[[y]])
    Y <- Hmisc::capitalize(names(out)[y])
  }
  M <- Hmisc::capitalize(moderator)
  if(y %in% c("sig", "all")){
    Y <- c("Significant interaction effects", "All interactions")[which(c("sig", "all") %in% y)]
    main <- paste0(Y, " across levels of ", M)
  } else {
    main <- paste0("Effects on ", Y, " across levels of ", M)
  }
  colnames(dd)[2:3] <- c("lower", "upper")
  dd$pred <- paste0("_", rownames(dd))
  p <- ggplot(data = dd, aes(x = reorder(pred, b), y = b)) +
    geom_point() + geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) + 
    geom_hline(yintercept = 0, linetype = 2) + 
    scale_x_discrete(labels = function(x){sub("[^_]*_", "", x)}) + 
    ggtitle(label = main) + xlab(Y) + ylab(expression(hat(beta))) + 
    theme_bw() + coord_flip()
  p
}

##### plotMods: plot networks at different levels of moderator
plotMods <- function(nets, nodewise = FALSE, elsize = 2, vsize = NULL, 
                     elabs = TRUE, predict = NULL, layout = NULL,
                     which.net = "temporal", ...){
  if(any(c("SURnet", "temporal") %in% unlist(lapply(nets, names)))){
    if("SURnet" %in% names(nets[[1]])){nets <- lapply(nets, '[[', "SURnet")}
    which.net <- match.arg(tolower(
      which.net), c("temporal", "contemporaneous", "pdc"))
    if(is.null(vsize)){vsize <- 10}
    nodewise <- FALSE
    ggm <- FALSE
  } else {
    ggm <- TRUE
  }
  if(is.null(vsize)){
    vsize <- (exp(-ncol(nets[[1]]$data)/80) * 8) + 1
    vsize <- vsize + (exp(-length(nets)/80) * 8)/length(nets)
  }
  getLayout <- function(x, which.net = "temporal"){
    stopifnot(class(x) == "list")
    averageLayout(lapply(x, function(z){
      plotNet(z, plot = FALSE, which.net = which.net)}))
  }
  getMax <- function(x, n = FALSE, which.net = "temporal"){
    if(!"ggm" %in% names(attributes(x[[1]]))){
      which.net <- match.arg(
        tolower(which.net), c("temporal", "contemporaneous", "pdc"))
      x <- lapply(x, '[[', ifelse(which.net == "pdc", "temporal", which.net))
      if(which.net == "pdc"){x <- lapply(x, '[[', "PDC")}
    }
    if(n == FALSE){
      max(unlist(lapply(x, function(z) max(abs(z$adjMat)))))
    } else {
      max(unlist(lapply(x, function(z) max(abs(z$nodewise$adjNW)))))
    }
  }
  if(is.null(layout)){layout <- getLayout(nets, which.net)}
  mx <- getMax(nets, nodewise, which.net)
  if(ggm){
    moderator <- Hmisc::capitalize(attr(nets[[1]], "moderator"))
    vals <- unname(sapply(nets, attr, "mval"))
  } else {
    moderator <- Hmisc::capitalize(nets[[1]]$call$moderators)
    vals <- unname(sapply(lapply(nets, '[[', "call"), '[[', "mval"))
  }
  layout(t(1:length(vals)))
  for(i in 1:length(vals)){
    plotNet(nets[[i]], predict = predict, layout = layout, elabs = elabs,
            elsize = elsize, nodewise = nodewise, maximum = mx, 
            title = paste0(moderator, " = ", vals[i]), 
            vsize = vsize, which.net = which.net, ...)
  }
}

##### plotCovNet: plot network using output from 'covNet'
plotCovNet <- function(object, border = NULL, color = NULL, ...){
  p <- attr(object, "cov")
  px <- ncol(object$data) - p
  shape <- c(rep("circle", px), rep("square", p))
  if(!is.null(border)){border <- c(rep(1, px), rep(border, p))}
  if(!is.null(color)){color <- c(rep("white", px), rep(color, p))}
  plotNet(object = object, directed = object$d, shape = shape, 
          border.width = border, color = color, ...)
}


### ======================================================================== ###
### =========================== HELPER FUNCTIONS =========================== ###
### ======================================================================== ###
##### margCIs: retrieve CIs for the effect of Z on the coefficient relating X to Y
margCIs <- function(mods, data = NULL, modname = NULL, alpha = .05, nsims = 500, 
                    compare = NULL, tTests = FALSE){
  set.seed(666)
  if(alpha == FALSE){alpha <- .05}
  if("adjMat" %in% names(mods)){mods <- mods$mods0}
  if(is.null(modname) & "dat" %in% names(mods)){modname <- colnames(mods$dat)[ncol(mods$dat)]}
  if(is.null(data) & "dat" %in% names(mods)){data <- mods$dat[, -which(colnames(mods$dat) == modname)]}
  if("models" %in% names(mods)){mods <- mods$models}
  if(!any(grepl(":", unlist(sapply(mods, function(z) names(coef(z))))))){stop("No interaction terms in the models")}
  p <- length(mods)
  vs <- colnames(data)
  if(is.null(compare)){
    xmin <- min(mods[[1]]$model[modname])
    xmax <- max(mods[[1]]$model[modname])
  } else if(length(compare) > 1){
    xmin <- compare[1]
    xmax <- compare[2]
  }
  msims <- min_sims <- max_sims <- ci_diff <- list()
  for(i in 1:p){
    msims[[i]] <- arm::sim(mods[[i]], nsims)
    intvars <- colnames(msims[[i]]@coef)[grep(":", colnames(msims[[i]]@coef))]
    vars <- which(colnames(msims[[i]]@coef) %in% vs)
    vars <- vars[colnames(msims[[i]]@coef)[vars] %in% gsub(":.*", "", intvars)]
    min_sims[[i]] <- max_sims[[i]] <- ci_diff[[i]] <- vector("list", length = length(vars))
    if(length(intvars) != 0){
      for(j in 1:length(vars)){
        intTerm <- paste0(colnames(msims[[i]]@coef)[vars[j]], ":", modname)
        min_sims[[i]][[j]] <- msims[[i]]@coef[, vars[j]] + xmin * msims[[i]]@coef[, intTerm]
        max_sims[[i]][[j]] <- msims[[i]]@coef[, vars[j]] + xmax * msims[[i]]@coef[, intTerm]
        dm <- max_sims[[i]][[j]] - min_sims[[i]][[j]]
        if(tTests){
          se <- sqrt(2 * (((var(min_sims[[i]][[j]]) + var(max_sims[[i]][[j]]))/2)/nsims))
          tval <- mean(dm)/se
          ci_diff[[i]][[j]] <- c(quantile(dm, alpha/2), quantile(dm, 1 - alpha/2), b = mean(dm), se = se,
                                 t = tval, p = 2 * pt(abs(tval), df = mods[[i]]$df.residual, lower.tail = F))
        } else {
          ci_diff[[i]][[j]] <- c(b = mean(dm), quantile(dm, alpha/2), quantile(dm, 1 - alpha/2))
        }
      }
      names(min_sims[[i]]) <- names(max_sims[[i]]) <- names(ci_diff[[i]]) <- colnames(msims[[i]]@coef)[vars]
    }
  }
  ci_diff <- lapply(ci_diff, function(z) do.call(rbind, z))
  names(ci_diff) <- vs
  attributes(ci_diff)$moderator <- modname
  if(!is.null(compare)){attributes(ci_diff)$compare <- compare[1:2]}
  ci_diff
}

##### getInts: retrieve interactions effects; e.g., which ones apply to both variables
getInts <- function(x, allInts = FALSE, getNames = FALSE, ...){
  if("adjMat" %in% names(x)){x <- x$mods0}
  if("models" %in% names(x)){x <- margCIs(x, ...)}
  ones <- lapply(x, function(z) as.numeric(!(z[, 2] < 0 & z[, 3] > 0)))
  newx <- lapply(1:length(x), function(z) data.frame(x[[z]], ones[[z]]))
  for(i in 1:length(newx)){
    if(nrow(newx[[i]]) >= 1){
      if(any(newx[[i]] == 1)){
        newx[[i]] <- rownames(newx[[i]][newx[[i]][, 4] == 1, ])
      } else {
        newx[[i]] <- NA
      }
    } else {
      newx[[i]] <- NA
    }
  }
  ints <- matrix(0, length(x), length(x))
  for(i in 1:length(x)){ints[i, match(newx[[i]], names(x))] <- 1}
  dimnames(ints) <- rep(list(names(x)), 2)
  zints <- t(ints) * ints
  if(allInts == FALSE & any(zints == 1)){
    f1 <- function(x){
      n <- ncol(x) - 1
      vals <- function(z){
        z1 <- (z * (z - 1))/2 + z + 1
        z1:(z1 + z)
      }
      y <- list()
      for(i in 1:n){y[[i]] <- vals(i - 1)}
      lapply(y, function(z) x[upper.tri(x)][z])
    }
    f2 <- function(x){
      vs <- colnames(x)
      x2 <- f1(x)
      x3 <- list()
      for(i in 1:length(x2)){
        if(!any(x2[[i]] == 1)){
          x3[[i]] <- NA
        } else {
          x3[[i]] <- vs[which(x2[[i]] == 1)]
        }
      }
      x3 <- as.list(unlist(x3[!is.na(x3)]))
      x4 <- rep(2:(length(x2) + 1), sapply(x2, sum))
      lapply(1:length(x3), function(z) c(vs[x4[z]], x3[[z]]))
    }
    if(getNames){return(f2(zints))}
  }
  if(allInts){return(ints)} else {return(zints)}
}

##### condEffects: generate values of beta for X conditioned on certain values of Z
condEffects <- function(mods, x = NULL, xn = NULL, data = NULL, alpha = .05, 
                        adjCI = FALSE, saveMods = TRUE){
  if("SURnet" %in% c(names(mods), names(attributes(mods)))){
    if(!"SURnet" %in% names(mods)){stop("Need 'SURfit' mods")}
    if(!"mnet" %in% names(attributes(mods))){stop("Must have only one exogenous moderator")}
    fitobj <- mods$SURfit$eq
    net <- mods$SURnet
    ynames <- names(net$mods)
    mname <- net$call$moderators
    vars <- net$interactions$coefvars
    vars0 <- as.matrix(vars[, 5:7])
    vars <- as.matrix(vars[, 1:4])
    dfs <- unname(sapply(fitobj, '[[', "df.residual")[match(vars[, 1], ynames)])
    xr <- range(net$data$X[, mname])
    if(length(xn) == 1){x <- seq(xr[1], xr[2], length.out = xn)}
    if(is.null(x)){x <- seq(xr[1], xr[2], length.out = 100)}
    mats <- lapply(seq_along(x), function(i){
      out <- SURnet(fit = mods$SURfit, dat = net$data, m = mname, mval = x[i])
      out <- out$temporal$adjMat
      colnames(out) <- gsub("[.]lag1[.]$", "", colnames(out))
      return(out)
    })
    margSE <- function(x, vars){sqrt(vars[1] + ((x^2) * vars[2]) + (2 * x * vars[3]))}
    dats <- lapply(seq_len(nrow(vars)), function(i){
      D <- data.frame(x, matrix(NA, ncol = 4, nrow = length(x)))
      D[, 2] <- sapply(mats, function(z) z[vars[i, 1], vars[i, 2]])
      D[, 3] <- margSE(x = x, vars = vars0[i, ])
      outlog <- capture.output({fdr <- interactionTest::fdrInteraction(
        D[, 2], D[, 3], df = dfs[i], level = 1 - alpha)})
      D[, 4] <- D[, 2] - (ifelse(adjCI, fdr, qnorm(1 - alpha/2)) * D[, 3])
      D[, 5] <- D[, 2] + (ifelse(adjCI, fdr, qnorm(1 - alpha/2)) * D[, 3])
      colnames(D) <- c("x", "y", "se", "lower", "upper")
      return(D)
    })
    yints <- unique(vars[, 1])
    Y <- lapply(seq_along(yints), function(i){
      whichy <- which(vars[, 1] == yints[i])
      out <- dats[whichy]
      names(out) <- vars[whichy, 2]
      return(out)
    })
    names(Y) <- gsub("[.]y$", "", yints)
    attributes(Y)[c("moderator", "alpha", "SURnet")] <- list(mname, alpha, TRUE)
    if(saveMods){attributes(Y)$mods <- net$data$X[, mname]}
  } else {
    if("adjMat" %in% names(mods)){mods <- mods$mods0}
    if(!any(grepl(":", unlist(sapply(mods$models, function(z) names(coef(z))))))){stop("No interaction terms in the models")}
    if(is.null(data)){data <- mods$dat[, -ncol(mods$dat)]}
    xr <- range(mods$dat[, ncol(mods$dat)])
    if(!is.null(xn)){
      stopifnot(length(xn) == 1)
      x <- seq(xr[1], xr[2], length.out = xn)
    }
    if(is.null(x)){
      x <- ifelse(dim(table(mods$dat[, ncol(mods$dat)])) <= 2, list(c(0, 1)), ifelse(nrow(data) > 100, 100, nrow(data)))[[1]]
      if(length(x) == 1){x <- seq(xr[1], xr[2], length.out = x)}
    }
    mats <- list()
    for(i in seq_along(x)){mats[[i]] <- t(modNet(mods, data, mval = x[i], nsims = 0)$nodewise$adjNW)}
    net0 <- modNet(mods, data)
    p <- ncol(net0$adjMat)
    vars <- net0$interactions$coefvars
    if("varMods" %in% names(attributes(mods$models))){
      inds <- net0$interactions$inds
      inds[, 2] <- match(gsub(":.*", "", inds[, 2]), colnames(data))
      n <- nrow(inds)
      df <- sapply(mods$models, function(z) z$df.residual)[inds[, 1]]
    } else {
      n <- (p * (p - 1))
      inds1 <- net0$interactions$inds
      inds2 <- cbind(inds1[, 2], inds1[, 1])
      inds <- rbind(inds1, inds2)
      inds <- inds[order(inds[, 2]), ]
      inds <- inds[order(inds[, 1]), ]
      inds3 <- cbind(inds[, 1], rep(c(1:(p - 1)), p))
      df <- rep(nrow(data) - (2 * p), n)
    }
    margSE <- function(x, vars){sqrt(vars[1] + ((x^2) * vars[2]) + (2 * x * vars[3]))}
    fdr <- c()
    dats <- list()
    for(i in 1:n){
      dats[[i]] <- data.frame(matrix(NA, ncol = 5, nrow = length(x)))
      dats[[i]][, 1] <- x
      dats[[i]][, 2] <- sapply(mats, function(z) z[inds[i, 1], inds[i, 2]])
      if(n < (p * (p - 1))){
        dats[[i]][, 3] <- margSE(x = x, vars = vars[[i]])
      } else {
        dats[[i]][, 3] <- margSE(x = x, vars = vars[[inds3[i, 1]]][[inds3[i, 2]]])
      }
      log <- capture.output({fdr[i] <- interactionTest::fdrInteraction(dats[[i]][, 2], dats[[i]][, 3], df = df[i], level = 1 - alpha)})
      dats[[i]][, 4] <- dats[[i]][, 2] - (ifelse(adjCI, fdr[i], qnorm(1 - alpha/2)) * dats[[i]][, 3])
      dats[[i]][, 5] <- dats[[i]][, 2] + (ifelse(adjCI, fdr[i], qnorm(1 - alpha/2)) * dats[[i]][, 3])
      colnames(dats[[i]]) <- c("x", "y", "se", "lower", "upper")
    }
    Y <- list()
    for(i in 1:p){
      if(i %in% inds[, 1]){
        Y[[i]] <- dats[which(inds[, 1] == i)]
        if(n < (p * (p - 1))){
          names(Y[[i]]) <- colnames(data)[inds[inds[, 1] == i, 2]]
        } else {
          names(Y[[i]]) <- names(vars[[i]])
        }
      } else {
        Y[[i]] <- NA
      }
    }
    names(Y) <- colnames(data)
    attributes(Y)$moderator <- colnames(mods$dat)[ncol(mods$dat)]
    attributes(Y)$alpha <- alpha
    if(saveMods){attributes(Y)$mods <- mods}
  }
  return(Y)
}

