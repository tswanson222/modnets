################################################################################
##### nodewise: nodewise regression, with option to simulate/add moderators
nodewise <- function(data, mods = NULL, varMods = NULL, lambda = "min", 
                     center = TRUE, scale = FALSE){
  data <- data.frame(data)
  intMatrix <- function(data, mods = NULL){
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
    return(eqs)
  }
  if(!is.null(mods)){
    if(class(mods) == "numeric"){
      mod <- list(data[, mods])
      names(mod) <- colnames(data)[mods]
      data <- data[, -mods]
      dat <- data.frame(data, mod)
      mods <- mod
    }
    dat <- data.frame(data, mods)
    if(!is.null(varMods)){
      lambda <- ifelse(match.arg(tolower(lambda), c("min", "1se")) == "min", 1, 2)
      ints <- as.list(paste(names(varMods), "~", lapply(varMods, function(z) paste0(z[[lambda]], collapse = " + "))))
    } else {
      ints <- intMatrix(data, mods)
    }
  } else {
    dat <- data.frame(data)
    if(!is.null(varMods)){
      ints <- as.list(paste(names(varMods), "~", lapply(varMods, function(z) paste0(z[[1]], collapse = " + "))))
    } else {
      ints <- intMatrix(data)
    }
  }
  if(center != FALSE){
    if(!is.null(mods) & dim(table(mods[[1]])) <= 2 | center != TRUE){
      dat[, -ncol(dat)] <- apply(dat[, -ncol(dat)], 2, scale, TRUE, scale)
    } else {
      dat <- apply(dat, 2, scale, TRUE, scale)
    }
  }
  dat <- data.frame(dat)
  m <- lapply(ints, function(z) lm(z, dat))
  mm <- lapply(lapply(m, coef), function(z) z[which(names(z) %in% colnames(data))])
  ps1 <- lapply(m, function(z){
    z2 <- t(data.frame(summary(z)$coefficients))
    z3 <- z2[4, which(names(coef(z)) %in% colnames(data)), drop = FALSE]
    rownames(z3) <- NULL
    return(z3[1, ])
  })
  if(!is.null(mods)){
    m2 <- lapply(lapply(m, coef), function(z) z[grep(":", names(z))])
    ps2 <- lapply(m, function(z){
      z2 <- t(data.frame(summary(z)$coefficients))
      z3 <- z2[4, grep(":", colnames(z2)), drop = FALSE]
      rownames(z3) <- NULL
      return(z3[1, ])
    })
    mname <- colnames(dat)[-c(1:ncol(data))]
    mx <- paste0(colnames(data), ":", mname)
    psx <- bx <- suppressWarnings(diag(mx))
    rownames(psx) <- rownames(bx) <- colnames(data)
    colnames(psx) <- colnames(bx) <- mx
    diag(psx) <- diag(bx) <- 1
    m3 <- lapply(m, function(z) summary(z)$coefficients)
    names(m3) <- colnames(data)
    bm <- do.call(rbind, lapply(m3, function(z) z[rownames(z) == mname, ]))
  }
  b <- ps <- matrix(NA, nrow = ncol(data), ncol = ncol(data))
  for(i in 1:ncol(data)){
    b[i, match(names(mm[[i]]), colnames(data))] <- mm[[i]]
    ps[i, match(names(ps1[[i]]), colnames(data))] <- ps1[[i]]
    if(!is.null(mods)){
      bx[i, match(names(m2[[i]]), mx)] <- m2[[i]]
      psx[i, match(names(ps2[[i]]), mx)] <- ps2[[i]]
    }
  }
  diag(b) <- diag(ps) <- 1
  if(any(is.na(b))){b[is.na(b)] <- 0}
  if(any(is.na(ps))){ps[is.na(ps)] <- 0}
  out <- list(models = m, B = list(b = b, ps = ps))
  if(is.null(mods)){
    attributes(out$models)$noMods <- TRUE
  } else {
    if(nrow(bm) >= 1){out$Bm <- bm}
    out$Bx <- list(bx = bx, px = psx)
  }
  if(!is.null(varMods)){
    if(!any(grepl(":", unlist(sapply(out$models, function(z) names(coef(z))))))){
      attributes(out$models)$noMods <- TRUE
      out$Bx <- NULL
    }
    attributes(out$models)$varMods <- c("min", "1se")[lambda]
  }
  out$dat <- dat
  out
}

################################################################################
##### modNet: create moderated network from nodewise regression models
modNet <- function(models, data = NULL, threshold = FALSE, rule = "AND", 
                   mval = NULL, useCIs = FALSE, nsims = 5000, mlty = 2){
  if("models" %in% names(models)){
    mods0 <- models
    models <- models$models
    if(is.null(data)){
      if("noMods" %in% names(attributes(models)) & length(models) == ncol(mods0$dat)){
        data <- mods0$dat
      } else {
        data <- mods0$dat[, -ncol(mods0$dat)]
      }
    }
  }
  p <- ncol(data)
  vs <- colnames(data)
  mods <- lapply(models, function(z){
    z2 <- matrix(coef(z), ncol = 1)
    rownames(z2) <- names(coef(z))
    return(z2)
  })
  mods2 <- lapply(mods, function(z) z[which(rownames(z) %in% vs), ])
  pvals <- lapply(models, function(z){
    z2 <- t(data.frame(summary(z)$coefficients))
    z3 <- z2[4, which(names(coef(z)) %in% vs), drop = FALSE]
    rownames(z3) <- NULL
    return(z3[1, ])
  })
  ses <- lapply(models, function(z){
    z2 <- t(data.frame(summary(z)$coefficients))
    z3 <- z2[2, which(names(coef(z)) %in% vs), drop = FALSE]
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
  if("noMods" %in% names(attributes(models))){
    inds <- ints <- mval <- vars1 <- NULL
  } else if("varMods" %in% names(attributes(models))){
    inds0 <- lapply(mods, function(z) rownames(z)[grepl(":", rownames(z))])
    inds1 <- unlist(lapply(inds0, function(z) gsub(":.*", "", z)))
    inds <- data.frame(outcome = rep(1:p, sapply(inds0, length)), interaction = unlist(inds0))
    vars <- lapply(models, vcov)
    vars1 <- ints <- list()
    for(i in 1:nrow(inds)){
      ints[[i]] <- mods[[inds[i, 1]]][inds[i, 2], 1]
      vars1[[i]] <- c(vars[[inds[i, 1]]][inds1[i], inds1[i]], vars[[inds[i, 1]]][inds[i, 2], inds[i, 2]], vars[[inds[i, 1]]][inds1[i], inds[i, 2]])
      names(vars1[[i]]) <- c(inds1[i], inds[i, 2], "cov")
    }
  } else {
    inds <- t(combn(1:length(mods), 2))
    ints <- vector("list", nrow(inds))
    for(i in 1:nrow(inds)){
      ints[[i]][1] <- mods[[inds[i, 1]]][ncol(data) + inds[i, 2], ]
      ints[[i]][2] <- mods[[inds[i, 2]]][ncol(data) + inds[i, 1] + 1, ]
    }
    vars <- lapply(models, vcov)
    vars1 <- list()
    for(i in 1:length(models)){
      vars1[[i]] <- vector("list", p - 1)
      for(j in 1:(p - 1)){
        vars1[[i]][[j]] <- c(vars[[i]][j + 1, j + 1], vars[[i]][j + p + 1, j + p + 1], vars[[i]][j + 1, j + p + 1])
        names(vars1[[i]][[j]]) <- c(colnames(vars[[i]])[j + 1], colnames(vars[[i]])[j + p + 1], "cov")
      }
      names(vars1[[i]]) <- colnames(vars[[i]])[2:p]
    }
    names(vars1) <- colnames(data)
  }
  b2ggm <- function(b, rule = "AND"){
    rule <- match.arg(tolower(rule), c("and", "or"))
    if(rule == "and"){
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
    return((b + t(b))/2)
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
    dfs <- matrix(sapply(models, function(z) z$df), ncol = p, nrow = p)
    pvals2 <- (2 * pt(abs(b/ses2), df = dfs, lower.tail = FALSE))
    if(any(is.na(pvals2))){pvals2[is.na(pvals2)] <- 1}
  }
  if(threshold != FALSE){
    if(threshold == TRUE){threshold <- .05}
    b <- b * ifelse(pvals2 <= threshold, 1, 0)
  }
  bb <- b2ggm(b, rule = rule)
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
    threshold <- ifelse(threshold == FALSE, .05, threshold)
    if(useCIs & nsims > 0){
      cis <- margCIs(mods = mods0, alpha = threshold, nsims = nsims)
      modEdgesNW <- getInts(x = cis, allInts = TRUE)
    } else {
      modEdgesNW <- ifelse(mods0$Bx$px <= threshold, 1, 0)
      if(any(mods0$Bx$px == 0)){modEdgesNW <- modEdgesNW * ifelse(mods0$Bx$px == 0, 0, 1)}
      colnames(modEdgesNW) <- rownames(modEdgesNW)
    }
    modEdges <- t(modEdgesNW) * modEdgesNW
    modEdgesNW <- (modEdgesNW * (mlty - 1)) + 1
    modEdges <- (modEdges * (mlty - 1)) + 1
    rule <- match.arg(tolower(rule), c("and", "or"))
    if(rule == "or"){modEdges <- modEdgesNW}
  } else {
    modEdges <- modEdgesNW <- matrix(1, p, p)
  }
  diag(modEdges) <- diag(modEdgesNW) <- 0
  dimnames(b) <- dimnames(bb) <- dimnames(pvals2) <- rep(list(colnames(data)), 2)
  names(mods) <- names(models) <- colnames(data)
  out <- list(adjMat = bb, edgeColors = getEdgeColors(bb), modEdges = t(modEdges),
              nodewise = list(adjNW = t(b), edgeColsNW = getEdgeColors(t(b)), 
                              pvalsNW = t(pvals2), modEdgesNW = t(modEdgesNW)),
              interactions = list(inds = inds, ints = ints, coefvars = vars1), 
              mods = mods, fitobj = models, data = data)
  if("noMods" %in% names(attributes(models))){
    out[c("interactions", "modEdges")] <- NULL
    out[["nodewise"]]["modEdgesNW"] <- NULL
  } else {
    if(useCIs){out$interactions$cis <- cis}
    attributes(out)$moderator <- colnames(mods0$dat)[ncol(mods0$dat)]
  }
  if(!is.null(mval)){attributes(out)$mval <- mval}
  attributes(out)$ggm <- TRUE
  out
}

################################################################################
##### margCIs: retrieve CIs for the effect of Z on the coefficient relating X to Y
margCIs <- function(mods, data = NULL, modname = NULL, alpha = .05, nsims = 5000, 
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
                                 t = tval, p = 2 * pt(abs(tval), df = mods[[i]]$df, lower.tail = F))
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

################################################################################
##### getInts: retrieve interactions effects; particularly those that apply to both variables
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

################################################################################
################################################################################
##### condEffects: generate values of beta for X conditioned on certain values of Z
condEffects <- function(mods, x = NULL, xn = NULL, data = NULL, alpha = .05, 
                        adjCI = FALSE, saveMods = TRUE){
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
    df <- sapply(mods$models, function(z) z$df)[inds[, 1]]
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
  Y
}

################################################################################
##### condPlot: plot conditional effects
condPlot <- function(out, to, from, swap = FALSE, avg = FALSE, compare = NULL, 
                     hist = FALSE, xlab = NULL, mods = NULL, nsims = 5000, 
                     xn = NULL, getCIs = FALSE, discrete = FALSE, midline = TRUE){
  require(ggplot2)
  if("adjMat" %in% names(out)){out <- out$mods0}
  if("models" %in% names(out)){out <- condEffects(out, xn = xn)}
  if(is.null(xlab)){xlab <- attributes(out)$moderator}
  if(length(to) == 2){from <- to[2]; to <- to[1]}
  if(class(to) == "numeric"){to <- names(out)[to]}
  if(class(from) == "numeric"){from <- names(out)[from]}
  if(swap == TRUE){
    toX <- to
    to <- from
    from <- toX
  }
  data <- out[[to]][[from]]
  stopifnot(!is.null(data))
  stopifnot(nrow(data) > 1)
  if(avg == TRUE){
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
  if(avg == FALSE){
    ylab <- paste0(fr2, " ---> ", to2)
    mainlab <- paste0("Estimated Effect of ", fr2, " on ", to2, " across levels of ", xlab)
  } else {
    ylab <- paste0(fr2, " ~ ", to2)
    mainlab <- paste0("Mean Relation between ", fr2, " and ", to2, " across levels of ", xlab)
  }
  alpha <- attributes(out)$alpha
  if("mods" %in% names(attributes(out))){mods <- attributes(out)$mods}
  if(!is.null(mods)){
    ci_diff <- margCIs(mods = mods, alpha = alpha, nsims = nsims, compare = compare)
    if(avg){ci_diff2 <- ci_diff[[from]][to, ]}
    ci_diff <- ci_diff[[to]][from, ]
  } else {
    ci_diff <- NULL
  }
  if(nrow(data) != 2 & discrete == FALSE){
    if(hist == FALSE){
      pp <- ggplot(data = data, aes(x = x, y = y)) + 
        geom_line(color = "red") +
        geom_ribbon(alpha = .2, aes(ymin = lower, ymax = upper))
      if(min(data$lower) < 0 & max(data$upper) > 0){pp <- pp + geom_hline(yintercept = 0, linetype = 2)}
      pp <- pp + xlab(xlab) + ylab(ylab) + ggtitle(mainlab) + theme_bw()
    } else {
      var2_dt <- mods$dat[, ncol(mods$dat)]
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

################################################################################
##### intsPlot: plot CI(max - min) for each predictor with interactions against Y
intsPlot <- function(out, y, moderator = NULL, nsims = 5000, alpha = .05){
  require(ggplot2)
  if("adjMat" %in% names(out)){out <- out$mods0}
  if("models" %in% names(out)){out <- margCIs(out, modname = moderator, nsims = nsims, alpha = alpha)}
  if(is.null(moderator)){moderator <- attributes(out)$moderator}
  if(class(y) == "character"){
    if(y == "all"){
      y0 <- as.vector(sapply(seq_along(out), function(z) 
        paste0(names(out)[z], ":", rownames(out[[z]]))))
      dd <- suppressWarnings(data.frame(do.call(rbind, out)))
      rownames(dd) <- y0
    } else if(y == "sig"){
      out2 <- lapply(seq_along(out), function(z){
        z2 <- data.frame(out[[z]], sig = as.numeric(!(out[[z]][, 2] < 0 & out[[z]][, 3] > 0)), 
                         check.names = FALSE)
        return(z2[z2$sig == 1, -4])
      })
      names(out2) <- names(out)
      if(any(sapply(out2, nrow) == 0)){
        out2 <- out2[-which(sapply(out2, nrow) == 0)]
      }
      y0 <- unlist(as.vector(sapply(seq_along(out2), function(z)
        paste0(names(out2)[z], ":", rownames(out2[[z]])))))
      dd <- suppressWarnings(data.frame(do.call(rbind, out2)))
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

################################################################################
##### plotMods: plot networks at different levels of moderator
plotMods <- function(nets, nodewise = FALSE, elsize = 2, vsize = NULL, 
                     elabs = TRUE, predict = TRUE, layout = NULL, ...){
  if(is.null(vsize)){
    vsize <- (exp(-ncol(nets[[1]]$data)/80) * 8) + 1
    vsize <- vsize + (exp(-length(nets)/80) * 8)/length(nets)
  }
  getLayout <- function(x){
    stopifnot(class(x) == "list")
    averageLayout(lapply(x, function(z) plotNet(z, plot = F)))
  }
  getMax <- function(x, n = FALSE){
    if(n == FALSE){
      max(unlist(lapply(x, function(z) max(abs(z$adjMat)))))
    } else {
      max(unlist(lapply(x, function(z) max(abs(z$nodewise$adjNW)))))
    }
  }
  if(is.null(layout)){layout <- getLayout(nets)}
  mx <- getMax(nets, nodewise)
  moderator <- Hmisc::capitalize(attr(nets[[1]], "moderator"))
  vals <- unname(sapply(nets, attr, "mval"))
  layout(t(1:length(vals)))
  for(i in 1:length(vals)){
    plotNet(nets[[i]], predict = predict, layout = layout, elabs = elabs,
            elsize = elsize, nodewise = nodewise, maximum = mx, 
            title = paste0(moderator, " = ", vals[i]), vsize = vsize, ...)
  }
}
