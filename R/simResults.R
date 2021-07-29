##### compareMods: compare simulated networks with originals
compareMods <- function(fits, Sim, threshold = FALSE, Res = NULL, ind = "correlation",
                        gamma = .5, mats = NULL, rule = "OR"){
  ind <- match.arg(tolower(ind), c(
    "cosine", "mse", "rmse", "ssd", "mae", "msd", "correlation"))
  atts <- names(attributes(fits))
  if("between" %in% names(Sim)){
    if(any(diag(Sim$between) == 1)){diag(Sim$between) <- 0}
    if(any(c("mlGVAR", "lmerVAR") %in% atts)){fits <- list(fits)}
    if(is.null(mats)){mats <- c("kappa", "pcc", "beta", "pdc", "between")}
    out <- sapply(fits, function(z){
      K <- lapply(mats, function(k){
        k0 <- net(z, k, threshold, rule)
        if(k == "pdc"){k0 <- t(k0)}
        return(k0)
      })
      out2 <- sapply(seq_along(K), function(i){
        matrixDist(K[[i]], Sim[[i + 1]], ind = ind,
                   directed = isTRUE(mats[i] %in% c("beta", "pdc")))
      })
      names(out2) <- names(Sim)[2:6]
      return(out2)
    })
    if(isTRUE(Res)){
      Res <- graphicalVAR::mlGraphicalVAR(
        Sim$data, Sim$vars, idvar = Sim$idvar,
        subjectNetworks = FALSE, gamma = gamma)
      assign("Res", Res, envir = .GlobalEnv)
    }
    if(is(Res, "mlGraphicalVAR")){
      r1 <- matrixDist(Res$fixedResults$kappa, Sim$fixedKappa, ind = ind)
      r2 <- matrixDist(Res$fixedPCC, Sim$fixedPCC, ind = ind)
      r3 <- matrixDist(Res$fixedResults$beta[, -1], Sim$fixedBeta, ind = ind)
      r4 <- matrixDist(Res$fixedPDC, Sim$fixedPDC, ind = ind)
      r5 <- matrixDist(Res$betweenNet, Sim$between, ind = ind)
      out <- cbind(out, Res = c(r1, r2, r3, r4, r5))
    }
    out <- t(out)
    if(nrow(out) <= 2){
      if(length(fits) == 2){
        rownames(out) <- if(!is.null(names(fits))){names(fits)} else {c("fit1", "fit2")}
      } else {
        rownames(out) <- switch(nrow(out), "fit", c("fit", "Res"))
      }
    }
  } else if(is(Sim, "mlVARsim")){
    if(any(c("mlGVAR", "lmerVAR") %in% atts)){fits <- list(fits)}
    x <- c("temporal", "contemporaneous", "between")
    Sim <- lapply(x, function(z){
      #z1 <- t(as.matrix(mlVAR::getNet(Sim, z)))
      z1 <- net(Sim, z)
      dimnames(z1) <- NULL
      return(z1)
    })
    out <- sapply(fits, function(z){
      nets <- lapply(x, function(n){
        X <- switch(2 - is(z, "mlVAR"), mlVAR::getNet, net)
        return(X(z, n, threshold = threshold, rule = rule))
      })
      out2 <- setNames(mapply(
        matrixDist, nets, Sim, ind = ind,
        directed = c(TRUE, FALSE, FALSE)
      ), x)
      return(out2)
    })
    out <- t(out)
  } else {
    mats <- c("kappa", "beta", "PCC", "PDC")
    nets <- lapply(mats, function(z){
      z0 <- net(fits, z, threshold, rule)
      if(z == "PDC"){z0 <- t(z0)}
      return(z0)
    })
    names(nets) <- mats
    if(length(Sim) > length(mats)){
      Sim <- list(kappa = Sim$fixedKappa, beta = Sim$fixedBeta,
                  PCC = Sim$fixedPCC, PDC = Sim$fixedPDC)
    }
    out <- mapply(matrixDist, nets, Sim, ind = ind,
                  directed = rep(c(FALSE, TRUE), 2))
    if(is(Res, "graphicalVAR")){
      Res1 <- Res[mats]
      Res1$beta <- Res1$beta[, -1]
      out1 <- mapply(matrixDist, Res1, Sim, ind = ind,
                     directed = rep(c(FALSE, TRUE), 2))
      out <- rbind(out, out1)
      rownames(out) <- c("fit", "Res")
    }
  }
  output <- list(out = out)
  if(ifelse(is(out, "matrix"), ifelse(nrow(out) > 1, TRUE, FALSE), FALSE)){
    rnks <- apply(out, 2, function(z){
      order(order(z, decreasing = ifelse(startsWith(ind, "c"), TRUE, FALSE)))
    })
    rownames(rnks) <- rownames(out)
    output$ranks <- rnks
  }
  names(output)[1] <- ifelse(ind == "correlation", "cor", ind)
  return(output)
}

##### performance: sensitivity, specificity, precision, accuracy
performance <- function(est, trueMod, threshold = FALSE, combine = FALSE,
                        inds = "all", rule = "OR", getVals = FALSE,
                        mcc = TRUE, rmNAs = TRUE){
  if(combine){
    nn <- switch(2 - is.null(names(est)), paste0("fit", 1:length(est)), names(est))
    xx <- lapply(lapply(est, performance, trueMod, threshold,
                        inds = inds, rule = rule), function(z) data.frame(t(z)))
    nn2 <- colnames(xx[[1]])
    xx <- lapply(seq_along(nn2), function(z){
      data.frame(t(do.call(cbind, lapply(xx, '[', z))), row.names = nn)})
    names(xx) <- nn2
    return(xx)
  }
  if(all(inds == "all")){
    inds <- c("kappa", "beta", "PCC", "PDC")
    atts <- names(attributes(est))
    if(any(c("mlGVAR", "lmerVAR") %in% atts) | is(est, "mlGraphicalVAR")){
      inds <- c(inds, "between")
    }
    if(is(est, "matrix") & is(trueMod, "matrix")){inds <- "beta"}
  }
  inds <- match.arg(tolower(inds), c('kappa', 'beta', 'pcc', 'pdc', 'between'), several.ok = TRUE)
  if(any(grepl('^p', inds))){inds[grepl('^p', inds)] <- toupper(inds[grepl('^p', inds)])}
  tt <- FALSE
  if(length(inds) > 1){
    if(isTRUE(attr(trueMod, "simMLgvar")) | is(trueMod, "simMLgvar")){
      trueMod <- setNames(lapply(inds, function(z){
        n2 <- as.matrix(net(trueMod, z))
        if(z == "between"){diag(n2) <- 0}
        dimnames(n2) <- NULL
        return(n2)
      }), inds)
      #tt <- !is(est, "mlGraphicalVAR")
    } else if(is(trueMod, "mlVARsim")){
      inds <- c("temporal", "contemporaneous", "between")
      trueMod <- setNames(lapply(inds, function(z){
        #n2 <- t(as.matrix(mlVAR::getNet(trueMod, z)))
        n2 <- net(trueMod, z)
        dimnames(n2) <- NULL
        return(n2)
      }), inds)
    }
    nets1 <- setNames(lapply(inds, function(z){
      n1 <- as.matrix(net(est, z, threshold, rule))
      if(tt & z == "PDC"){n1 <- t(n1)}
      dimnames(n1) <- NULL
      return(n1)
    }), inds)
  } else {
    if(!is(est, "matrix")){est <- net(est, inds)}
    nets1 <- setNames(list(est), inds)
    trueMod <- setNames(list(trueMod), inds)
  }
  p <- unique(sapply(trueMod, nrow))
  trueVals <- lapply(inds, function(z){
    if(z %in% c("PCC", "between", "contemporaneous", "kappa")){
      z1 <- trueMod[[z]][lower.tri(trueMod[[z]])]
    } else {
      z1 <- as.vector(trueMod[[z]])
    }
    return(list(truePos = which(z1 != 0), trueNeg = which(z1 == 0)))
  })
  estVals <- lapply(inds, function(z){
    if(z %in% c("PCC", "between", "contemporaneous", "kappa")){
      z1 <- nets1[[z]][lower.tri(nets1[[z]])]
    } else {
      z1 <- as.vector(nets1[[z]])
    }
    return(list(estPos = which(z1 != 0), estNeg = which(z1 == 0)))
  })
  tp <- mapply(function(x1, x2){
    length(intersect(x1[[1]], x2[[1]]))}, estVals, trueVals)
  tn <- mapply(function(x1, x2){
    length(intersect(x1[[2]], x2[[2]]))}, estVals, trueVals)
  fp <- sapply(lapply(estVals, '[[', 1), length) - tp
  fn <- sapply(lapply(estVals, '[[', 2), length) - tn
  sensitivity <- tp/(tp + fn)
  specificity <- tn/(tn + fp)
  precision <- tp/(tp + fp)
  accuracy <- (tp + tn)/(tp + fp + tn + fn)
  out <- cbind.data.frame(sensitivity, specificity, precision, accuracy)
  if(mcc){
    numerator <- ((tp * tn) - (fp * fn))
    denom <- (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    mcc <- numerator/ifelse(denom == 0, 1, sqrt(denom))
    out <- cbind.data.frame(out, mcc = mcc)
  }
  if(any(is.na(out)) & rmNAs){out[is.na(out)] <- 0}
  rownames(out) <- inds
  if(getVals){
    out2 <- cbind.data.frame(tp, tn, fp, fn)
    rownames(out2) <- inds
    out <- list(performance = out, indices = out2,
                estVals = estVals, trueVals = trueVals)
  }
  if(tt){out <- t(out)}
  return(out)
}

##### getNets
getNets <- function(obj, interactions = TRUE, dat = NULL){
  if(!is.null(dat)){obj <- list(call = obj$call, dat = dat, fit = obj)}
  x <- switch(2 - all(c("dat", "fit") %in% names(obj)), 1, length(obj))
  if(x == 1){obj <- list(obj)}
  allNets <- c("temporal", "contemporaneous", "between")
  calls <- simnets <- fitnets <- list()
  for(i in seq_len(x)){
    out <- obj[[i]]
    calls[[i]] <- switch(2 - ("call" %in% names(out)), out$call, list())
    simnets[[i]] <- setNames(vector("list", 1), "nets0")
    simnets[[i]][[1]] <- setNames(lapply(allNets, function(z) net(out$dat, z)), allNets)
    if(is.null(out$call$m)){interactions <- FALSE}
    if(interactions){simnets[[i]][[1]]$interactions <- netInts(out$dat)}
    fitnets[[i]] <- setNames(vector("list", 3), paste0("nets", 0:2))
    for(j in 1:3){
      thresh <- c(FALSE, TRUE, TRUE); rule <- c("or", "or", "and")
      fitnets[[i]][[j]] <- setNames(lapply(allNets, function(z) net(out$fit, z, thresh[j], rule[j])), allNets)
      if(interactions){fitnets[[i]][[j]]$interactions <- netInts(out$fit, threshold = thresh[j])}
    }
  }
  if(!is.null(names(obj))){names(calls) <- names(simnets) <- names(fitnets) <- names(obj)}
  #if(x == 1){calls <- calls[[1]]; simnets <- simnets[[1]]; fitnets <- fitnets[[1]]}
  output <- list(calls = calls, simnets = simnets, fitnets = fitnets)
  return(output)
}

##### results
results <- function(obj, nets = 1, ind = "perf", X = NULL, abbv = 10){
  nn <- c("temporal", "contemporaneous", "between", "interactions")
  if(any(sapply(lapply(obj$calls, '[[', "m"), is.null))){nn <- nn[1:3]}
  if(isTRUE(attr(obj, "SURnet"))){nn <- setdiff(nn, "between")}
  if(ind %in% c("cosine", "cor", "mse", "mae")){
    Y1 <- lapply(1:3, function(z){
      yy <- do.call(rbind, lapply(seq_along(obj$calls), function(i){
        sapply(seq_along(nn), function(j){matrixDist(
          obj$fitnets[[i]][[z]][[j]], obj$simnets[[i]]$nets0[[j]],
          ind, isTRUE(nn[j] %in% c("temporal", "interactions")))})
      }))
      colnames(yy) <- nn
      if(!is.null(names(obj$calls))){
        rownames(yy) <- abbreviate(names(obj$calls), minlength = abbv)}
      return(yy)
    })
  } else {
    Y1 <- lapply(1:3, function(z){
      setNames(lapply(seq_along(nn), function(j){
        yy <- do.call(rbind, lapply(seq_along(obj$calls), function(i){
          performance(obj$fitnets[[i]][[z]][[j]], obj$simnets[[i]]$nets0[[j]])}))
        if(!is.null(names(obj$calls))){
          rownames(yy) <- abbreviate(names(obj$calls), minlength = abbv)}
        return(yy)
      }), nn)
    })
  }
  if(is.null(X)){
    return(Y1[[nets]])
  } else {
    k <- gsub("_.*", "", names(obj$calls)[1])
    X <- subset(X, FUN == k)
    if(ind %in% c("cosine", "cor", "mse", "mae")){
      Y1 <- cbind(X, Y1[[nets]])
    } else {
      Y1 <- setNames(lapply(seq_along(nn), function(z) cbind(X, Y1[[nets]][[z]])), nn)
    }
    return(Y1)
  }
}
