resample <- function(data, niter = 10, m = NULL, criterion = "AIC", sampMethod = "split", 
                     varMethod = "glinternet", rule = "OR", gamma = .25, nfolds = 10, 
                     nlam = 50, which.lam = "min", threshold = FALSE, bonf = TRUE, 
                     alpha = .05, exogenous = TRUE, split = .5, center = TRUE, 
                     varSeed = NULL, seeds = NULL, verbose = TRUE, binary = NULL, 
                     saveMods = FALSE, type = "gaussian", ...){
  t1 <- Sys.time()
  N <- nrow(data)
  data <- data.frame(data)
  criterion <- toupper(match.arg(
    tolower(criterion), 
    choices = c("cv", "aic", "bic", "ebic", "cp", "rss", "adjr2", "rsq", "r2")))
  if(!criterion %in% c("AIC", "BIC", "EBIC", "CV")){varMethod <- "subset"}
  varMethod <- match.arg(varMethod, c("glinternet", "hiernet", "subset"))
  if(!is.null(m)){if(all(is.na(m)) | all(m == 0)){m <- NULL}}
  if(!is.null(binary)){
    type <- rep("g", ifelse(is.null(m) | !exogenous, ncol(data), ncol(data) - 1))
    type[intersect(seq_along(type), binary)] <- "c"
  }
  if(varMethod == "glinternet" & is.null(m)){m <- 1:ncol(data)}
  if(!is.null(m)){if(length(m) > 1){exogenous <- FALSE}}
  sampMethod <- match.arg(tolower(sampMethod), c("split", "bootstrap", "stability"))
  preout <- list(niter = niter, criterion = criterion, sampMethod = sampMethod,
                 varMethod = varMethod, moderators = m, rule = rule, alpha = alpha, 
                 bonf = bonf, gamma = gamma, nfolds = nfolds, which.lam = which.lam,
                 split = split, center = center, exogenous = exogenous, 
                 varSeed = varSeed, type = type)
  #args <- list(...)
  #if(length(args) != 0){preout <- append(preout, args)}
  if(criterion != "CV"){preout[c("nfolds", "which.lam")] <- NULL; which.lam <- "min"}
  if(!is.null(seeds)){preout$seeds <- seeds} else {seeds <- 1:niter}
  if(nlam != 50 & varMethod != "subset"){preout$nlam <- nlam}
  if(sampMethod == "bootstrap"){preout$split <- NULL}
  if(criterion != "EBIC"){preout$gamma <- NULL}
  if(is.null(varSeed)){preout$varSeed <- NULL}
  lam <- ifelse(grepl("min", which.lam), 1, 2)
  sampInd <- samps <- vars <- vars1 <- fits <- list()
  if(sampMethod != "bootstrap"){
    if(split <= 0 | split >= 1){stop("split size must be between 0 and 1")}
    train <- list(); test <- list()
    n <- floor(N * split)
    for(i in 1:niter){
      set.seed(seeds[i])
      sampInd[[i]] <- sample(1:N, n, replace = FALSE)
      train[[i]] <- data[sampInd[[i]], ]
      test[[i]] <- data[-sampInd[[i]], ]
      samps[[i]] <- vector("list", length = 2)
      samps[[i]][[1]] <- train[[i]]
      samps[[i]][[2]] <- test[[i]]
      names(samps[[i]]) <- c("train", "test")
      if(verbose){
        if(sampMethod == "stability"){
          if(i == 1){pb <- txtProgressBar(min = 0, max = niter, char = "-", style = 3)}
        } else {
          if(i == 1){cat("\n")}
          cat("************ Sample Split: ", i, "/", niter, " ************\n", sep = "")
        }
      }
      set.seed(seeds[i])
      vars[[i]] <- varSelect(dat = train[[i]], m = m, criterion = criterion, 
                             method = varMethod, nfolds = nfolds, gamma = gamma, 
                             nlam = nlam, type = type, seed = varSeed, all = !exogenous,
                             verbose = ifelse(sampMethod == "stability", FALSE, verbose)) 
      if(sampMethod == "stability"){
        vars1[[i]] <- varSelect(dat = test[[i]], m = m, criterion = criterion,
                                method = varMethod, nfolds = nfolds, gamma = gamma,
                                nlam = nlam, type = type, seed = varSeed,
                                all = !exogenous, verbose = FALSE)
        if(verbose){setTxtProgressBar(pb, i); if(i == niter){close(pb)}}
      } else {
        fits[[i]] <- fitNetwork(data = test[[i]], moderators = m, type = vars[[i]], 
                                threshold = FALSE, which.lam = which.lam, 
                                gamma = gamma, center = center, rule = rule, 
                                exogenous = exogenous, verbose = FALSE) #, ...)
      }
    }
  }
  if(sampMethod == "bootstrap"){
    for(i in 1:niter){
      set.seed(seeds[i])
      sampInd[[i]] <- sample(1:N, N, replace = TRUE)
      samps[[i]] <- data[sampInd[[i]], ]
      if(verbose == TRUE){
        if(i == 1){cat("\n")}
        cat("************* Bootstrap: ", i, "/", niter, " **************\n", sep = "")
      }
      set.seed(seeds[i])
      vars[[i]] <- varSelect(dat = samps[[i]], m = m, criterion = criterion,
                             method = varMethod, verbose = verbose, gamma = gamma,
                             type = type, seed = varSeed, all = !exogenous)
      fits[[i]] <- fitNetwork(data = samps[[i]], moderators = m, type = vars[[i]],
                              threshold = FALSE, which.lam = which.lam,
                              gamma = gamma, center = center, rule = rule, 
                              exogenous = exogenous, verbose = FALSE) #, ...)
    }
  }
  fit0 <- fitNetwork(data = data, moderators = m, type = type, rule = rule,
                     threshold = FALSE, gamma = gamma, center = center, 
                     exogenous = exogenous, verbose = FALSE) #, ...)
  p <- length(fit0$mods)
  vs <- names(fit0$mods)
  allNames <- lapply(fit0$mods, function(z) rownames(z$model)[-1])
  varMods0 <- lapply(vars, function(z) lapply(z, function(zz) zz$mod0))
  if(sampMethod == "stability"){
    if(is.logical(threshold)){threshold <- .75}
    #dfmax in glmnet: q <- ceiling(sqrt(EV * p * (2 * threshold - 1)))
    varMods1 <- lapply(vars1, function(z) lapply(z, '[[', 1))
    simultvars <- lapply(seq_len(niter), function(z){
      lapply(seq_len(p), function(zz){
        varSamp1 <- !is.na(match(allNames[[zz]], varMods0[[z]][[zz]]))
        varSamp2 <- !is.na(match(allNames[[zz]], varMods1[[z]][[zz]]))
        return(varSamp1 & varSamp2)
      })
    })
    varFreqs <- suppressMessages(suppressWarnings(lapply(seq_len(p), function(z){
      modTerms <- t(sapply(seq_len(niter), function(zz) simultvars[[zz]][[z]]))
      colnames(modTerms) <- allNames[[z]]
      freqs <- colMeans(modTerms)
      coefs1 <- summary(fit0$fitobj[[z]])$coef[-1, , drop = FALSE]
      coefs2 <- confint(fit0$fitobj[[z]], level = 1 - alpha)[-1, , drop = FALSE]
      outFreqs <- data.frame(predictor = factor(allNames[[z]]), select = freqs >= threshold,
                             lower = coefs2[, 1], b = coefs1[, 1], upper = coefs2[, 2], 
                             Pvalue = coefs1[, 4], freq = freqs)
      rownames(outFreqs) <- 1:nrow(outFreqs)
      return(outFreqs)
    })))
    attributes(varFreqs) <- list(
      names = vs, type = type, niter = niter, criterion = criterion, 
      method = varMethod, threshold = threshold, select = "select")
    return(varFreqs)
  } else {
    mod0Freqs <- lapply(1:p, function(z){
      data.frame(table(unlist(lapply(varMods0, function(zz) zz[[z]]))))})
    mod0Freqs2 <- lapply(1:p, function(z) mod0Freqs[[z]][, 2])
    for(i in 1:p){
      names(mod0Freqs2[[i]]) <- as.character(mod0Freqs[[i]][, 1])
      mod0Freqs2[[i]] <- mod0Freqs2[[i]][match(allNames[[i]], names(mod0Freqs2[[i]]))]
      names(mod0Freqs2[[i]]) <- allNames[[i]]
      if(any(is.na(mod0Freqs2[[i]]))){mod0Freqs2[[i]][is.na(mod0Freqs2[[i]])] <- 0}
      mod0Freqs2[[i]] <- mod0Freqs2[[i]]/niter
    }
    mod0Freqs <- lapply(mod0Freqs2, as.matrix, ncol = 1)
    names(mod0Freqs) <- names(fit0$mods)
    mod0sizes <- sapply(varMods0, function(z) sapply(z, length))
    for(i in 1:niter){
      if(any(mod0sizes[, i] != max(mod0sizes[, i]))){
        ms0 <- which(mod0sizes[, i] != max(mod0sizes[, i]))
        for(j in seq_along(ms0)){
          varMods0[[i]][[ms0[j]]] <- c(
            varMods0[[i]][[ms0[j]]], 
            rep("", max(mod0sizes[, i]) - length(varMods0[[i]][[ms0[j]]])))
        }
      }
    }
    mod0coefs <- finalCoefs0 <- list()
    if(verbose){
      npb <- seq(43, 53, by = 2)[sum(niter >= c(10, 100, 1000, 10000, 100000)) + 1]
      cat("\n"); cat(paste0(rep("#", npb), collapse = ""), "\n")
      cat(paste0("Estimating ", (1 - alpha) * 100, "% CIs\n"))
      pb <- txtProgressBar(min = 0, max = niter, style = 1, char = "-", width = npb)
    }
    for(i in 1:niter){
      mod0coefs[[i]] <- suppressWarnings(suppressMessages(
        lapply(fits[[i]]$fitobj, function(z){
          coefs1 <- t(summary(z)$coef[-1, , drop = FALSE])
          coefs2 <- t(confint(z, level = 1 - alpha)[-1, , drop = FALSE])
          return(rbind(coefs1, coefs2))
        })))
      for(j in 1:p){
        rownames(mod0coefs[[i]][[j]]) <- c("b", "se", "t", "P", "lower", "upper")
        mod0coefs[[i]][[j]] <- mod0coefs[[i]][[j]][, match(allNames[[j]], colnames(mod0coefs[[i]][[j]]))]
        colnames(mod0coefs[[i]][[j]]) <- allNames[[j]]
        if(bonf){mod0coefs[[i]][[j]]["P", ] <- pmin(mod0coefs[[i]][[j]]["P", ] * mod0sizes[j, i], 1)}
      }
      if(verbose){setTxtProgressBar(pb, i); if(i == niter){close(pb)}}
    }
    n <- ifelse(sampMethod == "bootstrap", 0, n)
    for(j in 1:p){
      finalCoefs0[[j]] <- vector("list", length = 6)
      finalCoefs0[[j]][[1]] <- t(sapply(mod0coefs, function(z) z[[j]]["b", ]))
      finalCoefs0[[j]][[2]] <- t(sapply(mod0coefs, function(z) z[[j]]["se", ]))
      if(any(is.na(finalCoefs0[[j]][[2]]))){finalCoefs0[[j]][[2]][is.na(finalCoefs0[[j]][[2]])] <- Inf}
      finalCoefs0[[j]][[3]] <- t(sapply(mod0coefs, function(z) z[[j]]["P", ]))
      if(any(is.na(finalCoefs0[[j]][[3]]))){finalCoefs0[[j]][[3]][is.na(finalCoefs0[[j]][[3]])] <- 1}
      finalCoefs0[[j]][[4]] <- t(sapply(mod0coefs, function(z) z[[j]]["lower", ]))
      if(any(is.na(finalCoefs0[[j]][[4]]))){finalCoefs0[[j]][[4]][is.na(finalCoefs0[[j]][[4]])] <- -Inf}
      finalCoefs0[[j]][[5]] <- t(sapply(mod0coefs, function(z) z[[j]]["upper", ]))
      if(any(is.na(finalCoefs0[[j]][[5]]))){finalCoefs0[[j]][[5]][is.na(finalCoefs0[[j]][[5]])] <- Inf}
      finalCoefs0[[j]][[6]] <- unname((N - n - mod0sizes - 1)[j, ])
      names(finalCoefs0[[j]]) <- c("b", "se", "P", "lower", "upper", "df.res")
    }
    names(finalCoefs0) <- names(fit0$mods)
  }
  gammas <- seq(ceiling(alpha * niter)/niter, 1 - 1/niter, by = 1/niter)
  penalty <- ifelse(length(gammas) > 1, (1 - log(min(gammas))), 1)
  adjCIs0 <- list()
  for(i in 1:p){
    k <- ncol(finalCoefs0[[i]]$P)
    pvals0 <- gamma0 <- c()
    for(j in 1:k){
      quantGam0 <- quantile(finalCoefs0[[i]][["P"]][, j], gammas)/gammas
      pvals0[j] <- pmin((min(quantGam0) * penalty), 1)
      gamma0[j] <- gammas[which.min(quantGam0)]
    }
    sInd <- rep(1:k, each = niter)
    adjCIs0[[i]] <- suppressMessages(
      t(mapply(adjustedCI, lci = split(finalCoefs0[[i]][["lower"]], sInd), 
               rci = split(finalCoefs0[[i]][["upper"]], sInd), 
               centers = split(finalCoefs0[[i]][["b"]], sInd),
               ses = split(finalCoefs0[[i]][["se"]], sInd),
               df.res = list(df.res = finalCoefs0[[i]][["df.res"]]),
               gamma.min = min(gammas), ci.level = (1 - alpha), var = 1:k))
    )
    colnames(adjCIs0[[i]]) <- c("lower", "upper")
    rownames(adjCIs0[[i]]) <- rownames(mod0Freqs[[i]])
    err1 <- unname(apply(adjCIs0[[i]], 1, function(z) all(z == 0)))
    if(any(err1)){
      message(paste0("Errors in ", sum(err1), " adjusted CI", ifelse(sum(err1) > 1, "s", ""), "\n"))
      err2 <- which(err1)
      for(e in 1:length(err2)){
        centers <- finalCoefs0[[i]]$b[, err2[e]]
        centers <- centers[!is.infinite(centers)]
        margerrs <- sd(centers) * qt(c(alpha/2, 1 - alpha/2), length(centers) - 1)
        adjCIs0[[i]][err2[e], ] <- mean(centers) + margerrs
      }
    }
    adjCIs0[[i]] <- data.frame(adjCIs0[[i]])
    adjCIs0[[i]][adjCIs0[[i]]$lower == -Inf, ] <- NA
    selectp <- ifelse(is.na(pvals0), FALSE, ifelse(pvals0 <= alpha, TRUE, FALSE))
    select_ci <- ifelse(
      is.na(adjCIs0[[i]]$lower), FALSE, ifelse(
        adjCIs0[[i]]$lower <= 0 & adjCIs0[[i]]$upper >= 0, FALSE, TRUE))
    adjCIs0[[i]] <- data.frame(predictor = factor(rownames(adjCIs0[[i]])), 
                               select_ci = select_ci, lower = adjCIs0[[i]]$lower, 
                               b = rowMeans(adjCIs0[[i]]), upper = adjCIs0[[i]]$upper, 
                               Pvalue = pvals0, select = selectp, gamma = gamma0, 
                               freq = mod0Freqs[[i]][, 1])
    rownames(adjCIs0[[i]]) <- 1:nrow(adjCIs0[[i]])
    if(any(err1)){attr(adjCIs0[[i]], "err") <- as.character(adjCIs0[[i]]$predictor)[err2]}
  }
  if("err" %in% unlist(lapply(adjCIs0, function(z) names(attributes(z))))){
    attr(adjCIs0, "err") <- which(sapply(lapply(adjCIs0, function(z){
      names(attributes(z))}), function(zz) "err" %in% zz))
  }
  names(adjCIs0) <- names(fit0$mods)
  for(i in 1:p){
    finalCoefs0[[i]]$se[is.infinite(finalCoefs0[[i]]$se)] <- NA
    finalCoefs0[[i]]$lower[is.infinite(finalCoefs0[[i]]$lower)] <- NA
    finalCoefs0[[i]]$upper[is.infinite(finalCoefs0[[i]]$upper)] <- NA
  }
  varMods0 <- lapply(varMods0, function(z) do.call(cbind.data.frame, z))
  ########
  nodeVars <- lapply(lapply(adjCIs0, '[[', "predictor"), as.character)
  selector <- c("select", "select_ci")
  finalFits <- finalCoefs <- vector("list", 2)
  if(verbose){cat("Estimating final models\n")
    pb <- txtProgressBar(max = (p + 1) * 2, char = "-", width = npb)}
  for(j in 1:2){
    finalVars <- lapply(lapply(lapply(adjCIs0, '[[', selector[j]), function(z){
      ifelse(!z, 0, 1)}), function(zz) which(zz == 1))
    varMods <- lapply(1:p, function(z) list(mod0 = nodeVars[[z]][finalVars[[z]]]))
    if(any(sapply(lapply(varMods, '[[', "mod0"), length) == 0)){
      v0 <- which(sapply(lapply(varMods, '[[', "mod0"), length) == 0)
      for(jj in seq_along(v0)){varMods[[v0[jj]]]$mod0 <- "1"}
    }
    names(varMods) <- names(allNames)
    attributes(varMods)[c("criterion", "method", "select")] <- c(criterion, varMethod, selector[j])
    for(i in 1:p){
      if(any(grepl(":", varMods[[i]]$mod0))){
        ints <- varMods[[i]]$mod0[grepl(":", varMods[[i]]$mod0)]
        mains <- unique(unlist(strsplit(ints, ":")))
        mains <- mains[order(match(mains, colnames(data)))]
        varMods[[i]]$mod0 <- union(varMods[[i]]$mod0, union(mains, ints))
        varMods[[i]]$mod0 <- varMods[[i]]$mod0[match(allNames[[i]], varMods[[i]]$mod0)]
        varMods[[i]]$mod0 <- varMods[[i]]$mod0[!is.na(varMods[[i]]$mod0)]
      }
      attr(varMods[[i]], "family") <- ifelse(
        length(type) == 1, ifelse(type == "gaussian", "g", "c"), 
        ifelse(type[i] %in% c("g", "gaussian"), "g", "c")
      )
    }
    finalFits[[j]] <- fitNetwork(data = data, moderators = m, type = varMods, rule = rule, 
                                 threshold = threshold, gamma = gamma, center = center, 
                                 exogenous = exogenous, verbose = FALSE) #, ...)
    if(verbose){setTxtProgressBar(pb, ifelse(j == 1, 1, 2 + p))}
    finalCoefs[[j]] <- suppressWarnings(suppressMessages(lapply(1:p, function(z){
      coefs1 <- summary(finalFits[[j]]$fitobj[[z]])$coef[-1, , drop = FALSE]
      if(nrow(coefs1) == 0){return(NA)}
      coefs2 <- confint(finalFits[[j]]$fitobj[[z]], level = 1 - alpha)[-1, , drop = FALSE]
      coefs <- cbind.data.frame(coefs1, coefs2)
      colnames(coefs) <- c("b", "se", "t", "P", "lower", "upper")
      coefs <- coefs[match(allNames[[z]], rownames(coefs)), ]
      coefs <- data.frame(predictor = factor(allNames[[z]]), lower = coefs$lower, 
                          b = coefs$b, upper = coefs$upper, Pvalue = coefs$P,
                          select = ifelse(is.na(coefs$b), FALSE, TRUE))
      rownames(coefs) <- 1:nrow(coefs)
      if(any(is.na(coefs))){coefs[is.na(coefs)] <- NA}
      if(verbose){setTxtProgressBar(pb, ifelse(j == 1, 1 + z, 2 + p + z))}
      return(coefs)
    })))
    names(finalCoefs[[j]]) <- colnames(data)[1:p]
  }
  names(finalFits) <- c("fit0", "fit.5")
  names(finalCoefs) <- c("fitCIs", "fitCIs.5")
  if(verbose){close(pb)}
  ########
  boots <- lapply(1:niter, function(z){
    iter_p1 <- list(vars = varMods0[[z]], fit = fits[[z]])
    iter_p2 <- rev(list(sample = list(data = samps[[z]], inds = sampInd[[z]]), varMods = vars[[z]]))
    if(!saveMods){iter_p1 <- iter_p1[-2]; iter_p2 <- iter_p2[-1]}
    return(append(iter_p1, iter_p2))
  })
  names(boots) <- paste0("iter", 1:niter)
  out <- list(call = preout, fits = append(finalFits, append(list(fit1 = fit0), finalCoefs)))
  out$modLLs <- list(omnibus = modLL(out$fits[1:3]), 
                     modLL1 = modLL(out$fits[1:2], all = TRUE),
                     modLL2 = modLL(out$fits[2:3], all = TRUE))
  out$adjCIs <- adjCIs0
  out$samples <- list(coefs = finalCoefs0, iters = boots)
  if(verbose){cat("\n"); print(Sys.time() - t1)}
  out
}
