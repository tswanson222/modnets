##### resample_parallel: bootstrapping or multi-sample splits for variable selection
resample_parallel <- function(data, m = NULL, niter = 10, sampMethod = "split", criterion = "AIC",
                              method = "glmnet", rule = "OR", gamma = .5, nfolds = 10, 
                              nlam = 50, which.lam = "min", threshold = FALSE, bonf = FALSE, 
                              alpha = .05, exogenous = TRUE, split = .5, center = TRUE, 
                              scale = FALSE, varSeed = NULL, seed = NULL, verbose = TRUE, 
                              lags = NULL, binary = NULL, type = 'g', saveMods = TRUE,
                              saveData = FALSE, saveVars = FALSE, fitit = TRUE, 
                              nCores = 1, cluster = 'mclapply', ...){
  t1 <- Sys.time() # RESAMPLE START
  if(!is.null(lags)){
    if(!is.null(m)){
      if(all(m == 0)){
        m <- NULL
        if(!method %in% c('glmnet', 'subset')){method <- 'subset'}
      }
    }
    lags <- 1
  }
  if(!is.null(m)){method <- 'glinternet'}
  N <- nrow(data) - ifelse(is.null(lags), 0, lags)
  data <- data.frame(data)
  criterion <- toupper(match.arg(tolower(criterion), c(
    "cv", "aic", "bic", "ebic", "cp", "rss", "adjr2", "rsq", "r2")))
  if(!criterion %in% c("AIC", "BIC", "EBIC", "CV")){method <- "subset"}
  method <- match.arg(method, c("glinternet", "hiernet", "subset", "glmnet"))
  if(is.character(m)){
    m <- ifelse(all(m == "all"), list(1:ncol(data)), list(which(colnames(data) %in% m)))[[1]]
  }
  if(!is.null(binary)){
    type <- rep("g", ifelse(is.null(m) | !exogenous, ncol(data), ncol(data) - 1))
    type[intersect(seq_along(type), binary)] <- "c"
  }
  if(method == "glinternet" & is.null(m)){m <- 1:ncol(data)}
  if(!is.null(m)){if(length(m) > 1){exogenous <- FALSE}}
  if(!is.character(sampMethod)){sampMethod <- c("split", "bootstrap", "stability")[sampMethod]}
  sampMethod <- match.arg(tolower(sampMethod), c("split", "bootstrap", "stability"))
  preout <- list(niter = niter, criterion = criterion, sampMethod = sampMethod,
                 method = method, moderators = m, rule = rule, alpha = alpha, 
                 bonf = bonf, gamma = gamma, nfolds = nfolds, which.lam = which.lam,
                 split = split, center = center, scale = scale, exogenous = exogenous, 
                 varSeed = varSeed, type = type, lags = lags)
  args <- tryCatch({list(...)}, error = function(e){list()}) ##### RESAMPLE ARGS
  if(length(args) != 0){preout <- append(preout, args)}
  if(is.null(lags)){preout$lags <- NULL} else {preout$rule <- NULL}
  if(criterion != "CV"){preout[c("nfolds", "which.lam")] <- NULL; which.lam <- "min"}
  if(nlam != 50 & method != "subset"){preout$nlam <- nlam}
  if(sampMethod == "bootstrap"){preout$split <- NULL}
  if(criterion != "EBIC"){preout$gamma <- NULL}
  if(is.null(varSeed)){preout$varSeed <- NULL}
  lam <- ifelse(grepl("min", which.lam), 1, 2)
  if(length(seed) <= 1){
    if(length(seed) == 1){set.seed(seed)}
    seeds <- sample(1:10000000, niter, replace = FALSE)
  } else {
    seeds <- seed
  }
  sampInd <- samps <- vars <- vars1 <- fits <- train <- test <- list()
  if(nCores > 1 | isTRUE(nCores)){
    if(isTRUE(nCores)){nCores <- parallel::detectCores()}
    if(grepl('Windows', sessionInfo()$running)){cluster <- 'SOCK'}
    if(sampMethod != 'bootstrap'){
      if(split <= 0 | split >= 1){stop('split size must be between 0 and 1')}
      n <- floor(N * split)
    } else {
      n <- N
    }
    for(i in 1:niter){
      set.seed(seeds[i])
      sampInd[[i]] <- sample(1:N, n, replace = (sampMethod == 'bootstrap'))
      if(sampMethod != 'bootstrap'){
        if(is.null(lags)){
          train[[i]] <- data[sampInd[[i]], ]
          test[[i]] <- data[-sampInd[[i]], ]
        } else {
          train[[i]] <- test[[i]] <- data
          attr(train[[i]], 'samp_ind') <- sampInd[[i]]
          attr(test[[i]], 'samp_ind') <- (1:N)[-sampInd[[i]]]
        }
        samps[[i]] <- vector('list', length = 2)
        samps[[i]][[1]] <- train[[i]]
        samps[[i]][[2]] <- test[[i]]
        names(samps[[i]]) <- c('train', 'test')
      } else {
        if(is.null(lags)){
          samps[[i]] <- data[sampInd[[i]], ]
        } else {
          samps[[i]] <- data
          attr(samps[[i]], 'samp_ind') <- sampInd[[i]]
        }
      }
    }
    staticArgs <- list(m = m, criterion = criterion, center = center, method = method, 
                       nfolds = nfolds, gamma = gamma, lags = lags, nlam = nlam, 
                       type = type, varSeed = varSeed, exogenous = exogenous, 
                       scale = scale, verbose = FALSE, saveMods = saveMods, 
                       rule = rule, which.lam = which.lam)
    sampFun <- function(samp, seed, sampMethod, args){
      set.seed(seed)
      criterion <- args$criterion
      saveMods <- args$saveMods
      vars <- vars1 <- fits <- list()
      if(sampMethod != 'bootstrap'){
        vargs1 <- append(list(data = samp$train), args[intersect(names(args), formalArgs('varSelect'))])
        vars <- do.call(varSelect, vargs1)
        if(!saveMods){
          for(j in seq_len(length(vars))){
            if(criterion != 'CV'){vars[[j]]$fitobj$fit$fitted <- NA}
            if(criterion == 'CV'){
              vars[[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
              vars[[j]]$fitobj$fit0 <- vars[[j]]$fitobj$fit1se <- NA
            }
          }
        }
        if(sampMethod == 'stability'){
          vargs2 <- replace(vargs1, 'data', list(data = samp$test))
          vars1 <- do.call(varSelect, vargs2)
          if(!saveMods){
            for(j in seq_len(length(vars1))){
              if(criterion != 'CV'){vars1[[j]]$fitobj$fit$fitted <- NA}
              if(criterion == 'CV'){
                vars1[[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
                vars1[[j]]$fitobj$fit0 <- vars1[[j]]$fitobj$fit1se <- NA
              }
            }
          }
        } else {
          fitargs <- append(list(data = samp$test, moderators = args$m, type = vars, threshold = FALSE), 
                            args[setdiff(names(args), c('saveMods', 'm', 'type'))])
          fitargs <- fitargs[intersect(names(fitargs), formalArgs('fitNetwork'))]
          fits <- do.call(fitNetwork, fitargs)
          if(!saveMods){fits$mods0 <- NULL}
        }
      } else {
        vargs1 <- append(list(data = samp), args[intersect(names(args), formalArgs('varSelect'))])
        vars <- do.call(varSelect, vargs1)
        if(!saveMods){
          for(j in seq_len(length(vars))){
            if(criterion != 'CV'){vars[[j]]$fitobj$fit$fitted <- NA}
            if(criterion == 'CV'){
              vars[[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
              vars[[j]]$fitobj$fit0 <- vars[[j]]$fitobj$fit1se <- NA
            }
          }
        }
        fitargs <- append(list(data = samp, moderators = args$m, type = vars, threshold = FALSE), 
                          args[setdiff(names(args), c('saveMods', 'm', 'type'))])
        fitargs <- fitargs[intersect(names(fitargs), formalArgs('fitNetwork'))]
        fits <- do.call(fitNetwork, fitargs)
        if(!saveMods){fits$mods0 <- NULL}
      }
      out <- list(vars = vars, vars1 = vars1, fits = fits)
      return(out)
    }
    obj0 <- c('sampFun', 'varSelect', 'fitNetwork', 'Matrix', 'net', 'netInts')
    obj1 <- switch(2 - is.null(lags), c('nodewise', 'modNet', 'modLL'), 
                   c('lagMat', 'SURfit', 'SURnet', 'SURll', 'surEqs', 'getCoefs', 'systemfit'))
    obj2 <- switch(2 - (method == 'glinternet'), c('fitHierLASSO', ifelse(criterion == 'CV', 'glinternet.cv', 'glinternet')),
                   ifelse(method == 'glmnet', ifelse(criterion == 'CV', 'cv.glmnet', 'glmnet'), 'regsubsets'))
    if(method == 'glmnet'){obj2 <- c(obj2, 'lassoSelect')}
    objects <- c(obj0, obj1, obj2)
    if(tolower(cluster) != 'mclapply'){
      cluster <- match.arg(toupper(cluster), c('SOCK', 'FORK'))
      cl <- parallel::makeCluster(nCores, type = cluster)
      if(cluster == 'SOCK'){parallel::clusterExport(cl, objects, envir = environment())}
    } else {
      cl <- nCores
    }
    if(verbose){
      pbapply::pboptions(type = 'timer', char = '-')
      cl_out <- pbapply::pblapply(1:niter, function(i){
        sampFun(samp = samps[[i]], seed = seeds[i], 
                sampMethod = sampMethod, args = staticArgs)
      }, cl = cl)
    } else if(tolower(cluster) == 'mclapply'){
      cl_out <- parallel::mclapply(1:niter, function(i){
        sampFun(samp = samps[[i]], seed = seeds[i], 
                sampMethod = sampMethod, args = staticArgs)
      }, mc.cores = nCores)
    } else {
      cl_out <- parallel::parLapply(cl, 1:niter, function(i){
        sampFun(samp = samps[[i]], seed = seeds[i], 
                sampMethod = sampMethod, args = staticArgs)
      })
    }
    if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
    vars <- lapply(cl_out, '[[', 1)
    vars1 <- lapply(cl_out, '[[', 2)
    fits <- lapply(cl_out, '[[', 3)
    rm(cl_out, staticArgs, objects, obj1, obj2, obj0, sampFun, cl)
  } else {
    if(sampMethod != "bootstrap"){
      if(split <= 0 | split >= 1){stop("split size must be between 0 and 1")}
      train <- list(); test <- list()
      n <- floor(N * split)
      for(i in 1:niter){
        set.seed(seeds[i]) # RESAMPLE SPLIT
        sampInd[[i]] <- sample(1:N, n, replace = FALSE)
        if(is.null(lags)){
          train[[i]] <- data[sampInd[[i]], ]
          test[[i]] <- data[-sampInd[[i]], ]
        } else {
          train[[i]] <- test[[i]] <- data
          attr(train[[i]], "samp_ind") <- sampInd[[i]]
          attr(test[[i]], "samp_ind") <- (1:N)[-sampInd[[i]]]
        }
        samps[[i]] <- vector("list", length = 2)
        samps[[i]][[1]] <- train[[i]]
        samps[[i]][[2]] <- test[[i]]
        names(samps[[i]]) <- c("train", "test")
        if(verbose){
          if(i == 1){cat("\n")}
          t2 <- Sys.time()
          cat("************ Sample Split: ", i, "/", niter, " ************\n", sep = "")
        }
        set.seed(seeds[i])
        vars[[i]] <- varSelect(data = train[[i]], m = m, criterion = criterion, center = center,
                               method = method, nfolds = nfolds, gamma = gamma, lags = lags,
                               nlam = nlam, type = type, varSeed = varSeed, exogenous = exogenous,
                               scale = scale, verbose = ifelse(sampMethod == "stability", ifelse(
                                 verbose, "pbar", FALSE), verbose))
        if(!saveMods){
          for(j in seq_len(length(vars[[i]]))){
            if(criterion != "CV"){vars[[i]][[j]]$fitobj$fit$fitted <- NA}
            if(criterion == "CV"){
              vars[[i]][[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
              vars[[i]][[j]]$fitobj$fit0 <- vars[[i]][[j]]$fitobj$fit1se <- NA
            }
          }
        } # RESAMPLE STABILITY
        if(sampMethod == "stability"){
          vars1[[i]] <- varSelect(data = test[[i]], m = m, criterion = criterion, center = center,
                                  method = method, nfolds = nfolds, gamma = gamma, lags = lags,
                                  nlam = nlam, type = type, varSeed = varSeed, exogenous = exogenous, 
                                  scale = scale, verbose = ifelse(verbose, "pbar", FALSE))
          if(!saveMods){
            for(j in seq_len(length(vars1[[i]]))){
              if(criterion != "CV"){vars1[[i]][[j]]$fitobj$fit$fitted <- NA}
              if(criterion == "CV"){
                vars1[[i]][[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
                vars1[[i]][[j]]$fitobj$fit0 <- vars1[[i]][[j]]$fitobj$fit1se <- NA
              }
            }
          }
          if(verbose & !method %in% c("subset", "glmnet")){
            t3 <- Sys.time() - t2
            cat(paste0("Time: ", round(t3, 2), " ", attr(t3, "units")), "\n\n")
          } # RESAMPLE SPLIT FIT
        } else {
          fits[[i]] <- fitNetwork(data = test[[i]], moderators = m, type = vars[[i]], 
                                  threshold = FALSE, which.lam = which.lam, lags = lags,
                                  gamma = gamma, center = center, rule = rule, scale = scale,
                                  exogenous = exogenous, verbose = FALSE, ...)
          if(!saveMods){fits[[i]]$mods0 <- NULL}
        }
      }
    }
    if(sampMethod == "bootstrap"){
      for(i in 1:niter){
        set.seed(seeds[i]) # RESAMPLE BOOT
        sampInd[[i]] <- sample(1:N, N, replace = TRUE)
        if(is.null(lags)){
          samps[[i]] <- data[sampInd[[i]], ]
        } else {
          samps[[i]] <- data
          attr(samps[[i]], "samp_ind") <- sampInd[[i]]
        }
        if(verbose == TRUE){
          if(i == 1){cat("\n")}
          cat("************* Bootstrap: ", i, "/", niter, " **************\n", sep = "")
        }
        set.seed(seeds[i])
        vars[[i]] <- varSelect(data = samps[[i]], m = m, criterion = criterion, 
                               center = center, method = method, gamma = gamma, 
                               lags = lags, type = type, varSeed = varSeed, 
                               scale = scale, exogenous = exogenous,
                               verbose = verbose)
        if(!saveMods){
          for(j in seq_len(length(vars[[i]]))){
            if(criterion != "CV"){vars[[i]][[j]]$fitobj$fit$fitted <- NA}
            if(criterion == "CV"){
              vars[[i]][[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
              vars[[i]][[j]]$fitobj$fit0 <- vars[[i]][[j]]$fitobj$fit1se <- NA
            }
          }
        }
        fits[[i]] <- fitNetwork(data = samps[[i]], moderators = m, type = vars[[i]],
                                threshold = FALSE, which.lam = which.lam, lags = lags,
                                gamma = gamma, center = center, rule = rule, scale = scale,
                                exogenous = exogenous, verbose = FALSE, ...)
        if(!saveMods){fits[[i]]$mods0 <- NULL}
      }
    }
  }
  fit0 <- fitNetwork(data = data, moderators = m, type = type, rule = rule, scale = scale,
                     threshold = FALSE, gamma = gamma, center = center, lags = lags,
                     exogenous = exogenous, verbose = FALSE, ...)
  if(!is.null(lags)){
    if(sampMethod != "stability"){
      fits <- lapply(fits, function(fi){
        append(fi$SURnet[setdiff(names(fi$SURnet), 'data')], append(
          switch(2 - ('SURll' %in% names(fi)), list(SURll = fi$SURll), list()), 
          list(data = fi$SURnet$data, fitobj = fi$SURfit$eq)))
      })
    }
    fit0 <- append(fit0$SURnet[setdiff(names(fit0$SURnet), 'data')], append(
      switch(2 - ('SURll' %in% names(fit0)), list(SURll = fit0$SURll), list()), 
      list(data = fit0$SURnet$data, fitobj = fit0$SURfit$eq)))
  }
  p <- length(fit0$mods)
  vs <- names(fit0$mods)
  allNames <- lapply(fit0$mods, function(z) rownames(z$model)[-1])
  ##### RESAMPLE VARMODS #####
  varMods0 <- lapply(vars, function(z) lapply(z, function(zz) zz[[lam]]))
  if(sampMethod == "stability"){
    if(!is.numeric(threshold)){threshold <- .6}
    varMods1 <- lapply(vars1, function(z) lapply(z, '[[', lam))
    simultvars <- lapply(seq_len(niter), function(z){
      lapply(seq_len(p), function(zz){
        varSamp1 <- !is.na(match(allNames[[zz]], varMods0[[z]][[zz]]))
        varSamp2 <- !is.na(match(allNames[[zz]], varMods1[[z]][[zz]]))
        return(cbind(varSamp1, varSamp2, varSamp1 & varSamp2))
      })
    })
    varFreqs <- suppressMessages(suppressWarnings(lapply(seq_len(p), function(z){
      simultSel <- t(sapply(seq_len(niter), function(zz) simultvars[[zz]][[z]][, 3]))
      s1sel <- t(sapply(seq_len(niter), function(zz) simultvars[[zz]][[z]][, 1]))
      s2sel <- t(sapply(seq_len(niter), function(zz) simultvars[[zz]][[z]][, 2]))
      colnames(simultSel) <- colnames(s1sel) <- colnames(s2sel) <- allNames[[z]]
      freqs <- colMeans(simultSel)
      s1freqs <- colMeans(s1sel)
      s2freqs <- colMeans(s2sel)
      coefs1 <- summary(fit0$fitobj[[z]])$coefficients[-1, , drop = FALSE] 
      coefs2 <- confint(fit0$fitobj[[z]], level = 1 - alpha)[-1, , drop = FALSE]
      outFreqs <- data.frame(predictor = factor(allNames[[z]]), select = freqs >= threshold,
                             lower = coefs2[, 1], b = coefs1[, 1], upper = coefs2[, 2], 
                             Pvalue = coefs1[, 4], freq = freqs, split1 = s1freqs, split2 = s2freqs)
      rownames(outFreqs) <- 1:nrow(outFreqs)
      return(outFreqs)
    })))
    preout$nlam <- nlam
    attributes(varFreqs) <- list(
      names = vs, type = type, criterion = criterion, rule = rule, gamma = gamma, 
      center = center, scale = scale, exogenous = exogenous, moderators = m, 
      threshold = FALSE, lags = lags, method = unique(sapply(vars, attr, "method")))
    if(is.null(lags)){attr(varFreqs, "lags") <- NULL}
    split1 <- lapply(seq_len(niter), function(z){list(vars = varMods0[[z]], varMods = vars[[z]])})
    split2 <- lapply(seq_len(niter), function(z){list(vars = varMods1[[z]], varMods = vars1[[z]])})
    names(split1) <- names(split2) <- paste0("iter", 1:niter)
    out <- list(call = preout, freqs = varFreqs, split1 = split1, split2 = split2)
    if(!(method == "subset" & criterion == "CV")){
      out$stability <- suppressWarnings(tryCatch({
        append(stability(out), list(seeds = seeds))}, error = function(e){
          list(split1 = split1, split2 = split2, seeds = seeds)}))
    } else {
      out$seeds <- seeds
    }
    if(!saveVars & criterion != "CV"){out <- out[-c(3:4)]}
    if(fitit){
      tryfit <- tryCatch({modSelect(obj = out, data = data, fit = TRUE, saveMods = FALSE)},
                         error = function(e){TRUE})
      if(!isTRUE(tryfit)){out$fit0 <- tryfit}
    }
    out$results <- structure(cbind.data.frame(
      outcome = rep(gsub('[.]y$', '', vs), sapply(out$freqs, nrow)), 
      do.call(rbind, out$freqs)), row.names = 1:sum(sapply(out$freqs, nrow)))
    if(saveData){out$data <- data}
    if(verbose){cat("\n"); print(Sys.time() - t1)}
    return(out)
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
    if(nCores > 1 & FALSE){
      summaryFun <- function(i, fits, p, alpha, allNames, bonf, mod0sizes){
        mod0coefs <- suppressWarnings(suppressMessages(
          lapply(fits[[i]]$fitobj, function(z){
            coefs1 <- t(summary(z)$coefficients[-1, , drop = FALSE])
            if(length(coefs1) != 0){
              coefs2 <- t(confint(z, level = 1 - alpha)[-1, , drop = FALSE])
            } else {
              vsx <- allNames[[as.character(z$terms[[2]])]]
              coefs1 <- matrix(NA, nrow = 4, ncol = length(vsx))
              coefs2 <- matrix(NA, nrow = 2, ncol = length(vsx))
              colnames(coefs1) <- colnames(coefs2) <- vsx
            }
            return(rbind(coefs1, coefs2))
          })))
        for(j in 1:p){
          rownames(mod0coefs[[j]]) <- c("b", "se", "t", "P", "lower", "upper")
          mod0coefs[[j]] <- mod0coefs[[j]][, match(allNames[[j]], colnames(mod0coefs[[j]]))]
          colnames(mod0coefs[[j]]) <- allNames[[j]]
          if(bonf){mod0coefs[[j]]["P", ] <- pmin(mod0coefs[[j]]["P", ] * mod0sizes[j, i], 1)}
        }
        return(mod0coefs)
      }
      cl <- parallel::makeCluster(nCores, type = 'SOCK')
      parallel::clusterExport(cl, 'summaryFun', envir = environment())
      mod0coefs <- parallel::parLapply(cl, 1:niter, summaryFun, fits = fits, p = p, 
                                       alpha = alpha, allNames = allNames, 
                                       bonf = bonf, mod0sizes = mod0sizes)
      parallel::stopCluster(cl)
    } else {
      if(verbose){
        npb <- seq(43, 53, by = 2)[sum(niter >= c(10, 100, 1000, 10000, 100000)) + 1]
        cat("\n"); cat(paste0(rep("#", npb), collapse = ""), "\n")
        cat(paste0("Estimating ", (1 - alpha) * 100, "% CIs\n"))
        pb <- txtProgressBar(min = 0, max = niter, style = 1, char = "-", width = npb)
      }
      for(i in 1:niter){
        mod0coefs[[i]] <- suppressWarnings(suppressMessages(
          lapply(fits[[i]]$fitobj, function(z){
            coefs1 <- t(summary(z)$coefficients[-1, , drop = FALSE])
            if(length(coefs1) != 0){
              coefs2 <- t(confint(z, level = 1 - alpha)[-1, , drop = FALSE])
            } else {
              vsx <- allNames[[as.character(z$terms[[2]])]]
              coefs1 <- matrix(NA, nrow = 4, ncol = length(vsx))
              coefs2 <- matrix(NA, nrow = 2, ncol = length(vsx))
              colnames(coefs1) <- colnames(coefs2) <- vsx
            }
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
          margerrs <- sd(centers, na.rm = TRUE) * qt(c(alpha/2, 1 - alpha/2), length(centers) - 1)
          adjCIs0[[i]][err2[e], ] <- mean(centers, na.rm = TRUE) + margerrs
        }
      }
      adjCIs0[[i]] <- data.frame(adjCIs0[[i]])
      if(!all(is.na(adjCIs0[[i]]))){
        if(any(adjCIs0[[i]][!is.na(adjCIs0[[i]]$lower), "lower"] == -Inf)){
          adjCIs0[[i]][which(adjCIs0[[i]]$lower == -Inf), ] <- NA
        }
      }
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
    attributes(adjCIs0) <- list(
      names = vs, type = type, criterion = criterion, rule = rule, gamma = gamma, 
      center = center, scale = scale, exogenous = exogenous, moderators = m, 
      threshold = FALSE, lags = lags, method = unique(sapply(vars, attr, "method")))
    if(is.null(lags)){attr(adjCIs0, "lags") <- NULL}
    if("err" %in% unlist(lapply(adjCIs0, function(z) names(attributes(z))))){
      attr(adjCIs0, "err") <- names(which(sapply(lapply(adjCIs0, function(z){
        names(attributes(z))}), function(zz) "err" %in% zz)))
    }
    for(i in 1:p){
      finalCoefs0[[i]]$se[is.infinite(finalCoefs0[[i]]$se)] <- NA
      finalCoefs0[[i]]$lower[is.infinite(finalCoefs0[[i]]$lower)] <- NA
      finalCoefs0[[i]]$upper[is.infinite(finalCoefs0[[i]]$upper)] <- NA
    }
  }
  varMods0 <- lapply(varMods0, function(z) do.call(cbind.data.frame, z))
  boots <- lapply(1:niter, function(z){
    p1 <- list(vars = varMods0[[z]])
    if(saveMods){
      fits[[z]]$fitobj <- NULL
      if(is.null(lags)){fits[[z]]$mods0$models <- NULL}
      if(!is.null(lags)){attr(fits[[z]], 'SURnet') <- TRUE}
      if(!saveData){fits[[z]]$data <- NULL}
      p1$fit <- fits[[z]]
    }
    if(saveVars){p1$varMods <- vars[[z]]}
    return(append(p1, list(samp_inds = sampInd[[z]])))
  })
  names(boots) <- paste0("iter", 1:niter)
  out <- list(call = preout, adjCIs = adjCIs0)
  out$samples <- list(coefs = finalCoefs0, iters = boots, seeds = seeds)
  if(fitit){
    tryfit <- tryCatch({modSelect(obj = out, data = data, fit = TRUE, saveMods = FALSE)},
                       error = function(e){TRUE})
    if(!isTRUE(tryfit)){out$fit0 <- tryfit}
  }
  out$results <- structure(cbind.data.frame(
    outcome = rep(gsub('[.]y$', '', vs), sapply(out$adjCIs, nrow)), 
    do.call(rbind, out$adjCIs)), row.names = 1:sum(sapply(out$adjCIs, nrow)))
  out$data <- data
  if(verbose){
    if(nCores == 1){cat("\n")}
    print(Sys.time() - t1)
  }
  attr(out, "time") <- Sys.time() - t1
  attr(out, 'resample') <- TRUE
  out
}
