resample2 <- function(data, niter = 10, m = NULL, criterion = "CV", sampMethod = "multisplit", 
                      gamma = .25, nfolds = 10, nlam = 50, which.lam = "min", rule = "AND",
                      varMethod = "glinternet", threshold = FALSE, bonf = TRUE, alpha = .05,
                      split = .5, center = TRUE, varSeed = NULL, seeds = NULL, verbose = TRUE,
                      exogenous = TRUE, binary = NULL, saveMods = FALSE, type = "gaussian", ...){
  t1 <- Sys.time()
  N <- nrow(data)
  criterion <- toupper(match.arg(
    tolower(criterion), 
    choices = c("cv", "aic", "bic", "ebic", "cp", "rss", "adjr2", "rsq", "r2")))
  if(criterion != "CV"){which.lam <- "min"}
  if(is.null(seeds)){seeds <- 1:niter}
  if(!is.null(binary)){
    type <- rep("g", ifelse(is.null(m) | !exogenous, ncol(data), ncol(data) - 1))
    type[intersect(seq_along(type), binary)] <- "c"
  }
  if(varMethod == "glinternet" & is.null(m)){
    exogenous <- FALSE
    m <- 1:ncol(data)
  } else if(!is.null(m)){
    if(length(m) > 1){exogenous <- FALSE}
  }
  sampMethod <- match.arg(tolower(sampMethod), c("multisplit", "split", "bootstrap"))
  preout <- list(niter = niter, m = m, criterion = criterion, sampMethod = sampMethod,
                 gamma = gamma, nfolds = nfolds, nlam = nlam, which.lam = which.lam,
                 rule = rule, varMethod = varMethod, bonf = bonf, alpha = alpha, 
                 split = split, center = center, exogenous = exogenous,
                 varSeed = varSeed, seeds = seeds, type = type)
  sampInd <- list(); vars <- list(); fits <- list(); samps <- list()
  if(sampMethod %in% c("multisplit", "split")){
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
      if(verbose == TRUE){
        if(i == 1){cat("\n")}
        cat("************ Sample Split: ", i, "/", niter, " ************\n", sep = "")
      }
      set.seed(seeds[i])
      vars[[i]] <- varSelect(dat = train[[i]], m = m, criterion = criterion, 
                             method = varMethod, nfolds = nfolds, gamma = gamma, 
                             nlam = nlam, type = type, verbose = verbose, 
                             seed = varSeed, all = !exogenous)
      fits[[i]] <- fitNetwork(data = test[[i]], moderators = m, type = vars[[i]], 
                              threshold = threshold, which.lam = which.lam, 
                              gamma = gamma, center = center, rule = rule, 
                              exogenous = exogenous, verbose = FALSE, ...)
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
                              threshold = threshold, which.lam = which.lam,
                              gamma = gamma, center = center, rule = rule, 
                              exogenous = exogenous, verbose = FALSE, ...)
    }
  }
  fit0 <- fitNetwork(data = data, moderators = m, type = type, threshold = threshold, 
                     rule = rule, gamma = gamma, center = center, 
                     exogenous = exogenous, verbose = FALSE, ...)
  allNames <- lapply(fit0$mods, function(z) rownames(z$model)[-1])
  p <- length(allNames)
  varsMod0 <- lapply(vars, function(z) lapply(z, function(zz) zz$mod0))
  mod0Freqs <- lapply(1:p, function(z) data.frame(table(unlist(lapply(varsMod0, function(zz) zz[[z]])))))
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
  mod0sizes <- sapply(varsMod0, function(z) sapply(z, length))
  for(i in 1:niter){
    if(any(mod0sizes[, i] != max(mod0sizes[, i]))){
      ms0 <- which(mod0sizes[, i] != max(mod0sizes[, i]))
      for(j in seq_along(ms0)){
        varsMod0[[i]][[ms0[j]]] <- c(varsMod0[[i]][[ms0[j]]], rep("", max(mod0sizes[, i]) - length(varsMod0[[i]][[ms0[j]]])))
      }
    }
  }
  mod0coefs <- finalCoefs0 <- list()
  if(verbose){
    npb <- seq(43, 53, by = 2)[sum(niter >= c(10, 100, 1000, 10000, 100000)) + 1]
    cat(paste0("\nEstimating adjusted ", (1 - alpha) * 100, "% CIs:\n"))
    pb <- txtProgressBar(min = 0, max = niter, style = 1, char = "-", width = npb)
  }
  for(i in 1:niter){
    mod0coefs[[i]] <- lapply(fits[[i]]$fitobj, function(z){
      coefs1 <- t(summary(z)$coef[-1, , drop = FALSE])
      coefs2 <- t(suppressMessages(confint(z, level = 1 - alpha)[-1, , drop = FALSE]))
      return(rbind(coefs1, coefs2))
    })
    for(j in 1:p){
      rownames(mod0coefs[[i]][[j]]) <- c("b", "se", "t", "P", "lower", "upper")
      mod0coefs[[i]][[j]] <- mod0coefs[[i]][[j]][, match(allNames[[j]], colnames(mod0coefs[[i]][[j]]))]
      colnames(mod0coefs[[i]][[j]]) <- allNames[[j]]
      if(bonf){mod0coefs[[i]][[j]]["P", ] <- pmin(mod0coefs[[i]][[j]]["P", ] * mod0sizes[j, i], 1)}
    }
    if(verbose){
      setTxtProgressBar(pb, i)
      if(i == niter){close(pb)}
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
  adjCIs0 <- NULL
  if(niter > 1){
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
      err1 <- unname(apply(adjCIs0[[i]], 1, function(z) all(z == 0)))
      if(any(err1 == TRUE)){
        message("Errors in some adjusted CIs")
        err2 <- which(err1)
        for(e in 1:length(err2)){
          centers <- finalCoefs0[[i]]$b[, err2[e]]
          centers <- centers[!is.infinite(centers)]
          ses <- finalCoefs0[[i]]$se[, err2[e]]
          ses <- ses[!is.infinite(ses)]
          adjCIs0[[i]][err2[e], ] <- c(mean(centers) - (mean(ses) * qnorm(1 - alpha/2)),
                                       mean(centers) + (mean(ses) * qnorm(1 - alpha/2)))
        }
      }
      colnames(adjCIs0[[i]]) <- c("lower", "upper")
      rownames(adjCIs0[[i]]) <- rownames(mod0Freqs[[i]])
      adjCIs0[[i]] <- data.frame(adjCIs0[[i]])
      adjCIs0[[i]][adjCIs0[[i]]$lower == -Inf, ] <- NA
      adjCIs0[[i]]$b <- rowMeans(adjCIs0[[i]])
      adjCIs0[[i]]$predictor <- factor(rownames(adjCIs0[[i]]))
      adjCIs0[[i]] <- data.frame(predictor = adjCIs0[[i]]$predictor, lower = adjCIs0[[i]]$lower, 
                                 b = adjCIs0[[i]]$b, upper = adjCIs0[[i]]$upper, Pvalue = pvals0, 
                                 gammaMin = gamma0, freq = mod0Freqs[[i]][, 1])
      rownames(adjCIs0[[i]]) <- 1:nrow(adjCIs0[[i]])
      if(any(err1 == TRUE)){adjCIs0[[i]]$err <- err1}
      if(any(rowSums(is.na(adjCIs0[[i]])) != 0)){
        adjCIs0[[i]][rowSums(is.na(adjCIs0[[i]])) != 0, c("Pvalue", "gammaMin")] <- NA
      }
    }
    names(adjCIs0) <- names(fit0$mods)
  }
  for(i in 1:p){
    finalCoefs0[[i]]$se[is.infinite(finalCoefs0[[i]]$se)] <- NA
    finalCoefs0[[i]]$lower[is.infinite(finalCoefs0[[i]]$lower)] <- NA
    finalCoefs0[[i]]$upper[is.infinite(finalCoefs0[[i]]$upper)] <- NA
  }
  varsMod0 <- lapply(varsMod0, function(z) do.call(cbind.data.frame, z))
  nodeVars <- lapply(lapply(adjCIs0, '[[', "predictor"), as.character)
  finalVars <- lapply(lapply(lapply(adjCIs0, '[[', "Pvalue"), function(z) 
    ifelse(is.na(z) | z > alpha, 0, 1)), function(zz) which(zz == 1))
  varsMod1 <- lapply(1:p, function(z) list(mod0 = nodeVars[[z]][finalVars[[z]]]))
  attr(varsMod1, "criterion") <- criterion
  attr(varsMod1, "method") <- varMethod
  for(i in 1:p){
    attr(varsMod1[[i]], "family") <- ifelse(
      length(type) == 1, ifelse(type == "gaussian", "g", "c"), 
      ifelse(type[i] %in% c("g", "gaussian"), "g", "c")
    )
  }
  fit1 <- fitNetwork(data = data, moderators = m, type = varsMod1, threshold = threshold, 
                     rule = rule, gamma = gamma, center = center, exogenous = exogenous, 
                     verbose = FALSE, ...)
  boots <- lapply(1:niter, function(z){
    iter_p1 <- list(vars = varsMod0[[z]], fit = fits[[z]])
    if(!saveMods){iter_p1 <- iter_p1[-2]}
    iter_p2 <- rev(list(sample = list(data = samps[[z]], inds = sampInd[[z]]), varMods = vars[[z]]))
    return(append(iter_p1, iter_p2))
  })
  names(boots) <- paste0("iter", 1:niter)
  out <- list(call = preout, fit1 = fit1, adjCIs = adjCIs0, coefs = finalCoefs0, boots = boots, fit0 = fit0)
  if(verbose){print(Sys.time() - t1)}
  out
}


