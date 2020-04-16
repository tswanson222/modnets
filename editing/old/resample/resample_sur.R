##### resampleOLD: perform either bootstrapping or multi-sample splits for variable selection
resample <- function(dat, niter = 10, criterion = "CV", sampMethod = "multisplit",
                     splitSize = .5, m = 1, varMethod = "glinternet", varSeed = NULL, 
                     alpha = .05, seeds = NULL, verbose = TRUE, bonf = TRUE, ...){
  method <- match.arg(tolower(varMethod), c("glinternet", "hiernet", "subset", 
                                            "exhaustive", "backward", "forward", "seqrep"))
  sampMethod <- match.arg(tolower(sampMethod), c("bootstrap", "multisplit", "split"))
  if(is.null(seeds)){seeds <- 1:niter}
  N <- nrow(dat$full); p <- ncol(dat$Y)
  sampInd <- list(); vars <- list(); fits <- list(); samps <- list()
  if(sampMethod == "bootstrap"){
    for(i in 1:niter){
      set.seed(seeds[i])
      sampInd[[i]] <- sample(1:N, N, replace = TRUE)
      samps[[i]] <- vector("list", length = 3)
      samps[[i]][[1]] <- dat$Y[sampInd[[i]],]
      samps[[i]][[2]] <- dat$X[sampInd[[i]],]
      samps[[i]][[3]] <- dat$full[sampInd[[i]],]
      names(samps[[i]]) <- c("Y", "X", "full")
      if(verbose == TRUE){
        if(i == 1){cat("\n")}
        cat("************* Bootstrap: ", i, "/", niter, " **************\n", sep = "")}
      set.seed(seeds[i])
      vars[[i]] <- varSelect(dat = samps[[i]], method = method, criterion = criterion, 
                             m = m, verbose = verbose, seed = varSeed, ...)
      if(method %in% c("glinternet", "hiernet")){
        fits[[i]] <- vector("list", length = 2)
        fits[[i]][[1]] <- SURfit(dat = samps[[i]], varMods = vars[[i]], mod = "min")
        fits[[i]][[2]] <- SURfit(dat = samps[[i]], varMods = vars[[i]], mod = "1se")
        names(fits[[i]]) <- c("mod0", "mod1se")
      } else {
        fits[[i]] <- SURfit(dat = samps[[i]], varMods = vars[[i]])
      }
    }
  }
  if(sampMethod %in% c("multisplit", "split")){
    if(splitSize <= 0 | splitSize >= 1){stop("splitSize must be between 0 and 1")}
    train <- list(); test <- list()
    n <- floor(N * splitSize)
    for(i in 1:niter){
      set.seed(seeds[i])
      sampInd[[i]] <- sample(1:N, n, replace = FALSE)
      train[[i]] <- vector("list", length = 3)
      train[[i]][[1]] <- dat$Y[sampInd[[i]],]
      train[[i]][[2]] <- dat$X[sampInd[[i]],]
      train[[i]][[3]] <- dat$full[sampInd[[i]],]
      test[[i]] <- vector("list", length = 3)
      test[[i]][[1]] <- dat$Y[-sampInd[[i]],]
      test[[i]][[2]] <- dat$X[-sampInd[[i]],]
      test[[i]][[3]] <- dat$full[-sampInd[[i]],]
      names(train[[i]]) <- names(test[[i]]) <- c("Y", "X", "full")
      samps[[i]] <- vector("list", length = 2)
      samps[[i]][[1]] <- train[[i]]
      samps[[i]][[2]] <- test[[i]]
      names(samps[[i]]) <- c("train", "test")
      if(verbose == TRUE){
        if(i == 1){cat("\n")}
        cat("************ Sample Split: ", i, "/", niter, " ************\n", sep = "")}
      set.seed(seeds[i])
      vars[[i]] <- varSelect(dat = train[[i]], method = method, criterion = criterion, 
                             m = m, verbose = verbose, seed = varSeed, ...)
      if(method %in% c("glinternet", "hiernet")){
        fits[[i]] <- vector("list", length = 2)
        fits[[i]][[1]] <- SURfit(dat = test[[i]], varMods = vars[[i]], mod = "min")
        fits[[i]][[2]] <- SURfit(dat = test[[i]], varMods = vars[[i]], mod = "1se")
        names(fits[[i]]) <- c("mod0", "mod1se")
      } else {
        fits[[i]] <- SURfit(dat = test[[i]], varMods = vars[[i]])
      }
    }
  }
  names(fits) <- names(vars) <- names(samps) <- names(sampInd) <- paste0("iter", 1:niter)
  if(method %in% c("glinternet", "hiernet")){
    data <- dat$X
    if(any(grepl(":", colnames(data)))){data <- data[,!grepl(":", colnames(data))]}
    allNames <- c(colnames(data), apply(combn(colnames(data), 2), 2, paste, collapse = ":"))
    if(!is.null(m)){
      allNames <- c(allNames[!grepl(":", allNames)], 
                    allNames[grepl(":", allNames)][grepl(paste0("X", m), allNames[grepl(":", allNames)])])
    }
    varsMod0 <- lapply(vars, function(iterMods) lapply(iterMods, function(yMod) yMod$mod0))
    varsMod1se <- lapply(vars, function(iterMods) lapply(iterMods, function(yMod) yMod$mod1se))
    mod0Freqs <- lapply(1:p, function(eqs) 
      data.frame(table(unlist(lapply(varsMod0, function(iterMods) iterMods[[eqs]])))))
    mod0Freqs2 <- lapply(1:p, function(freqs) mod0Freqs[[freqs]][,2])
    mod1Freqs <- lapply(1:p, function(eqs) 
      data.frame(table(unlist(lapply(varsMod1se, function(iterMods) iterMods[[eqs]])))))
    mod1Freqs2 <- lapply(1:p, function(freqs) mod1Freqs[[freqs]][,2])
    for(i in 1:p){
      names(mod0Freqs2[[i]]) <- as.character(mod0Freqs[[i]][,1])
      names(mod1Freqs2[[i]]) <- as.character(mod1Freqs[[i]][,1])
      mod0Freqs2[[i]] <- mod0Freqs2[[i]][match(allNames, names(mod0Freqs2[[i]]))]
      names(mod0Freqs2[[i]]) <- allNames
      if(any(is.na(mod0Freqs2[[i]]))){mod0Freqs2[[i]][is.na(mod0Freqs2[[i]])] <- 0}
      mod0Freqs2[[i]] <- mod0Freqs2[[i]]/niter
      mod1Freqs2[[i]] <- mod1Freqs2[[i]][match(allNames, names(mod1Freqs2[[i]]))]
      names(mod1Freqs2[[i]]) <- allNames
      if(any(is.na(mod1Freqs2[[i]]))){mod1Freqs2[[i]][is.na(mod1Freqs2[[i]])] <- 0}
      mod1Freqs2[[i]] <- mod1Freqs2[[i]]/niter
    }
    mod0Freqs <- do.call(cbind, mod0Freqs2)
    mod1Freqs <- do.call(cbind, mod1Freqs2)
    colnames(mod0Freqs) <- colnames(mod1Freqs) <- colnames(dat$Y)
    mod0sizes <- sapply(varsMod0, function(yMod) sapply(yMod, length))
    mod1sizes <- sapply(varsMod1se, function(yMod) sapply(yMod, length))
    for(i in 1:niter){
      if(any(mod0sizes[,i] != max(mod0sizes[,i]))){
        ms0 <- which(mod0sizes[,i] != max(mod0sizes[,i]))
        for(j in seq_along(ms0)){
          varsMod0[[i]][[ms0[j]]] <- c(varsMod0[[i]][[ms0[j]]], rep("", max(mod0sizes[,i]) - length(varsMod0[[i]][[ms0[j]]])))
        }
      }
      if(any(mod1sizes[,i] != max(mod1sizes[,i]))){
        ms1 <- which(mod1sizes[,i] != max(mod1sizes[,i]))
        for(j in seq_along(ms1)){
          varsMod1se[[i]][[ms1[j]]] <- c(varsMod1se[[i]][[ms1[j]]], rep("", max(mod1sizes[,i]) - length(varsMod1se[[i]][[ms1[j]]])))
        }
      }
    }
    if(!is.null(m)){
      mod0coefs <- list(); mod1coefs <- list()
      mod0coefs2 <- list(); mod1coefs2 <- list()
      for(i in 1:niter){
        mod0coefs[[i]] <- t(summary(fits[[i]][[1]])$coefficients[!grepl("Intercept", rownames(summary(fits[[i]][[1]])$coefficients)),])
        mod0coefs2[[i]] <- t(confint(fits[[i]][[1]], level = (1 - alpha))[!grepl("Intercept", rownames(confint(fits[[i]][[1]], level = (1 - alpha)))),])
        mod1coefs[[i]] <- t(summary(fits[[i]][[2]])$coefficients[!grepl("Intercept", rownames(summary(fits[[i]][[2]])$coefficients)),])
        mod1coefs2[[i]] <- t(confint(fits[[i]][[2]], level = (1 - alpha))[!grepl("Intercept", rownames(confint(fits[[i]][[2]], level = (1 - alpha)))),])
        mod0coefs[[i]] <- rbind(mod0coefs[[i]], mod0coefs2[[i]])
        mod1coefs[[i]] <- rbind(mod1coefs[[i]], mod1coefs2[[i]])
        g0 <- matrix(rep(rep(1:p, mod0sizes[,i]), each = 6), nrow = 6)
        g1 <- matrix(rep(rep(1:p, mod1sizes[,i]), each = 6), nrow = 6)
        mod0coefs[[i]] <- lapply(t(split(mod0coefs[[i]], g0)), matrix, nrow = 6)
        mod1coefs[[i]] <- lapply(t(split(mod1coefs[[i]], g1)), matrix, nrow = 6)
        for(j in 1:p){
          colnames(mod0coefs[[i]][[j]]) <- varsMod0[[i]][[j]][1:ncol(mod0coefs[[i]][[j]])]
          colnames(mod1coefs[[i]][[j]]) <- varsMod1se[[i]][[j]][1:ncol(mod1coefs[[i]][[j]])]
          rownames(mod0coefs[[i]][[j]]) <- rownames(mod1coefs[[i]][[j]]) <- c("b", "se", "t", "P", "lower", "upper")
          mod0coefs[[i]][[j]] <- mod0coefs[[i]][[j]][, match(allNames, colnames(mod0coefs[[i]][[j]]))]
          mod1coefs[[i]][[j]] <- mod1coefs[[i]][[j]][, match(allNames, colnames(mod1coefs[[i]][[j]]))]
          colnames(mod0coefs[[i]][[j]]) <- colnames(mod1coefs[[i]][[j]]) <- allNames
          if(bonf == TRUE){
            mod0coefs[[i]][[j]]["P",] <- pmin(mod0coefs[[i]][[j]]["P",] * mod0sizes[j,i], 1)
            mod1coefs[[i]][[j]]["P",] <- pmin(mod1coefs[[i]][[j]]["P",] * mod1sizes[j,i], 1)
          }
        }
      }
      finalCoefs0 <- list(); finalCoefs1 <- list()
      n <- ifelse(sampMethod == "bootstrap", 0, n)
      for(j in 1:p){
        finalCoefs0[[j]] <- vector("list", length = 6)
        finalCoefs0[[j]][[1]] <- t(sapply(mod0coefs, function(z) z[[j]]["b",]))
        finalCoefs0[[j]][[2]] <- t(sapply(mod0coefs, function(z) z[[j]]["se",]))
        if(any(is.na(finalCoefs0[[j]][[2]]))){finalCoefs0[[j]][[2]][is.na(finalCoefs0[[j]][[2]])] <- Inf}
        finalCoefs0[[j]][[3]] <- t(sapply(mod0coefs, function(z) z[[j]]["P",]))
        if(any(is.na(finalCoefs0[[j]][[3]]))){finalCoefs0[[j]][[3]][is.na(finalCoefs0[[j]][[3]])] <- 1}
        finalCoefs0[[j]][[4]] <- t(sapply(mod0coefs, function(z) z[[j]]["lower",]))
        if(any(is.na(finalCoefs0[[j]][[4]]))){finalCoefs0[[j]][[4]][is.na(finalCoefs0[[j]][[4]])] <- -Inf}
        finalCoefs0[[j]][[5]] <- t(sapply(mod0coefs, function(z) z[[j]]["upper",]))
        if(any(is.na(finalCoefs0[[j]][[5]]))){finalCoefs0[[j]][[5]][is.na(finalCoefs0[[j]][[5]])] <- Inf}
        finalCoefs0[[j]][[6]] <- unname((N - n - mod0sizes - 1)[j,])
        finalCoefs1[[j]] <- vector("list", length = 6)
        finalCoefs1[[j]][[1]] <- t(sapply(mod1coefs, function(z) z[[j]]["b",]))
        finalCoefs1[[j]][[2]] <- t(sapply(mod1coefs, function(z) z[[j]]["se",]))
        if(any(is.na(finalCoefs1[[j]][[2]]))){finalCoefs1[[j]][[2]][is.na(finalCoefs1[[j]][[2]])] <- Inf}
        finalCoefs1[[j]][[3]] <- t(sapply(mod1coefs, function(z) z[[j]]["P",]))
        if(any(is.na(finalCoefs1[[j]][[3]]))){finalCoefs1[[j]][[3]][is.na(finalCoefs1[[j]][[3]])] <- 1}
        finalCoefs1[[j]][[4]] <- t(sapply(mod1coefs, function(z) z[[j]]["lower",]))
        if(any(is.na(finalCoefs1[[j]][[4]]))){finalCoefs1[[j]][[4]][is.na(finalCoefs1[[j]][[4]])] <- -Inf}
        finalCoefs1[[j]][[5]] <- t(sapply(mod1coefs, function(z) z[[j]]["upper",]))
        if(any(is.na(finalCoefs1[[j]][[5]]))){finalCoefs1[[j]][[5]][is.na(finalCoefs1[[j]][[5]])] <- Inf}
        finalCoefs1[[j]][[6]] <- unname((N - n - mod1sizes - 1)[j,])
        names(finalCoefs0[[j]]) <- names(finalCoefs1[[j]]) <- c("b", "se", "P", "lower", "upper", "df.res")
      }
      names(finalCoefs0) <- names(finalCoefs1) <- colnames(dat$Y)
      if(niter > 1){
        gamma <- seq(ceiling(.05 * niter)/niter, 1 - 1/niter, by = 1/niter)
        penalty <- ifelse(length(gamma) > 1, (1 - log(min(gamma))), 1)
        k <- ncol(finalCoefs0$X1.y$P)
        pvals0 <- pvals1 <- gamma0 <- gamma1 <- matrix(0, ncol = k, nrow = p)
        for(i in 1:p){
          for(j in 1:k){
            quantGam0 <- quantile(finalCoefs0[[i]][["P"]][,j], gamma)/gamma
            pvals0[i,j] <- pmin((min(quantGam0) * penalty), 1)
            gamma0[i,j] <- gamma[which.min(quantGam0)]
            quantGam1 <- quantile(finalCoefs1[[i]][["P"]][,j], gamma)/gamma
            pvals1[i,j] <- pmin((min(quantGam1) * penalty), 1)
            gamma1[i,j] <- gamma[which.min(quantGam1)]
          }
        }
        colnames(pvals0) <- colnames(pvals1) <- colnames(gamma0) <- colnames(gamma1) <- colnames(finalCoefs0$X1.y$P)
        rownames(pvals0) <- rownames(pvals1) <- rownames(gamma0) <- rownames(gamma1) <- colnames(dat$Y)
        adjCIs0 <- list(); adjCIs1 <- list()
        sInd <- rep(1:k, each = niter)
        for(i in 1:p){
          adjCIs0[[i]] <- t(mapply(adjustedCI, 
                                   lci = split(finalCoefs0[[i]][["lower"]], sInd),
                                   rci = split(finalCoefs0[[i]][["upper"]], sInd), 
                                   centers = split(finalCoefs0[[i]][["b"]], sInd),
                                   ses = split(finalCoefs0[[i]][["se"]], sInd),
                                   df.res = list(df.res = finalCoefs0[[i]][["df.res"]]),
                                   gamma.min = min(gamma), ci.level = (1 - alpha), var = 1:k))
          colnames(adjCIs0[[i]]) <- c("lower", "upper")
          rownames(adjCIs0[[i]]) <- colnames(pvals0)
          adjCIs0[[i]] <- data.frame(adjCIs0[[i]])
          adjCIs0[[i]][adjCIs0[[i]]$lower == -Inf, ] <- NA
          adjCIs0[[i]]$b <- rowMeans(adjCIs0[[i]])
          adjCIs0[[i]]$predictor <- factor(rownames(adjCIs0[[i]]))
          rownames(adjCIs0[[i]]) <- 1:nrow(adjCIs0[[i]])
          adjCIs0[[i]] <- data.frame(predictor = adjCIs0[[i]]$predictor, lower = adjCIs0[[i]]$lower, 
                                     b = adjCIs0[[i]]$b, upper = adjCIs0[[i]]$upper, 
                                     Pvalue = pvals0[i,], gammaMin = gamma0[i,])
          if(any(rowSums(is.na(adjCIs0[[i]])) != 0)){
            adjCIs0[[i]][rowSums(is.na(adjCIs0[[i]])) != 0, c("Pvalue", "gammaMin")] <- NA}
          adjCIs1[[i]] <- t(mapply(adjustedCI, 
                                   lci = split(finalCoefs1[[i]][["lower"]], sInd),
                                   rci = split(finalCoefs1[[i]][["upper"]], sInd), 
                                   centers = split(finalCoefs1[[i]][["b"]], sInd),
                                   ses = split(finalCoefs1[[i]][["se"]], sInd),
                                   df.res = list(df.res = finalCoefs1[[i]][["df.res"]]),
                                   gamma.min = min(gamma), ci.level = (1 - alpha), var = 1:k))
          colnames(adjCIs1[[i]]) <- c("lower", "upper")
          rownames(adjCIs1[[i]]) <- colnames(pvals1)
          adjCIs1[[i]] <- data.frame(adjCIs1[[i]])
          adjCIs1[[i]][adjCIs1[[i]]$lower == -Inf, ] <- NA
          adjCIs1[[i]]$b <- rowMeans(adjCIs1[[i]])
          adjCIs1[[i]]$predictor <- factor(rownames(adjCIs1[[i]]))
          rownames(adjCIs1[[i]]) <- 1:nrow(adjCIs1[[i]])
          adjCIs1[[i]] <- data.frame(predictor = adjCIs1[[i]]$predictor, lower = adjCIs1[[i]]$lower, 
                                     b = adjCIs1[[i]]$b, upper = adjCIs1[[i]]$upper, 
                                     Pvalue = pvals1[i,], gammaMin = gamma1[i,])
          if(any(rowSums(is.na(adjCIs1[[i]])) != 0)){
            adjCIs1[[i]][rowSums(is.na(adjCIs1[[i]])) != 0, c("Pvalue", "gammaMin")] <- NA}
        }
        names(adjCIs0) <- names(adjCIs1) <- colnames(dat$Y)
      } else {
        adjCIs0 <- adjCIs1 <- NA
      }
      for(i in 1:p){
        finalCoefs0[[i]]$se[is.infinite(finalCoefs0[[i]]$se)] <- NA
        finalCoefs0[[i]]$lower[is.infinite(finalCoefs0[[i]]$lower)] <- NA
        finalCoefs0[[i]]$upper[is.infinite(finalCoefs0[[i]]$upper)] <- NA
        finalCoefs1[[i]]$se[is.infinite(finalCoefs1[[i]]$se)] <- NA
        finalCoefs1[[i]]$lower[is.infinite(finalCoefs1[[i]]$lower)] <- NA
        finalCoefs1[[i]]$upper[is.infinite(finalCoefs1[[i]]$upper)] <- NA
      }
    } else {
      finalCoefs0 <- finalCoefs1 <- adjCIs0 <- adjCIs1 <- NA
    }
    varsMod0 <- lapply(varsMod0, function(modIter) do.call(cbind.data.frame, modIter))
    varsMod1se <- lapply(varsMod1se, function(modIter) do.call(cbind.data.frame, modIter))
    names(varsMod0) <- names(varsMod1se) <- paste0("iter", 1:niter)
    mod0fits <- lapply(fits, function(modIter) modIter$mod0)
    mod1fits <- lapply(fits, function(modIter) modIter$mod1se)
    mod0 <- list(CIs = adjCIs0, freqs = mod0Freqs, coefs = finalCoefs0, 
                 fits = mod0fits, iterVars = varsMod0)
    mod1se <- list(CIs = adjCIs1, freqs = mod1Freqs, coefs = finalCoefs1,
                   fits = mod1fits, iterVars = varsMod1se)
    out <- list(mod0 = mod0, mod1se = mod1se, samps = samps, vars = vars, sampInd = sampInd)
  } else {
    out <- list(fits = fits, samps = samps, vars = vars, sampInd = sampInd)
  }
  attributes(out)$call <- match.call()
  out
}