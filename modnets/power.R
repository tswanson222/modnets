##### mnetPowerSim: power simulator for cross-sectional and idiographic networks
mnetPowerSim <- function(niter = 10, N = 100, p = 5, m = FALSE, m1 = 0, m2 = .1, sparsity = .5,
                         lags = NULL, trueNet = NULL, threshold = TRUE, rule = 'OR', avg = TRUE, 
                         maxiter = 100, saveFits = TRUE, saveData = FALSE, intercepts = NULL, 
                         mbinary = FALSE, select = NULL, vargs = list(), type = 'g', 
                         gibbs = TRUE, ordinal = FALSE, mord = FALSE, nLevels = 5, 
                         minOrd = 3, div = 1000, modType = 'none', m1_range = NULL, 
                         m2_range = c(.1, .3), time = TRUE, skewErr = FALSE, 
                         nCores = 1, cluster = 'mclapply', fixedPar = NULL, V2 = 1, ...){
  t1 <- Sys.time()
  runagain <- FALSE
  if(is.null(threshold)){threshold <- is.null(select)}
  if(length(p) > 1 | length(m2) > 1 | length(sparsity) > 1){
    runagain <- TRUE
    rerun <- expand.grid(p, m2, sparsity, stringsAsFactors = FALSE)
    colnames(rerun) <- c('p', 'm2', 'sparsity')
    p <- rerun[1, 'p']
    m2 <- rerun[1, 'm2']
    sparsity <- rerun[1, 'sparsity']
    rerun <- rerun[-1, ]
    runs <- nrow(rerun)
  }
  if(!is.null(select)){
    select <- ifelse(isTRUE(select), 'varSelect', match.arg(select, c('varSelect', 'resample')))
    vargs$m <- switch(2 - identical(m, FALSE), NULL, p + 1)
    vargs$rule <- rule
    vargs$lags <- lags
    vargs$exogenous <- TRUE
    vargs$verbose <- FALSE
    vargs <- vargs[intersect(names(vargs), formalArgs(select))]
    if(select == 'resample'){select <- function(...){modSelect(resample(...))}}
  }
  parms <- list(N = N)[length(N) >= 2]
  if(!is.null(lags) & length(parms) == 1){names(parms) <- 'Nt'}
  args <- tryCatch({list(...)}, error = function(e){list()})
  trueNet2 <- vars <- NULL
  if(is.null(trueNet)){
    if(is.null(lags)){
      mat <- b2 <- diag(0, p)
      pn <- (p * (p - 1))/2
      if(V2 == 1 & m2 < 1){m2 <- ceiling(m2 * pn)}
      while(all(b2 == 0) & modType != 'zero'){
        if(m2 >= 0 & m2 < 1){
          b2 <- sample(0:1, pn, TRUE, prob = c(1 - m2, m2))
        } else if(m2 >= 1){
          b2 <- numeric(pn)
          b2[sample(1:pn, ifelse(m2 > pn, pn, round(m2)))] <- 1
        }
        b2 <- b2 * sample(c(-1, 1), pn, TRUE, prob = c(.5, .5))
        mat[upper.tri(mat)] <- b2 * runif(pn, min(m2_range), max(m2_range))
        b2 <- as.matrix(Matrix::forceSymmetric(mat))
        if(!is.null(fixedPar)){
          b2[b2 != 0] <- switch(length(fixedPar), fixedPar, fixedPar[2]) * sign(b2[b2 != 0])
        }
      }
      args$b2 <- trueNet2 <- b2
      args1 <- list(Nvar = p, sparsity = sparsity, finalRange = m1_range)
      if(!identical(m, FALSE)){
        if(is.numeric(modType)){modType <- switch(modType + 1, 'none', 'full', 'partial')}
        modType <- match.arg(tolower(modType), c('none', 'full', 'partial'))
      } else {
        modType <- 'none'
      }
      trueNet1 <- do.call(simPcor, args1)
      if(!is.null(fixedPar)){trueNet1[trueNet1 != 0] <- fixedPar[1] * sign(trueNet1[trueNet1 != 0])}
      if(FALSE){
        en1 <- eigen(diag(p) - trueNet1)$values
        while(any(round(en1, 14) == 1)){
          trueNet1 <- do.call(simPcor, args1)
          if(!is.null(fixedPar)){trueNet1[trueNet1 != 0] <- fixedPar[1] * sign(trueNet1[trueNet1 != 0])}
          en1 <- eigen(diag(p) - trueNet1)$values
        }
      }
      if(modType != 'none'){
        while(ifelse(modType == 'full', !all(trueNet1[trueNet2 != 0] == 0), any(trueNet1[trueNet2 != 0] == 0))){
          trueNet1 <- do.call(simPcor, args1)
          if(!is.null(fixedPar)){trueNet1[trueNet1 != 0] <- fixedPar[1] * sign(trueNet1[trueNet1 != 0])}
          if(FALSE){
            en1 <- eigen(diag(p) - trueNet1)$values
            while(any(round(en1, 14) == 1)){
              trueNet1 <- do.call(simPcor, args1)
              if(!is.null(fixedPar)){trueNet1[trueNet1 != 0] <- fixedPar[1] * sign(trueNet1[trueNet1 != 0])}
              en1 <- eigen(diag(p) - trueNet1)$values
            }
          }
        }
      }
      args$b1 <- trueNet1
    } else {
      avg <- FALSE
      m0 <- ifelse(identical(m, FALSE), FALSE, ifelse(mbinary, 'binary', ifelse(mord, 'ordinal', m)))
      out1 <- mlGVARsim(nTime = N[1], nPerson = 1, nNode = p, m = m0, m2 = m2, 
                        onlyNets = TRUE, ordinal = ordinal, nLevels = nLevels, 
                        GGMsparsity = sparsity, m1_range = m1_range, 
                        m2_range = m2_range, minOrd = minOrd, 
                        skewErr = skewErr, m1 = m1, modType = modType)
      if(any(eigen(out1$residuals)$values == 1) & FALSE){
        while(any(eigen(out1$residuals)$values == 1)){
          out1 <- mlGVARsim(nTime = N[1], nPerson = 1, nNode = p, m = m0, m2 = m2, 
                            onlyNets = TRUE, ordinal = ordinal, nLevels = nLevels, 
                            GGMsparsity = sparsity, m1_range = m1_range, 
                            m2_range = m2_range, minOrd = minOrd, 
                            skewErr = skewErr, m1 = m1, modType = modType)
        }
      }
      args <- append(args, out1)
      trueNet1 <- list(beta = out1$parms[[1]], 
                       kappa = round(corpcor::pseudoinverse(out1$residuals), 14),
                       PCC = round(pcor2(out1$residuals), 14))
      trueNet2 <- out1$mb2
    }
  } else if(is(trueNet, 'list')){
    args$b1 <- trueNet1 <- trueNet[[1]]
    args$b2 <- trueNet2 <- trueNet[[2]]
  } else {
    args$b1 <- trueNet1 <- trueNet
    args$b2 <- diag(0, p)
  }
  FUN <- switch(2 - is.null(lags), 'simNet', 'trevSimulateVAR')
  args2 <- args[intersect(names(args), setdiff(c(formalArgs(FUN), 'divup', 'nostop'), '...'))]
  if(!is.null(intercepts)){args2$intercepts <- intercepts}
  if(is.null(lags)){
    args2 <- append(args2, list(
      N = N, p = p, m = m, m2 = m2, m1 = m1, onlyDat = TRUE, pbar = FALSE, 
      gibbs = gibbs, ordinal = ordinal, mord = mord, time = FALSE,
      mbinary = mbinary, nLevels = nLevels, skewErr = skewErr))
  }
  if(nCores > 1 | isTRUE(nCores)){
    pbapply::pboptions(type = 'timer', char = '-')
    if(isTRUE(nCores)){nCores <- parallel::detectCores()}
    if(grepl('Windows', sessionInfo()$running)){cluster <- 'SOCK'}
    if(tolower(cluster) != 'mclapply'){
      cluster <- match.arg(toupper(cluster), c('SOCK', 'FORK'))
      cl <- parallel::makeCluster(nCores, type = cluster)
      if(cluster == 'SOCK'){
        if(!is.null(lags)){
          suppressMessages(invisible(sapply(c('mvtnorm', 'sn'), require, character.only = TRUE)))
          objects <- c('rmvnorm', 'rmsn', 'lagMat', 'SURfit', 'SURnet', 
                       'SURll', 'surEqs', 'getCoefs', 'systemfit')
        } else {
          objects <- c('simPcor', 'nodewise', 'modNet', 'modLL')
        }
        if(is.character(select)){
          objects <- c(objects, 'varSelect', 'Matrix', ifelse(!identical(m, FALSE), 'glinternet', ifelse(
            'method' %in% names(vargs), switch(vargs$method, subset = 'regsubsets', vargs$method), 'glmnet')))
          objects <- c(objects, ifelse('glinternet' %in% objects, 'fitHierLASSO', ifelse(
            'glmnet' %in% objects, 'lassoSelect', '')))
          if('criterion' %in% names(vargs)){
            if(toupper(vargs$criterion) == 'CV'){
              objects <- gsub(ifelse('glmnet' %in% objects, 'glmnet', 'glinternet'), ifelse(
                'glmnet' %in% objects, 'cv.glmnet', 'glinternet.cv'), objects)
            }
          }
        }
        parallel::clusterExport(cl, c(FUN, 'fitNetwork', objects, 'simNet', 'simNet2', 'getFitCIs'), envir = environment())
      }
    } else {
      cl <- nCores
    }
  }
  if(FUN == 'simNet'){FUN <- 'simNet2'}
  Data <- Fits <- list()
  if(length(parms) == 0){
    if(nCores > 1){
      OUT <- pbapply::pblapply(seq_len(niter), function(i){
        t0 <- 0
        Data <- TRUE
        while(isTRUE(Data)){
          if(!is.null(lags) & !identical(m, FALSE)){
            args2$m <- rnorm(args2$Nt + args2$burnin, 0, 1)
            if(mord){
              ord <- c()
              while(length(unique(ord)) < minOrd){
                ord <- as.numeric(cut(args2$m, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
              }
              args2$m <- ord
            }
          }
          args3 <- switch(2 - (FUN == 'simNet2'), list(nets = args2), args2)
          Data <- tryCatch({do.call(match.fun(FUN), args3)}, error = function(e){TRUE})
          t0 <- t0 + 1
          if(t0 >= maxiter | any(Data > div)){
            if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
            stop('Failed to simulate dataset')
          }
        }
        if(!is.null(lags)){
          if(ordinal){
            Data <- data.frame(apply(Data, 2, function(z){
              ord <- c()
              while(length(unique(ord)) < minOrd){
                ord <- as.numeric(cut(z, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
              }
              return(ord)
            }))
          }
          if(!identical(m, FALSE)){Data$M <- args2$m[-(1:args2$burnin)]}
        }
        if(!is.null(select)){
          vargs$data <- Data
          type <- do.call(match.fun(select), vargs)
        }
        Fits <- fitNetwork(Data, moderators = switch(
          2 - identical(m, FALSE), NULL, p + 1), type = type, 
          rule = rule, saveMods = FALSE, threshold = threshold, 
          lags = lags, fitCoefs = TRUE)
        return(list(Data = Data, Fits = Fits))
      }, cl = cl)
      Data <- lapply(OUT, '[[', 'Data')
      Fits <- lapply(OUT, '[[', 'Fits')
      if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
      rm(cl, OUT)
    } else {
      pb <- txtProgressBar(max = niter, style = 3)
      for(i in seq_len(niter)){
        t0 <- 0
        Data[[i]] <- TRUE
        while(isTRUE(Data[[i]])){
          if(!is.null(lags) & !identical(m, FALSE)){
            args2$m <- rnorm(args2$Nt + args2$burnin, 0, 1)
            if(mord){
              ord <- c()
              while(length(unique(ord)) < minOrd){
                ord <- as.numeric(cut(args2$m, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
              }
              args2$m <- ord
            }
          }
          args3 <- switch(2 - (FUN == 'simNet2'), list(nets = args2), args2)
          Data[[i]] <- tryCatch({do.call(match.fun(FUN), args3)}, error = function(e){TRUE})
          t0 <- t0 + 1
          if(t0 >= maxiter | any(Data[[i]] > div)){
            close(pb)
            stop('Failed to simulate dataset')
          }
        }
        if(!is.null(lags)){
          if(ordinal){
            Data[[i]] <- data.frame(apply(Data[[i]], 2, function(z){
              ord <- c()
              while(length(unique(ord)) < minOrd){
                ord <- as.numeric(cut(z, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
              }
              return(ord)
            }))
          }
          if(!identical(m, FALSE)){Data[[i]]$M <- args2$m[-(1:args2$burnin)]}
        }
        if(!is.null(select)){
          vargs$data <- Data[[i]]
          type <- do.call(match.fun(select), vargs)
        }
        Fits[[i]] <- fitNetwork(Data[[i]], moderators = switch(
          2 - identical(m, FALSE), NULL, p + 1), type = type, 
          rule = rule, saveMods = FALSE, threshold = threshold, 
          lags = lags, fitCoefs = TRUE)
        setTxtProgressBar(pb, i)
      }
      close(pb)
    }
    names(Data) <- names(Fits) <- paste0('iter', 1:niter)
    Data <- list(Data); Fits <- list(Fits)
  } else {
    vars <- expand.grid(parms, stringsAsFactors = FALSE)
    if(nCores > 1){
      vars0 <- rep(N, each = niter)
      OUT <- pbapply::pblapply(seq_len(niter * nrow(vars)), function(i){
        t0 <- 0
        Data <- TRUE
        while(isTRUE(Data)){
          if(!is.null(lags) & 'Nt' %in% colnames(vars)){args2$Nt <- vars0[i]}
          if(!is.null(lags) & !identical(m, FALSE)){
            args2$m <- rnorm(args2$Nt + args2$burnin, 0, 1)
            if(mord){
              ord <- c()
              while(length(unique(ord)) < minOrd){
                ord <- as.numeric(cut(args2$m, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
              }
              args2$m <- ord
            }
          }
          args3 <- replace(args2, colnames(vars), as.list(vars0[i]))
          args3 <- switch(2 - (FUN == 'simNet2'), list(nets = args3), args3)
          Data <- tryCatch({
            do.call(match.fun(FUN), args3)},
            error = function(e){TRUE})
          t0 <- t0 + 1
          if(t0 >= maxiter | any(Data > div)){
            if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
            stop('Failed to simulate dataset')
          }
        }
        if(!is.null(lags)){
          if(ordinal){
            Data <- data.frame(apply(Data, 2, function(z){
              ord <- c()
              while(length(unique(ord)) < minOrd){
                ord <- as.numeric(cut(z, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
              }
              return(ord)
            }))
          }
          if(!identical(m, FALSE)){Data$M <- args2$m[-(1:args2$burnin)]}
        }
        if(!is.null(select)){
          vargs$data <- Data
          type <- do.call(match.fun(select), vargs)
        }
        Fits <- fitNetwork(Data, moderators = switch(
          2 - identical(m, FALSE), NULL, p + 1), type = type, 
          rule = rule, saveMods = FALSE, threshold = threshold, 
          lags = lags, fitCoefs = TRUE)
        return(list(Data = Data, Fits = Fits))
      }, cl = cl)
      Data <- lapply(seq_len(niter), function(i){
        lapply(OUT, '[[', 'Data')[i + ((seq_len(nrow(vars)) - 1) * niter)]
      })
      Fits <- lapply(seq_len(niter), function(i){
        lapply(OUT, '[[', 'Fits')[i + ((seq_len(nrow(vars)) - 1) * niter)]
      })
      if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
      rm(cl, OUT)
    } else {
      pb <- txtProgressBar(max = niter, style = 3)
      for(i in seq_len(niter)){
        Data[[i]] <- lapply(seq_len(nrow(vars)), function(j){
          t0 <- 0
          out0 <- TRUE
          while(isTRUE(out0)){
            if(!is.null(lags) & 'Nt' %in% colnames(vars)){args2$Nt <- vars[j, 'Nt']}
            if(!is.null(lags) & !identical(m, FALSE)){
              args2$m <- rnorm(args2$Nt + args2$burnin, 0, 1)
              if(mord){
                ord <- c()
                while(length(unique(ord)) < minOrd){
                  ord <- as.numeric(cut(args2$m, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
                }
                args2$m <- ord
              }
            }
            args3 <- replace(args2, colnames(vars), as.list(vars[j, ]))
            args3 <- switch(2 - (FUN == 'simNet2'), list(nets = args3), args3)
            out0 <- tryCatch({
              do.call(match.fun(FUN), args3)}, 
              error = function(e){TRUE})
            t0 <- t0 + 1
            if(t0 >= maxiter | any(out0 > div)){
              close(pb)
              stop('Failed to simulate dataset')
            }
          }
          if(!is.null(lags)){
            if(ordinal){
              out0 <- data.frame(apply(out0, 2, function(z){
                ord <- c()
                while(length(unique(ord)) < minOrd){
                  ord <- as.numeric(cut(z, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
                }
                return(ord)
              }))
            }
            if(!identical(m, FALSE)){out0$M <- args2$m[-(1:args2$burnin)]}
          }
          return(out0)
        })
        Fits[[i]] <- lapply(seq_len(nrow(vars)), function(j){
          if(!is.null(select)){
            vargs$data <- Data[[i]][[j]]
            type <- do.call(match.fun(select), vargs)
          }
          fitNetwork(data = Data[[i]][[j]], moderators = switch(
            2 - identical(m, FALSE), NULL, p + 1), type = type, 
            rule = rule, saveMods = FALSE, threshold = threshold, 
            lags = lags, fitCoefs = TRUE)
        })
        setTxtProgressBar(pb, i)
      }
      close(pb)
    }
    Data <- lapply(seq_len(nrow(vars)), function(z){
      setNames(lapply(Data, '[[', z), paste0('iter', 1:niter))
    })
    Fits <- lapply(seq_len(nrow(vars)), function(z){
      setNames(lapply(Fits, '[[', z), paste0('iter', 1:niter))
    })
    args2 <- append(args2, list(vars = vars))
  }
  Nets <- lapply(seq_along(Fits), function(z){
    if(!is(trueNet1, 'list')){trueNet1 <- list(trueNet1)}
    which.net <- c('beta', 'kappa', 'PCC')
    if(length(trueNet1) == 1){which.net <- 'between'}
    net1 <- net0 <- strei1 <- list()
    for(n1 in seq_along(trueNet1)){
      tru <- trueNet1[[n1]]
      z0 <- lapply(Fits[[z]], net, which.net[n1], threshold = threshold, rule = rule)
      z1 <- lapply(z0, function(z){
        cbind(cor = matrixDist(z, tru, 'cor', directed = isTRUE(which.net[n1] == 'beta')), 
              mae = matrixDist(z, tru, 'mae', directed = isTRUE(which.net[n1] == 'beta')),
              performance(z, tru, inds = which.net[n1]))
      })
      net1[[n1]] <- do.call(rbind, z1)
      strei1[[n1]] <- do.call(rbind, lapply(z0, function(z){
        #cent0 <- centAuto(tru)[[1]]
        #cent1 <- centAuto(z)[[1]]
        #btn0 <- cor0(cent0[, 'Betweenness'], cent1[, 'Betweenness'])
        #clo0 <- cor0(cent0[, 'Closeness'], cent1[, 'Closeness'])
        diag(tru) <- diag(z) <- 0
        str0 <- cor0(colSums(abs(tru)), colSums(abs(z)))
        ei0 <- cor0(colSums(tru), colSums(z))
        z2 <- c(Strength = str0, EI = ei0)
        if(which.net[n1] != 'between'){
          str1 <- cor0(rowSums(abs(tru)), rowSums(abs(z)))
          ei1 <- cor0(rowSums(tru), rowSums(z))
          z2 <- c(OutStrength = str0, InStrength = str1, 
                  OutEI = ei0, InEI = ei1)
        }
        #z2 <- c(Betweenness = btn0, Closeness = clo0, z2)
        return(z2)
      }))
      net0[[n1]] <- setNames(do.call(cbind.data.frame, lapply(z0, function(zz){
        switch(2 - (which.net[n1] != 'beta'), zz[lower.tri(zz)], c(zz))
      })), paste0('boot', 1:length(z0)))
    }
    if(length(net0) == 1){
      net0 <- net0[[1]]
    } else {
      names(net0) <- which.net
    }
    net1 <- do.call(rbind, net1)
    strei1 <- do.call(rbind, strei1)
    nstr1 <- cbind(net1, strei1)
    return(list(nstr1, net0))
  })
  nets1 <- lapply(Nets, '[[', 2)
  Nets <- lapply(Nets, '[[', 1)
  allvars1 <- c('N', 'p', 'sparsity')
  for(i in seq_along(Fits)){
    Nets[[i]] <- cbind.data.frame(Nets[[i]], N = NA, p = NA, sparsity = NA)
    if(is.null(vars)){
      Nets[[i]][, allvars1] <- list(N, p, sparsity)
    } else {
      if('Nt' %in% colnames(vars)){colnames(vars) <- gsub('Nt', 'N', colnames(vars))}
      Nets[[i]][, match(colnames(vars), colnames(Nets[[i]]))] <- as.list(as.numeric(vars[i, ]))
      if(length(setdiff(allvars1, colnames(vars))) != 0){
        xx <- list(N = N, p = p, sparsity = sparsity)
        Nets[[i]][, setdiff(allvars1, colnames(vars))] <- xx[setdiff(allvars1, colnames(vars))]
      }
    }
    if(!identical(m, FALSE)){
      if(ifelse(!is.null(vars), 'm2' %in% colnames(vars), FALSE)){
        Nets[[i]]$m2 <- vars[i, 'm2']
      } else {
        Nets[[i]]$m2 <- m2
      }
    }
    Nets[[i]]$index <- i
    Nets[[i]]$iter <- 1:niter
  }
  if(length(Fits) != 1){Nets <- do.call(rbind.data.frame, Nets)}
  if(!is.null(trueNet2) & !identical(m, FALSE)){
    Nets2 <- lapply(seq_along(Fits), function(z){
      z0 <- lapply(Fits[[z]], netInts, threshold = threshold, rule = rule, avg = avg)
      z1 <- lapply(z0, function(zz){
        cbind(cor = matrixDist(zz, trueNet2, 'cor', directed = !avg),
              mae = matrixDist(zz, trueNet2, 'mae', directed = !avg),
              performance(zz, trueNet2, inds = ifelse(avg, 'between', 'beta')))
      })
      z2 <- do.call(rbind, lapply(z0, function(z){
        #cent0 <- centAuto(trueNet2)[[1]]
        #cent1 <- centAuto(z)[[1]]
        #btn0 <- cor0(cent0[, 'Betweenness'], cent1[, 'Betweenness'])
        #clo0 <- cor0(cent0[, 'Closeness'], cent1[, 'Closeness'])
        diag(trueNet2) <- diag(z) <- 0
        str0 <- cor0(colSums(abs(trueNet2)), colSums(abs(z)))
        ei0 <- cor0(colSums(trueNet2), colSums(z))
        z2 <- c(Strength = str0, EI = ei0)
        if(!is.null(lags)){
          str1 <- cor0(rowSums(abs(trueNet2)), rowSums(abs(z)))
          ei1 <- cor0(rowSums(trueNet2), rowSums(z))
          z2 <- c(OutStrength = str0, InStrength = str1, 
                  OutEI = ei0, InEI = ei1)
        }
        #z2 <- c(Betweenness = btn0, Closeness = clo0, z2)
        return(z2)
      }))
      z3 <- cbind(do.call(rbind, z1), z2)
      z4 <- setNames(do.call(cbind.data.frame, lapply(z0, function(zz){
        switch(2 - avg, zz[lower.tri(zz)], c(zz))
      })), paste0('boot', 1:length(z0)))
      return(list(z3, z4))
    })
    nets2 <- lapply(Nets2, '[[', 2)
    Nets2 <- lapply(Nets2, '[[', 1)
    for(i in seq_along(Fits)){
      Nets2[[i]] <- cbind.data.frame(Nets2[[i]], N = NA, p = NA, sparsity = NA)
      if(is.null(vars)){
        Nets2[[i]][, allvars1] <- list(N, p, sparsity)
      } else {
        Nets2[[i]][, match(colnames(vars), colnames(Nets2[[i]]))] <- as.list(as.numeric(vars[i, ]))
        if(length(setdiff(allvars1, colnames(vars))) != 0){
          xx <- list(N = N, p = p, sparsity = sparsity)
          Nets2[[i]][, setdiff(allvars1, colnames(vars))] <- xx[setdiff(allvars1, colnames(vars))]
        }
      }
      if(ifelse(!is.null(vars), 'm2' %in% colnames(vars), FALSE)){
        Nets2[[i]]$m2 <- vars[i, 'm2']
      } else {
        Nets2[[i]]$m2 <- m2
      }
      Nets2[[i]]$index <- i
      Nets2[[i]]$iter <- 1:niter
    }
    if(length(Fits) == 1){
      Nets2 <- Nets2[[1]]
      nets2 <- nets2[[1]]
    } else {
      Nets2 <- do.call(rbind.data.frame, Nets2)
      #names(nets2) <- paste0('N', N)
    }
    rownames(Nets2) <- 1:nrow(Nets2)
  }
  if(length(Data) == 1){
    Data <- Data[[1]]; Fits <- Fits[[1]]; Nets <- Nets[[1]]; nets1 <- nets1[[1]]
  } else if(FALSE){
    names(Data) <- names(Fits) <- names(Nets) <- names(nets1) <- paste0('N', N)
  }
  rownames(Nets) <- 1:nrow(Nets)
  Nets$type <- 'Pairwise'
  Nets$network <- switch(2 - is.null(lags), rep('Between', niter), 
                         rep(c('Beta', 'Kappa', 'PCC'), each = niter))
  if(!identical(m, FALSE)){
    Nets2$type <- 'Interactions'
    Nets2$network <- rep(ifelse(avg, 'Between', 'Beta'), niter)
    Nets <- rbind.data.frame(Nets, Nets2)
  }
  Nets$type <- factor(Nets$type)
  Nets$network <- factor(Nets$network)
  output <- list(Results = Nets, Fits = Fits, Data = Data)
  if(any(is.na(output$Results))){output$Results[is.na(output$Results)] <- 0}
  if(!saveFits){output$Fits <- NULL}
  if(!saveData){output$Data <- NULL}
  output$args <- append(list(call = as.list(match.call())[-1]), args2)
  if(runagain){
    call <- as.list(match.call())[-1]
    call$time <- FALSE
    mpowsim2 <- function(...){
      tryCatch({suppressWarnings(mnetPowerSim(...))}, 
               error = function(e){mpowsim2(...)})
    }
    output2 <- lapply(seq_len(runs), function(z){
      call[c('p', 'm2', 'sparsity')] <- rerun[z, ]
      do.call(mpowsim2, call)
    })
    output2 <- lapply(output2, '[[', 'Results')
    if(length(output2) == 1){
      output2 <- output2[[1]]
    } else {
      output2 <- do.call(rbind, output2)
    }
    output$Results <- do.call(rbind.data.frame, list(output$Results, output2))
  }
  output$nets1 <- nets1
  if(!identical(m, FALSE)){output$nets2 <- nets2}
  class(output) <- c('list', 'mnetPower')
  class(output$Results) <- c('data.frame', 'mnetPower')
  attr(output, 'time') <- Sys.time() - t1
  if(time){print(Sys.time() - t1)}
  return(output)
}

##### plotPower: plot results from power simulation
plotPower <- function(x, by = 'type', yvar = 'default', yadd = NULL, hline = .8,
                      xlab = 'Number of cases', title = NULL, ...){
  args <- list(...)
  if(is(x, 'list')){x <- x$Results}
  if(identical(yvar, 'default')){yvar <- c('sensitivity', 'specificity', 'correlation')}
  if(any(colnames(x) == 'cor')){colnames(x)[colnames(x) == 'cor'] <- 'correlation'}
  if(!is.null(yadd)){yvar <- c(yvar, yadd)}
  cents <- c('Strength', 'EI')
  xvar0 <- c('N', 'p', 'm2', 'index', 'iter', 'type', 'network')
  yvar0 <- setdiff(c(colnames(x), cents, 'fdr'), xvar0)
  yvar <- match.arg(tolower(yvar), unique(tolower(yvar0)), several.ok = TRUE)
  if('fdr' %in% tolower(yvar)){
    x$FDR <- 1 - x$precision
    yvar[tolower(yvar) == 'fdr'] <- 'FDR'
  }
  if(any(tolower(cents) %in% tolower(yvar))){
    for(i in which(tolower(cents) %in% tolower(yvar))){
      yvar <- c(yvar[-which(tolower(yvar) == tolower(cents[i]))], 
                colnames(x)[grepl(cents[i], colnames(x))])
    }
  }
  yvar <- gsub('EI', 'ExpectedInfluence', Hmisc::capitalize(yvar))
  colnames(x) <- gsub('EI', 'ExpectedInfluence', Hmisc::capitalize(colnames(x)))
  if(!is.null(by)){
    by <- match.arg(tolower(by), c(tolower(xvar0), 'pairwise', 'interactions', 'beta', 'kappa', 'pcc'), several.ok = TRUE)
    if('pcc' %in% by){by <- gsub('pcc', 'PCC', by)}
    by <- Hmisc::capitalize(by)
    if(by %in% c('Pairwise', 'Interactions')){
      x <- subset(x, Type == by)
      by <- switch(2 - (length(unique(x$Network)) > 1), 'Network', NULL)
    }
  }
  if(length(unique(x$Type)) != 1 & !identical(by, 'Type')){warning('Interactions and pairwise parameters aggregated')}
  if(length(unique(x$Network)) != 1 & !identical(by, 'Network')){warning('Multiple networks aggregated')}
  if(any(colnames(x) == 'N')){colnames(x)[colnames(x) == 'N'] <- 'nCases'}
  if(length(unique(x$nCases)) == 1){stop('Must have used more than one sample size to plot')}
  FUN <- bootnet:::plot.netSimulator
  args0 <- setdiff(formalArgs(FUN), '...')
  allargs <- formals(FUN)[intersect(names(formals(FUN)), args0)]
  allargs <- replace(allargs, c('x', 'yvar', 'xlab'), list(x = x, yvar = yvar, xlab = xlab))
  if(!is.null(by)){allargs$yfacet <- by}
  args <- args[intersect(names(args), args0[-1])]
  if(length(args) > 0){allargs <- replace(allargs, names(args), args)}
  g <- do.call(FUN, replace(allargs, 'print', FALSE))
  if(!is.null(hline) & !identical(hline, FALSE)){
    g <- g + ggplot2::geom_hline(yintercept = hline, linetype = 2, 
                                 colour = 'red', alpha = .3)
  }
  if(!is.null(title)){g <- g + ggplot2::ggtitle(title)}
  return(g)
}

##### summary.mnetPower: descriptive statistics for power simulation results
summary.mnetPower <- function(x, ind = 'all', order = NULL, decreasing = FALSE){
  if(is(x, 'list')){x <- x$Results}
  inds <- list(N = unique(x$N), p = unique(x$p), sparsity = unique(x$sparsity),
               network = unique(x$network), type = unique(x$type))
  if('m2' %in% colnames(x)){inds$m2 <- unique(x$m2)}
  inds <- expand.grid(inds, stringsAsFactors = FALSE)
  if(length(unique(inds$network)) > 1){
    if(length(unique(inds$type)) > 1){
      inds <- inds[!(inds$type == 'Interactions' & inds$network != 'Beta'), ]
    }
  }
  nn <- colnames(inds)
  x$fdr <- 1 - x$precision
  x$type <- as.character(x$type)
  x$network <- as.character(x$network)
  niter <- max(x$iter)
  ys <- c('cor', 'mae', 'sensitivity', 'specificity', 'precision', 'accuracy', 'fdr')
  means <- sds <- medians <- ses <- list()
  for(i in seq_len(nrow(inds))){
    tdat <- x
    for(j in seq_len(ncol(inds))){
      tdat <- tdat[tdat[, nn[j]] == inds[i, j], ]
    }
    means[[i]] <- apply(tdat[, ys], 2, mean, na.rm = TRUE)
    medians[[i]] <- apply(tdat[, ys], 2, median, na.rm = TRUE)
    sds[[i]] <- apply(tdat[, ys], 2, sd, na.rm = TRUE)
    ses[[i]] <- apply(tdat[, ys], 2, function(z) sd(z, na.rm = TRUE)/sqrt(nrow(tdat)))
  }
  means <- do.call(rbind, means)
  sds <- do.call(rbind, sds)
  medians <- do.call(rbind, medians)
  ses <- do.call(rbind, ses)
  means <- cbind.data.frame(inds, means)
  sds <- cbind.data.frame(inds, sds)
  medians <- cbind.data.frame(inds, medians)
  ses <- cbind.data.frame(inds, ses)
  out <- list(means = means, medians = medians, sds = sds, ses = ses)
  out <- lapply(out, function(z){
    z <- cbind.data.frame(sims = niter, z[order(z$N), ])
    rownames(z) <- 1:nrow(z)
    return(z)
  })
  ind <- switch(2 - identical(ind, 'all'), names(out), 
                match.arg(tolower(ind), names(out), several.ok = TRUE))
  out <- out[ind]
  if(is(out, 'list') & length(out) == 1){out <- out[[1]]}
  if(!is.null(order) & length(ind) == 1){
    if(length(decreasing) != length(order)){decreasing <- rep(decreasing[1], length(order))}
    for(i in seq_along(order)){out <- out[order(out[, order[i]], decreasing = decreasing[i]), ]}
  }
  return(out)
}

##### summary.bootNet: descriptive statistics for bootNet
summary.bootNet <- function(x, centrality = TRUE){
  inds1 <- c('edges', 'strength', 'EI')
  inds2 <- paste0('boot_', inds1)
  inds2[2:3] <- paste0(inds2[2:3], 's')
  if('interactions' %in% names(x)){
    cat('BLANK')
  } else {
    samps <- lapply(inds1, function(z) x[[z]]$sample)
    cors <- lapply(1:3, function(i){
      boots <- x$boots[[inds2[i]]]
      if(!grepl('edge', inds1[i])){boots <- t(boots)}
      sapply(seq_len(ncol(boots)), function(z) cor0(samps[[i]], boots[, z]))
    })
    mae <- lapply(1:3, function(i){
      boots <- x$boots[[inds2[i]]]
      if(!grepl('edge', inds1[i])){boots <- t(boots)}
      sapply(seq_len(ncol(boots)), function(z) mean(abs(samps[[i]] - boots[, z])))
    })
    boots <- x$boots$boot_edges
    sens <- sapply(seq_len(ncol(boots)), function(z){
      sum(boots[, z] != 0 & samps[[1]] != 0)/sum(samps[[1]] != 0)
    })
    spec <- sapply(seq_len(ncol(boots)), function(z){
      sum(boots[, z] == 0 & samps[[1]] == 0)/sum(samps[[1]] == 0)
    })
    prec <- sapply(seq_len(ncol(boots)), function(z){
      sum(boots[, z] != 0 & samps[[1]] != 0)/sum(boots[, z] != 0)
    })
    accu <- sapply(seq_len(ncol(boots)), function(z){
      tp <- sum(boots[, z] != 0 & samps[[1]] != 0)
      tn <- sum(boots[, z] == 0 & samps[[1]] == 0)
      (tp + tn)/length(samps[[1]])
    })
    if(any(is.na(sens))){sens[is.na(sens)] <- 0}
    if(any(is.na(spec))){spec[is.na(spec)] <- 0}
    if(any(is.na(prec))){prec[is.na(prec)] <- 0}
    if(any(is.na(accu))){accu[is.na(accu)] <- 0}
    out <- cbind.data.frame(cor = cors[[1]], mae = mae[[1]], sensitivity = sens, 
                            specificity = spec, precision = prec, accuracy = accu, 
                            N = attr(x, 'n'), p = length(samps[[2]]), 
                            sparsity = sum(samps[[1]] == 0)/length(samps[[1]]),
                            index = 1, iter = 1:ncol(boots), type = 'Pairwise',
                            network = 'Between')
    if(centrality){out <- cbind.data.frame(out, strength = cors[[2]], EI = cors[[3]])}
    class(out) <- c('data.frame', 'mnetPower')
  }
  return(out)
}

##### sampleSize: tells the minimum possible sample size given some fixed p, m, and lags
sampleSize <- function(p, m = 0, lags = 0, print = TRUE){
  if(identical(as.numeric(lags), 0) | is.null(lags)){
    params <- p * (m + 1) + (m * (m - 1))/2
  } else {
    if(!identical(as.numeric(lags), 1)){stop('Lags greater than 1 not supported yet')}
    params <- p * (m + 2) + 1
    #Nmin <- (p + 1) * (m + 2)
  }
  if(isTRUE(print)){
    x0 <- paste(rep('#', 6), collapse = '')
    x1 <- paste0('\n', x0, ' Min. N =')
    x2 <- ifelse(length(strsplit(as.character(params + 1), '')[[1]]) > 2, substring(x0, 1, 5), x0)
    cat(paste0('P = ', p, ' | M = ', m, ' | lags = ', lags), x1, params + 1, x2)
  } else {
    return(params + 1)
  }
}
