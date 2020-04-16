### ======================================================================== ###
### ======================================================================== ###
##### params: set global parameters for testing functions
params <- function(p = "varselect"){
  p <- match.arg(tolower(p), c("varselect", "fitnetwork", "resample", "eqs", 
                               "lagmat", "plotnet", "modnet", "condplot", 
                               "surll", "mlgvar", "gvar", "trevmlvarsim"))
  if(p == "trevmlvarsim"){
    out <- list(nTime = 50, nVar = 3, nPerson = 20, propPositive = 0.5, 
                kappaRange = c(0.25, 0.5), betaRange = c(0.25, 0.5), 
                betweenRange = c(0.25, 0.5), rewireWithin = 0, 
                betweenVar = 1, withinVar = 0.25, temporalOffset = 2)
  }
  if(p == "mlgvar"){
    out <- list(vars = NULL, idvar = NULL, subjectNetworks = TRUE, useTrevor = TRUE, 
                beepvar = NULL, dayvar = NULL, scale = TRUE, centerWithin = TRUE, 
                gamma = 0.5, verbose = TRUE, lambda_min_kappa_fixed = 0.001,
                lambda_min_beta_fixed = 0.001, lambda_min_kappa = 0.05, 
                lambda_min_beta = 0.05, lambda_min_glasso = 0.01)
  }
  if(p == "gvar"){
    out <- list(m = NULL, idvar = "ID", subjectNets = FALSE, exogenous = TRUE,
                center = TRUE, scale = TRUE, fixedType = "g", betweenType = "g",
                centerWithin = TRUE, scaleWithin = FALSE, selectFUN = NULL,
                verbose = TRUE)
  }
  if(p == "varselect"){
    out <- list(m = NULL, criterion = "AIC", method = "glinternet",
                lags = NULL, exogenous = TRUE, type = "g", center = TRUE,
                scale = FALSE, gamma = .5, nfolds = 10, varSeed = NULL, 
                useSE = TRUE, nlam = NULL, covs = NULL, verbose = TRUE,
                outMsgs = FALSE, diag = FALSE)
  }
  if(p == "eqs"){
    out <- list(varMods = NULL, mod = "min", type = "g", m = NULL, center = TRUE,
                scale = FALSE, exogenous = TRUE, covs = NULL)
  }
  if(p == "lagmat"){
    out <- list(type = "g", m = NULL, covariates = NULL, center = TRUE,
                scale = FALSE, exogenous = TRUE, lags = 1, checkType = FALSE)
  }
  if(p == "fitnetwork"){
    out <- list(moderators = 1, type = "gaussian", seed = NULL,
                lags = 1, lambda = NULL, folds = 10, gamma = 0.25,
                which.lam = 'lambda.min', rule = "OR", threshold = FALSE,
                scale = FALSE, std = TRUE, group = NULL, grPenalty = "gel",
                adaGam = NULL, adaMethod = "ridge", measure = "deviance",
                alpha = 1, saveData = TRUE, center = TRUE, covariates = NULL, 
                verbose = FALSE, exogenous = TRUE, binary = NULL, 
                residMat = "sigma", medges = 1, mval = NULL)
  }
  if(p == "resample"){
    out <- list(niter = 3, m = 1, criterion = "AIC", sampMethod = "split", 
                method = "glinternet", rule = "OR", gamma = .25, nfolds = 10, 
                nlam = 50, which.lam = "min", threshold = FALSE, bonf = FALSE, 
                alpha = .05, exogenous = TRUE, split = .5, center = TRUE, 
                varSeed = NULL, seed = NULL, verbose = TRUE, binary = NULL, 
                saveMods = FALSE, type = "g", lags = 1, scale = FALSE)
  }
  if(p == "plotnet"){
    out <- list(which.net = "temporal", predict = NULL, names = TRUE, 
                scale = FALSE, lag = NULL, nodewise = FALSE, con = "adjR2", 
                cat = "nCC", layout = "spring", elabs = FALSE, elsize = 1, 
                mnet = FALSE, plot = TRUE)
  }
  if(p == "modnet"){
    out <- list(threshold = FALSE, rule = "AND", mval = NULL, 
                pcor = FALSE, useCIs = FALSE, nsims = 5000, mlty = 2)
  }
  if(p == "condplot"){
    out <- list(swap = FALSE, avg = FALSE, compare = NULL, 
                hist = FALSE, xlab = NULL, mods = NULL, nsims = 500, 
                xn = NULL, getCIs = FALSE, discrete = FALSE, 
                ylab = NULL, main = NULL, midline = TRUE)
  }
  if(p == "surll"){
    out <- list(net1 = NULL, all = FALSE, lrt = NULL, nodes = FALSE, 
                d = 4, alpha = .05, s = "res", full = FALSE)
  }
  list2env(out, envir = .GlobalEnv)
}


pars <- function(fun, ..., add = FALSE, eval = FALSE, argList = FALSE){
  args0 <- list(...)
  out <- formals(fun)
  out <- out[setdiff(names(out), c("dat", "data", "..."))]
  if(length(args0) > 0){
    args <- args0[intersect(names(args0), names(out))]
    out[names(args)] <- args
    if(isTRUE(add)){out <- append(out, args0[setdiff(names(args0), names(out))])}
  }
  if(any(sapply(out, class) == 'call') & isTRUE(eval)){
    oc <- which(sapply(out, class) == 'call')
    for(i in seq_along(oc)){out[[oc[i]]] <- eval(out[[oc[i]]])}
  }
  list2env(out, .GlobalEnv)
  if(argList){list2env(list(args = out), .GlobalEnv)}
}

### ======================================================================== ###
### ======================================================================== ###
##### makeFits
makeFits <- function(data, m = NULL, exo = TRUE, covs = NULL, crit = "CV", ...){
  fit0 <- fitNetwork(data, m, lags = 1, covariates = covs, exogenous = exo, ...)
  if(!is.null(m)){
    vars <- varSelect(data, m, crit, lags = 1, covs = covs, exogenous = exo, ...)
    fit1 <- fitNetwork(data, m, vars, lags = 1, covariates = covs, exogenous = exo, ...)
    if(crit == "CV"){
      fit2 <- fitNetwork(data, m, vars, lags = 1, covariates = covs, 
                         exogenous = exo, which.lam = "1se", ...)
    } else {fit2 <- fit1}
  } else {fit2 <- fit1 <- fit0}
  out <- list(fits = list(fit0, fit1, fit2))
  if(!is.null(m)){out$vars <- vars}
  list2env(out, .GlobalEnv)
}

##### makeFits2
makeFits2 <- function(x = c("dep1", "sur2")){
  x <- match.arg(x)
  settings(dat = x)
  if(x == "dep1"){
    fit0 <- fitNetwork(data, pcor = TRUE)
    fit1 <- fitNetwork(data, 5, exogenous = FALSE)
    vars <- varSelect(data, 5, "EBIC", gamma = .5, exogenous = FALSE)
    fit2 <- fitNetwork(data, 5, vars, exogenous = FALSE)
    fits0 <- list(fit0, fit1, fit2)
    names(fits0) <- paste0("fit_dep", 0:2)
    out <- list(dat0 = data, fits0 = fits0)
  } else {
    fit0 <- fitNetwork(data, lags = 1)
    fit1 <- fitNetwork(data, 1, lags = 1, exogenous = FALSE)
    vars <- varSelect(data, 1, "EBIC", lags = 1, exogenous = FALSE, gamma = .5)
    fit2 <- fitNetwork(data, 1, vars, lags = 1, exogenous = FALSE)
    fits1 <- list(fit0, fit1, fit2)
    names(fits1) <- paste0("fit_sur", 0:2)
    out <- list(dat1 = data, fits1 = fits1)
  }
  list2env(out, .GlobalEnv)
  msg <- paste0("dat", ifelse(x == "dep1", 0, 1), " and fits", 
                ifelse(x == "dep1", 0, 1), " in .GlobalEnv")
  return(message(msg))
}


### ======================================================================== ###
### ======================================================================== ###
##### setupFit: run lines of fitNetwork
setupFit <- function(network = c("temporal", "cross"), nodewise = FALSE){
  network <- match.arg(network)
  if(is.null(lags)){network <- "cross"}
  ff <- readLines("~/C: Desktop/COMPS/METHODS/CODE/modnets/functions.R")
  fstart <- which(ff == "  t1 <- Sys.time() # START")
  fsetup <- which(ff == "  ##### SETUP #####")
  fargs <- which(ff == "    args <- list(...) # FITNETWORK ARGS")
  fcross1 <- which(ff == "    ##### CROSS-SECTIONAL #####")
  slines("functions.R", fstart:(fsetup - 1))
  slines("functions.R", c((fsetup + 2):(fargs - 1), (fargs + 2):(fcross1 - 1)))
  if(network == "temporal"){
    ftemp1 <- which(ff == "      ##### TEMPORAL #####")
    ftemp2 <- which(ff == "      net <- SURnet(fit = fit, dat = dat, s = residMat, m = moderators, ")
    slines("functions.R", c((ftemp1 + 1):(ftemp2 - 1)))
    assign("s", "sigma", envir = .GlobalEnv)
    assign("m", moderators, envir = .GlobalEnv)
  } else {
    fcross2 <- which(ff == "      mods0 <- nodewise(data = data, mods = moderators, varMods = type, ")
    slines("functions.R", c((fcross1 + 2):(fcross2 - 1)))
    if(nodewise){
      slines("functions.R", c(fcross2:(fcross2 + 2)))
      assign("models", mods0, envir = .GlobalEnv)
      assign("varMods", type, envir = .GlobalEnv)
      assign("mods", moderators, envir = .GlobalEnv)
      assign("getEqs", FALSE, envir = .GlobalEnv)
      assign("lambda", "min", envir = .GlobalEnv)
    }
  }
}

##### setupRes: test iterator in resample function
setupRes <- function(ii, ss = c("split", "bootstrap", "stability")){
  ss <- match.arg(ss)
  assign("sampMethod", ss, envir = .GlobalEnv)
  ff <- readLines("~/C: Desktop/COMPS/METHODS/CODE/modnets/functions.R")
  fstart <- which(ff == "  t1 <- Sys.time() # RESAMPLE START")
  fargs <- which(ff == "  args <- list(...) ##### RESAMPLE ARGS")
  fsetup <- which(ff == "  ##### RESAMPLE #####")
  fend <- which(ff == "  ##### RESAMPLE END #####")
  slines("functions.R", c(fstart:(fargs - 1), (fargs + 2):(fsetup - 1)))
  if(length(ii) == 1 & all(ii != 0)){ii <- 1:ii}
  if(all(ii == 0)){return(ii)}
  if(ss != "bootstrap"){
    slines("functions.R", (fsetup + 2):(fsetup + 4))
    fsplit1 <- which(ff == "      set.seed(seeds[i]) # RESAMPLE SPLIT")
    fstab <- which(ff == "      } # RESAMPLE STABILITY")
    fsplit2 <- which(ff == "        } # RESAMPLE SPLIT FIT")
    for(jj in ii){
      assign("i", jj, envir = .GlobalEnv)
      slines("functions.R", fsplit1:fstab)
      if(ss == "stability"){
        slines("functions.R", (fstab + 2):fsplit2)
      } else {
        slines("functions.R", (fsplit2 + 2):(fsplit2 + 6))
      }
    }
  } else {
    fboot <- which(ff == "      set.seed(seeds[i]) # RESAMPLE BOOT")
    for(jj in ii){
      assign("i", jj, envir = .GlobalEnv)
      slines("functions.R", fboot:(fend - 3))
    }
  }
  fvars <- which(ff == "  ##### RESAMPLE VARMODS #####")
  slines("functions.R", (fend + 1):(fvars - 1))
}

##### testVars: test lines of varSelect
testVars <- function(ii = 0){
  ff <- readLines("~/C: Desktop/COMPS/METHODS/CODE/modnets/functions.R")
  if(ii == 0){
    assign("method", "glmnet", envir = .GlobalEnv)
    assign("m", NULL, envir = .GlobalEnv)
    fsetup <- which(ff == "  dmnames <- NULL # VARSELECT START")
    fstart <- which(ff == "  ##### VARSELECT #####")
    slines("functions.R", fsetup:(fstart - 1))
  } else {
    flasso <- which(ff == "  ##### LASSO SELECTION #####")
    assign("i", ii, envir = .GlobalEnv)
    assign("yvar", ii, envir = .GlobalEnv)
    assign("nlam", 100, envir = .GlobalEnv)
    assign("alpha", 1, envir = .GlobalEnv)
    slines("functions.R", (flasso + 4):(flasso + 8))
  }
}


### ======================================================================== ###
### ======================================================================== ###
##### blockBoot: split timepoints into blocks, then bootstrap within blocks
blockBoot <- function(data, boots, blocks = 10, lags = 1){
  n <- 1:(nrow(data) - lags)
  Qt <- quantile(n, probs = seq(0, 1, length = blocks + 1))
  ind_blocks <- cut(x = n, breaks = Qt, labels = FALSE)
  ind_blocks[1] <- 1
  ind <- list()
  for(i in 1:boots){
    ind_blocks2 <- list()
    for(j in 1:blocks){
      block2 <- which(ind_blocks == j)
      ind_blocks2[[j]] <- sample(x = block2, size = length(block2), replace = TRUE)
    }
    ind[[i]] <- unlist(ind_blocks2)
  }
  return(ind)
}

##### checkSel: check whether select == select_ci
checkSel <- function(x){sapply(x[["adjCIs"]], function(z) all(z$select == z$select_ci))}

