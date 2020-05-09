if(!exists("pkgs")){pkgs <- c("glmnet", "mgm", "qgraph", "grpregOverlap", "systemfit", "leaps", "hierNet", "glinternet")}
if(!"none" %in% pkgs){invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))}
if(exists("clear")){if(clear == TRUE){rm(list = ls()); clear <- TRUE}} else {rm(pkgs)}
if(".modnets" %in% search()){detach(".modnets")}
#setwd('~/C: Desktop/COMPS/METHODS/CODE/modnets')
files <- paste0('modnets/', c('ggm', 'centrality', 'sim', 'mlGVAR', 'simGVAR', 'penalized', 'power', 'plots'),'.R')
invisible(sapply(files, source)); rm(files)
if(!'methods' %in% sessionInfo()$basePkgs){invisible(suppressMessages(require(methods)))}
options(stringsAsFactors = FALSE)


### ======================================================================== ###
### ==== MODNETS: Estimating and Analyzing Moderated Temporal Networks ===== ###
### ======================================================================== ###
##### setup: Prepare data for model fitting
setup <- function(data, type, y = NULL, lags = NULL, scale = TRUE, allNum = FALSE){
  if(length(type) == 1 & ncol(data) > 1){type <- rep(type, ncol(data))}
  data <- data.frame(data)
  for(i in which(type == "c")){
    data[,i] <- as.factor(data[,i])
    levels(data[,i]) <- paste(1:length(levels(data[,i])))
  }
  if(scale == TRUE){for(i in which(type == "g")){data[,i] <- scale(data[,i])}}
  data <- data.frame(data)
  names(data) <- paste0("V", 1:ncol(data), ".")
  if(!is.null(lags)){data <- lagIt(data = data, y = y, lags = lags)}
  if(allNum == TRUE & scale == TRUE){
    for(j in 1:ncol(data)){data[,j] <- as.numeric(data[,j])}
  }
  data
}

##### settings: load model settings and data into global env (used for debugging)
settings <- function(dat = NULL, moderators = NULL, lags = NULL, d = FALSE){
  whatData <- c("sur1", "sur2", "sur3", "mvar", "mgm", "gauss1", "gauss2", 
                "autism", "large", "symptom", "short", "resting", "big5", 
                "esm", "esm1", "esm2", "esm3", "esm4", "esm5", "m3", "m5", 
                "mod", "PTSD", "fried", "wichers", "new", "dep1", "dep2", 
                "bfi1", "bfi2", "bfi3", "obama", "bfiDat", "constantini",
                "covid", "covid2", "esm6", "esm7")
  if(d){return(data.frame(sort(whatData)))}
  setwd('~/C: Desktop/COMPS/METHODS/CODE/modnets')
  if(!is.null(dat)){
    dat <- match.arg(dat, c(whatData, "caseDrop"))
    if("caseDrop" %in% c(dat, moderators)){
      poop <- list(nboots = 10, m = NULL, threshold = FALSE, ci = .95, 
                   verbose = FALSE, caseDrop = TRUE, caseMin = .05, 
                   caseMax = .75, caseN = 10, saveMods = FALSE)
      list2env(poop, envir = .GlobalEnv)
      if(dat == "caseDrop"){return(message("settings loaded for case dropping"))}
    }
    if(grepl("sur", dat)){
      if(dat == "sur1"){
        S <- matrix(c(1, .4, .3, .4, 1, .6, .3, .6, 1), 3, 3)
        beta <- matrix(c(.01, .014, 0, .005, .07, .3, .01, .09, .1, .23, .3, .001), 
                       ncol = 4, nrow = 3, byrow = TRUE)
        beta2 <- matrix(c(.01, .05, .1, .1, .01, .005), ncol = 2, byrow = T)
        poop <- list(B = cbind(beta, beta2), S = S)
        poop$data <- SURsampler(B = poop$B, S = poop$S, n = 1000, seed = 666, allDat = FALSE)
        list2env(poop, envir = .GlobalEnv)
        return(message("B, S, and data in .GlobalEnv"))
      }
      if(dat == "sur2"){
        B <- list()
        B[[1]] <- c(0, .3, 0, .25, 0, 0, .07, 0)
        B[[2]] <- c(0, .4, 0, 0, .15, 0, 0, .09)
        B[[3]] <- c(0, .18, .27, 0, 0, .08, 0, 0)
        B[[4]] <- c(0, .5, 0, .31, 0, 0, .1, 0)
        B <- do.call(rbind, B)
        S <- matrix(c(1, .3, .7, .2, .3, 1, 0, .06, .7, 0, 1, .4, .2, .06, .4, 1), 
                    ncol = 4, byrow = T)
        poop <- list(B = B, S = S)
        poop$data <- SURsampler(B = poop$B, S = poop$S, n = 1000, seed = 666, allDat = FALSE)
        list2env(poop, envir = .GlobalEnv)
        return(message("B, S, and data in .GlobalEnv"))
      }
      if(dat == "sur3"){
        B <- list()
        B[[1]] <- c(0, .4, 0, .25, 0, 0, .07, 0)
        B[[2]] <- c(0, .3, 0, 0, .15, 0, 0, .09)
        B[[3]] <- c(0, .35, .27, 0, 0, .08, 0, 0)
        B[[4]] <- c(0, .5, 0, .31, 0, 0, .1, 0)
        B <- do.call(rbind, B)
        S <- matrix(c(1, .3, .7, .2, .3, 1, 0, .06, .7, 0, 1, .4, .2, .06, .4, 1), 
                    ncol = 4, byrow = T)
        poop <- list(B = B, S = S)
        poop$data <- SURsampler(B = poop$B, S = poop$S, n = 1000, seed = 666, allDat = FALSE)
        list2env(poop, envir = .GlobalEnv)
        return(message("B, S, and data in .GlobalEnv"))
      }
    }
    if(grepl("dep", dat)){
      if(dat == "dep1"){data <- na.omit(psych::msq[, c("hostile", "lonely", "nervous", "sleepy", "depressed")])}
      if(dat == "dep2"){data <- msq_p5}
      data0 <- data.frame(apply(data, 2, function(z) ifelse(z == 0, 0, 1)))
      poop <- list(data = data, data0 = data0)
      list2env(poop, envir = .GlobalEnv)
      return(message("data and data0 in .GlobalEnv"))
    }
    if(grepl('covid', dat)){
      data <- readRDS('data/covid.RDS')
      if(dat == 'covid2'){
        vs <- colnames(data)
        k <- attr(data, 'k')
        m <- vs[k][startsWith(vs[k], 'C')]
        data <- data[, setdiff(colnames(data), setdiff(vs[k], m))]
        attr(data, 'k') <- which(colnames(data) %in% m)
        attr(data, 'args') <- list(m = attr(data, 'k'), beepno = 'beepvar', dayno = 'dayvar')
      }
      assign('data', data, envir = .GlobalEnv)
      return(message('data in .GlobalEnv'))
    }
    if(dat == "constantini"){
      data0 <- readRDS("data/data_all.RDS")
      data1 <- readRDS("data/data_between.RDS")
      data2 <- readRDS("data/data_within.RDS")
      poop <- list(data0 = data0, data1 = data1, data2 = data2)
      list2env(poop, .GlobalEnv)
      return(message("data0, data1, and data2 in .GlobalEnv"))
    }
    if(dat == "mvar"){
      data <- mvar_data$data
      type <- c(rep("c", 4), rep("g", 2))
      level <- c(2, 2, 4, 4, 1, 1)
    }
    if(dat == "mgm"){
      data <- mgm_data$data
      type <- c("g", "c", "c", "g")
      level <- c(1, 2, 4, 1)
    }
    if(grepl("gauss", dat)){
      if(dat == "gauss1"){data <- readRDS("data/datG.RDS")$data}
      if(dat == "gauss2"){data <- readRDS("data/dataG.RDS")$data}
      type <- rep("g", 10)
      level <- rep(1, 10)
    }
    if(dat == "autism"){
      data <- autism_data$data[,-4]
      type <- autism_data$type[-4]
      level <- autism_data$lev[-4]
    }
    if(dat == "large"){
      data <- autism_data_large$data[,-c(8:13,15,18:21)]
      type <- autism_data_large$type[-c(8:13,15,18:21)]
      level <- autism_data_large$level[-c(8:13,15,18:21)]
    }
    if(grepl("symptom|short", dat)){
      if(dat == "symptom"){
        data <- symptom_data$data
        type <- symptom_data$type
        level <- symptom_data$level
      }
      if(dat == "short"){
        data <- symptom_data$data
        sel <- c(which(symptom_data$groups == "Mood"), 
                 which(symptom_data$colnames %in% c("Action", "Who with")))
        data <- data.frame(as.matrix(data[,sel]))
        type <- symptom_data$type[sel]
        level <- symptom_data$level[sel]
        colnames(data) <- symptom_data$colnames[sel]
        colnames(data)[13] <- "Who"
      }
      if(dat == "short2"){
        data <- symptom_data$data
        colnames(data) <- symptom_data$colnames
        wichrz = c(1, 3, 4, 7, 8, 9, 13:15)
        data <- data[, wichrz]
      }
      assign("data", data, envir = .GlobalEnv)
      return(message("data in .GlobalEnv"))
    }
    if(dat == "resting"){
      data <- restingstate_data$data
      type <- rep("g", 68)
      level <- rep(1, 68)
    }
    if(dat == "big5"){
      data <- B5MS
      type <- rep("g", 5)
      level <- rep(1, 5)
    }
    if(dat == "m3"){
      data <- msq_p3
      type <- rep("g", 3)
      level <- rep(1, 3)
    }
    if(dat == "m5"){
      data <- msq_p5
      type <- rep("g", 5)
      level <- rep(1, 5)
    }
    if(dat == "mod"){
      data <- modnw
      type <- rep("g", 13)
      level <- rep(1, 13)
    }
    if(dat == "PTSD"){
      data <- PTSD_data$data
      type <- PTSD_data$type
      level <- PTSD_data$level
    }
    if(dat == "fried"){
      data <- readRDS("data/Fried2015_nD.RDS")$data
      colnames(data) <- readRDS("data/Fried2015_nD.RDS")$names
      data <- data.frame(data)
      assign("data", data, envir = .GlobalEnv)
      return(message("data in .GlobalEnv"))
    }
    if(dat == "obama"){
      path <- "data/"
      file <- paste0(path, "anes_timeseries_2012_Stata12.dta")
      ANES2012 <- foreign::read.dta(file)
      ObamaCog <- data.frame(Mor = as.numeric(ANES2012$ctrait_dpcmoral),
                             Led = as.numeric(ANES2012 $ ctrait_dpclead),
                             Car = as.numeric(ANES2012$ctrait_dpccare),
                             Kno = as.numeric(ANES2012$ctrait_dpcknow),
                             Int = as.numeric(ANES2012$ctrait_dpcint),
                             Hns = as.numeric(ANES2012$ctrait_dpchonst))
      ObamaCog[ObamaCog < 3] <- NA
      ObamaCog2 <- suppressWarnings(bootnet::binarize(ObamaCog, 5, removeNArows = FALSE))
      ObamaAff <- data.frame(Ang = as.numeric(ANES2012$candaff_angdpc),
                             Hop = as.numeric(ANES2012$candaff_hpdpc), 
                             Afr = as.numeric(ANES2012$candaff_afrdpc), 
                             Prd = as.numeric(ANES2012$candaff_prddpc))
      ObamaAff[ObamaAff < 3] <- NA
      ObamaAff2 <- suppressWarnings(bootnet::binarize(ObamaAff, 4, removeNArows = FALSE))
      RomneyCog <- data.frame(Mor = as.numeric(ANES2012$ctrait_rpcmoral),
                              Led = as.numeric(ANES2012 $ ctrait_rpclead),
                              Car = as.numeric(ANES2012$ctrait_rpccare),
                              Kno = as.numeric(ANES2012$ctrait_rpcknow),
                              Int = as.numeric(ANES2012$ctrait_rpcint),
                              Hns = as.numeric(ANES2012$ctrait_rpchonst))
      RomneyCog[RomneyCog < 3] <- NA
      RomneyCog2 <- suppressWarnings(bootnet::binarize(RomneyCog, 5, removeNArows = FALSE))
      RomneyAff <- data.frame(Ang = as.numeric(ANES2012$candaff_angrpc),
                              Hop = as.numeric(ANES2012$candaff_hprpc), 
                              Afr = as.numeric(ANES2012$candaff_afrrpc), 
                              Prd = as.numeric(ANES2012$candaff_prdrpc))
      RomneyAff[RomneyAff < 3] <- NA
      RomneyAff2 <- suppressWarnings(bootnet::binarize(RomneyAff, 4, removeNArows = FALSE))
      ObaRom <- data.frame(ObamaCog, ObamaAff, RomneyCog, RomneyAff)
      ObaRom2 <- data.frame(ObamaCog2, ObamaAff2, RomneyCog2, RomneyAff2)
      ObaRom <- na.omit(ObaRom)
      ObaRom2 <- na.omit(ObaRom2)
      Obama <- ObaRom[, 1:10]
      Obama2 <- ObaRom2[, 1:10]
      Romney <- ObaRom[, 11:20]
      Romney2 <- ObaRom2[, 11:20]
      colnames(Romney) <- colnames(Obama)
      colnames(Romney2) <- colnames(Obama2)
      data <- data.frame(rbind(Obama, Romney))
      data0 <- data.frame(rbind(Obama2, Romney2))
      data$person <- rep(c(0, 1), each = 5509)
      data0$person <- rep(c(0, 1), each = 5509)
      rownames(data) <- rownames(data0) <- 1:nrow(data)
      list2env(list(data = data, data0 = data0), envir = .GlobalEnv)
      return(message("data and data0 in .GlobalEnv"))
    }
    if(dat == "wichers"){
      data <- readRDS("data/Wicherts2016_Mood.RDS")$data_mood
      colnames(data) <- readRDS("data/Wicherts2016_Mood.RDS")$labels
      type <- rep("g", 9)
      level <- rep(1, 9)
    }
    if(grepl('esm', dat)){
      data <- readRDS('data/ESMdata.RDS')
      timevars <- c(attr(data, 'timevars'), 'period')
      if(dat == 'esm1'){
        data <- data[, c(attr(data, 'useVars'), timevars)]
      } else if(dat %in% c('esm2', 'esm3', 'esm4')){
        data <- data[, c(attr(data, 'modVars'), timevars)]
        if(dat == 'esm3'){data <- data[, setdiff(colnames(data), setdiff(timevars, 'period'))]}
        if(dat == 'esm4'){data <- data[, setdiff(colnames(data), timevars)]}
      } else if(dat %in% c('esm5', 'esm6', 'esm7')){
        k <- switch(2 - (dat != 'esm5'), c('concentrat', 'period'), 'concentrat')
        data <- data[, c(attr(data, 'modVars'), k, setdiff(timevars, 'period'))]
        if(dat != 'esm5'){data <- na.omit(data[, 1:18]); levels(data$period) <- 1:6}
        if(dat == 'esm7'){data <- subset(data, period == 3)[, setdiff(names(data), 'period')]}
      }
      assign('data', data, envir = .GlobalEnv)
      return(message('data in .GlobalEnv'))
    }
    if(dat == "new"){
      data <- readRDS("data/newData.RDS")
      type <- data$call$type
      level <- data$call$level
      data <- data$data
    }
    if(grepl("bfi", dat)){
      if(dat == "bfiDat"){
        data <- readRDS("data/bfiDat.RDS")
        assign("data", data, envir = .GlobalEnv)
        return(message("data in .GlobalEnv"))
      }
      data <- psych::bfi
      aceno <- c("A", "C", "E", "N", "O")
      if(dat == "bfi1"){
        data <- na.omit(data[, 1:26])
        data[, "gender"] <- data[, "gender"] - 1
      }
      if(dat == "bfi2"){data <- na.omit(data[, c(1:25, 27)])}
      if(dat == "bfi3"){data <- na.omit(data[, c(1:25, 28)])}
      rv <- c("A1", "C4", "C5", "E1", "E2", "O2", "O5")
      data[, rv] <- (data[, rv] * (-1)) + 7
      data2 <- do.call(cbind, lapply(aceno, function(z) rowMeans(data[, grep(z, colnames(data))])))
      datNames <- c(aceno, colnames(data)[26])
      data <- data.frame(data2, data[, 26])
      colnames(data) <- datNames
      rownames(data) <- 1:nrow(data)
      mod <- list(data[, 6])
      names(mod) <- colnames(data)[6]
      dat <- list(Y = data[, 1:5], X = data)
      data <- data[, 1:5]
      list2env(list(data = data, mod = mod, dat = dat), envir = .GlobalEnv)
      return(message("data, mod, and dat in .GlobalEnv"))
    }
    poop <- list(data = data, type = type, level = level)
    if(all(type == "g") & all(level == 1)){poop <- list(data = data)}
  } else {
    poop <- list(lambda = NULL, lags = lags, rule = "AND", threeWay = TRUE, std = TRUE,
                 which.lam = "lambda.min", moderators = moderators, gamma = 0.25, 
                 scale = FALSE, seed = 666, family = NULL, folds = 10, threshold = TRUE,
                 measure = "deviance", saveData = TRUE, alpha = 1, group = NULL,
                 adaMethod = "ridge", adaGam = NULL, grPenalty = "grLasso",
                 time = FALSE, covariates = NULL, type = NULL, center = TRUE, 
                 verbose = TRUE, exogenous = NULL)
  }
  list2env(poop, envir = .GlobalEnv)
}

##### pars: load default arguments of functions to global environment
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
##### fitNetwork: Top-level function that ties everything together
fitNetwork <- function(data, moderators = NULL, type = "gaussian", lags = NULL,
                       seed = NULL, lambda = NULL, folds = 10, gamma = 0.5,
                       which.lam = 'lambda.min', rule = "OR", threshold = FALSE,
                       scale = FALSE, std = TRUE, group = NULL, grPenalty = "gel",
                       adaGam = NULL, adaMethod = "ridge", measure = "deviance",
                       alpha = 1, saveData = TRUE, center = TRUE, covariates = NULL, 
                       verbose = FALSE, exogenous = TRUE, binary = NULL, mval = NULL,
                       residMat = "sigma", medges = 1, pcor = FALSE, maxiter = 100,
                       getLL = TRUE, saveMods = TRUE, binarize = FALSE, fitCoefs = FALSE, 
                       detrend = FALSE, beepno = NULL, dayno = NULL, ...){
  t1 <- Sys.time() # START
  fullcall <- match.call()
  if(!identical(detrend, FALSE) | (!is.null(beepno) & !is.null(dayno))){lags <- 1}
  if(any(is.na(data))){
    ww <- which(apply(data, 1, function(z) any(is.na(z))))
    if(is.null(lags) | identical(lags, 0)){
      data <- data[-ww, ]
      warning(paste0(length(ww), ' cases deleted due to missingness'))
    } else {
      stop(paste0(length(ww), ' rows contain missing values'))
    }
  }
  args <- tryCatch({list(...)}, error = function(e){list()})
  if(identical(type, 'varSelect')){
    vargs <- list(data = data, m = moderators, lags = lags, exogenous = exogenous,
                  center = center, scale = scale, gamma = gamma, verbose = verbose)
    vargs$method <- ifelse(!is.null(moderators), 'glinternet', 
                           ifelse('method' %in% names(args), args$method, 'glmnet'))
    otherargs <- c('criterion', 'nfolds', 'varSeed', 'useSE', 'nlam')
    if(any(otherargs %in% names(args))){
      vargs <- append(vargs, args[intersect(names(args), otherargs)])
    }
    type <- do.call(varSelect, vargs)
  }
  if(!is.null(lags)){if(all(lags == 0)){lags <- NULL}}
  if(!is.null(lags)){
    if(lags != FALSE){lags <- 1}
    if("samp_ind" %in% names(attributes(data))){samp_ind <- attr(data, "samp_ind")}
    data <- data.frame(data)
    if(exists("samp_ind", inherits = FALSE)){attr(data, "samp_ind") <- samp_ind}
  } else {
    data <- data.frame(data)
  }
  if(!is.null(moderators)){
    if(all(moderators == 0)){moderators <- NULL}
    if(all(moderators == "all")){moderators <- 1:ncol(data)}
  }
  if(ifelse(!is.null(covariates), ifelse(!is(covariates, 'list'), TRUE, FALSE), FALSE)){
    if(length(covariates) >= ncol(data) - 1){stop("Must have at least 2 outcome variables")}
  }
  if(is.null(lambda) & ncol(data) >= nrow(data)){lambda <- "EBIC"}
  if(is.null(lambda)){ ### SETUP
    if(grepl("min", which.lam)){which.lam <- "min"} else {which.lam <- "1se"}
    output <- list()
    output$call <- list(type = type, moderators = moderators, mval = mval, lags = lags,
                        which.lam = which.lam, rule = rule, threshold = threshold,
                        center = center, scale = scale)
    if(length(args) != 0){output$call <- append(output$call, args)}
    if(is.null(mval)){output$call$mval <- NULL}
    if(class(type) == "list"){
      output$call$type <- lapply(type, function(z) unlist(z[[ifelse(which.lam == "min", 1, 2)]]))
      if(attr(type, "criterion") != "CV"){output$call$which.lam <- attr(type, "criterion")}
      if(attr(type, "criterion") == "EBIC"){output$call$gamma <- attr(type, "gamma")}
    } else {
      output$call$which.lam <- NULL
      if(length(type) == 1){
        type <- rep("gaussian", ncol(data))
        bb <- unname(which(apply(data, 2, function(z) dim(table(z)) <= 2)))
        if(length(bb) > 0){type[bb] <- "binomial"}
      }
      if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
      if("binomial" %in% type){
        type[type == "binomial"] <- "c"
        if(all(type == "c")){output$call[c("center", "scale")] <- center <- scale <- FALSE}
      }
      output$call$type <- type
    }
  }
  if(is.null(lambda)){
    ### CROSS-SECTIONAL
    if(is.null(lags)){
      output$call["lags"] <- NULL
      if(is.character(threshold)){output$call$threshold <- threshold <- FALSE}
      if(!is.null(covariates)){
        if(class(covariates) == "list"){covs <- names(covariates)}
        if(class(covariates) %in% c("numeric", "integer")){covs <- colnames(data)[covariates]}
        output$call$covariates <- covs
      }
      if(!is.null(moderators)){
        output$call$moderators <- ifelse(is.list(moderators), list(names(moderators)), 
                                         list(colnames(data)[moderators]))[[1]]
        if("pcor" %in% names(output$call)){output$call$pcor <- NULL}
      } else if("pcor" %in% names(output$call)){
        if(pcor != FALSE){
          if(threshold == FALSE){output$call$pcor <- TRUE}
          output$call$rule <- NULL
        }
      }
      if(!is.null(covariates) & !is.null(moderators)){
        if(!is.list(covariates) & !is.list(moderators)){
          if(max(moderators) > min(covariates)){
            if(is.character(type)){output$call$type <- type <- type[-covariates]}
            covs <- list(data[, covariates])
            data <- data[, -covariates]
            names(covs) <- output$call$covariates
            covariates <- covs
            moderators <- which(colnames(data) %in% output$call$moderators)
          }
        }
      }
      mods0 <- nodewise(data = data, mods = moderators, varMods = type, 
                        lambda = which.lam, center = center, scale = scale,
                        covariates = covariates, exogenous = exogenous)
      if(any(sapply(mods0$models, function(z) length(coef(z))) >= nrow(data)) & verbose){
        warning('Model is overspecified; not enough cases to estimate all parameters')
      }
      out <- modNet(models = mods0, threshold = threshold, rule = rule, 
                    mval = mval, pcor = pcor, binarize = binarize)
      output <- append(output, out)
      output$mods0 <- mods0
      if("moderator" %in% names(attributes(out))){
        attr(output, "moderator") <- attr(out, "moderator")
        output$call$exogenous <- attr(mods0$models, "exogenous")
      }
      if("mval" %in% names(attributes(out))){attr(output, "mval") <- attr(out, "mval")}
      attributes(output)$ggm <- TRUE
      if(getLL){
        output$modLL <- tryCatch({modLL(output, all = TRUE)}, error = function(e){list()})
        if(length(output$modLL) == 0){
          output$modLL <- NULL
        } else {
          dd <- output$data
          output[c('data', 'mods0')] <- NULL
          output$data <- dd
          output$mods0 <- mods0
        }
      }
      class(output) <- c('list', 'ggm')
      if(fitCoefs | !saveMods){
        output$mods0$models <- NULL
        output$fitobj <- if(fitCoefs){getFitCIs(output)} else {NULL}
      }
      attr(output, 'fullcall') <- fullcall
      if(verbose){print(Sys.time() - t1)}
      return(output)
    } else {
      ### TEMPORAL
      if(!is.null(beepno) & !is.null(dayno)){
        beepday <- list(beepno, dayno)
        stopifnot(sum(sapply(c(1, nrow(data)), function(z) all(sapply(beepday, length) == z))) == 1)
        if(all(sapply(beepday, length) == 1)){
          if(is.character(beepno)){beepno <- which(colnames(data) == beepno)}
          if(is.character(dayno)){dayno <- which(colnames(data) == dayno)}
          data0 <- data[, -c(beepno, dayno)]
          beepno <- data[, beepno]
          dayno <- data[, dayno]
          data <- data0
          if(exists("samp_ind", inherits = FALSE)){
            attr(data, "samp_ind") <- samp_ind
          }
        }
      } else if(!is.null(beepno) & is.null(dayno) | !is.null(dayno) & is.null(beepno)){
        stop('Must specify both beepno AND dayno, or neither')
      }
      if(!identical(detrend, FALSE)){
        if(is(detrend, 'list')){
          stopifnot(length(detrend) %in% 1:2)
          if(is.null(names(detrend))){
            if(length(detrend) == 1){
              detrend <- detrend[[1]]
            } else {
              names(detrend)[order(sapply(detrend, length))] <- c('timevar', 'vars')
            }
          }
        }
        timevar <- switch(2 - isTRUE(detrend), NULL, ifelse(
          is(detrend, 'list'), detrend$timevar, ifelse(
            is(detrend, 'numeric'), colnames(data)[detrend], detrend)))
        dvars <- switch(2 - is(detrend, 'list'), detrend$vars, NULL)
        data <- detrender(data = data, timevar = timevar, vars = dvars, verbose = verbose)
        if(exists("samp_ind", inherits = FALSE)){
          attr(data, "samp_ind") <- samp_ind
        }
      }
      if(!is.null(beepno) & !is.null(dayno)){
        consec <- mgm:::beepday2consec(beepvar = beepno, dayvar = dayno)
        consec <- mgm:::lagData(data = data, lags = 1, consec = consec)[[3]][-1]
        output$call$consec <- which(consec)
      } else {
        consec <- NULL
      }
      output$call$rule <- NULL
      if(threshold != FALSE){output$call$pcor <- ifelse(is.logical(pcor), "none", pcor)}
      if(class(type) == "list"){
        if(length(type) == 1){stop("Need more than one outcome to construct network")}
        if(!attr(type, "method") %in% c("regsubsets", "glmnet")){
          exogenous <- attr(type, "exogenous")
          moderators <- attr(type, "moderators")
          if(!is.character(moderators)){
            attr(type, "moderators") <- moderators <- colnames(data)[moderators]
          }
          anymods <- any(sapply(moderators, function(z){
            any(grepl(paste0(z, ":|:", z), unlist(output$call$type)))
          }))
          moderators <- which(colnames(data) %in% moderators)
          if("covs" %in% names(attributes(type))){covariates <- attr(type, "covs")}
        }
      }
      if(!is.null(union(moderators, covariates))){
        mco <- sort(union(moderators, covariates))
        stopifnot(identical(mco, sort(c(moderators, covariates))))
        if(exogenous & length(mco) >= ncol(data) - 1){exogenous <- FALSE}
      }
      if(!is.null(moderators)){
        output$call$moderators <- colnames(data)[moderators]
        output$call$exogenous <- exogenous
      }
      if(exists("anymods", inherits = FALSE)){
        if(!anymods){
          covariates <- union(covariates, moderators)
          moderators <- NULL
        }
      }
      fit <- SURfit(data = data, varMods = type, m = moderators, mod = which.lam, 
                    center = center, scale = scale, exogenous = exogenous, 
                    covs = covariates, sur = std, maxiter = maxiter, consec = consec)
      dat <- lagMat(data = data, type = type, m = moderators, covariates = covariates, 
                    center = center, scale = scale, exogenous = exogenous, consec = consec)
      net <- SURnet(fit = fit, dat = dat, s = residMat, m = moderators, pcor = pcor,
                    threshold = threshold, mval = mval, medges = medges)
      if(!is.null(covariates)){
        output$call$covariates <- net$call$covariates <- colnames(data)[covariates]}
      output$SURnet <- net
      attributes(output$SURnet)$SURnet <- TRUE
      output$SURfit <- fit
      if("mnet" %in% names(net)){
        attr(output, "mnet") <- net$call$moderators
        if("mval" %in% names(net$call)){attr(output, "mval") <- net$call$mval}
      }
      attributes(output)$SURnet <- TRUE
      if(getLL){
        output$SURll <- tryCatch({SURll(output, all = TRUE, s = residMat)}, error = function(e){list()})
        if(length(output$SURll) == 0){output$SURll <- NULL}
      }
      attr(output, 'rank') <- sum(sapply(lapply(output$SURnet$mods, '[[', 'model'), nrow))
      class(output) <- c('list', 'SURnet')
      if(fitCoefs | !saveMods){
        output$SURfit <- if(fitCoefs){getFitCIs(output)} else {NULL}
      }
      attr(output, 'fullcall') <- fullcall
      if(verbose){print(Sys.time() - t1)}
      return(output)
    }
  }
  ### PENALIZED MODELS
  if(!is.null(moderators)){
    if(class(moderators) %in% c("numeric", "integer")){
      vnames <- colnames(data)[-moderators]
      mnames <- colnames(data)[moderators]
      data <- data.frame(data[, -moderators], data[, moderators])
      colnames(data) <- c(vnames, mnames)
    } else if(is.list(moderators)){
      data <- data.frame(data, moderators)
    }
    moderators <- ncol(data)
  } else if("glinternet" %in% lambda){
    moderators <- 1:ncol(data)
  }
  if(is.null(type)){type <- "gaussian"}
  if(length(type) == 1){type <- match.arg(type, c("gaussian", "binomial", "multinomial"))}
  output <- list()
  output$call <- list(type = type, moderators = moderators, lags = lags, 
                      lambda = lambda, gamma = gamma, folds = folds, seed = seed, 
                      rule = rule, group = group, adaGam = adaGam, alpha = alpha, 
                      which.lam = which.lam, measure = measure, scale = scale, 
                      adaMethod = adaMethod, grPenalty = grPenalty, 
                      threshold = threshold, std = std)
  if(all(type == "gaussian")){type <- rep("g", ncol(data))}
  if(!all(lambda %in% c("CV", "glinternet")) & is.null(adaGam)){
    output$call[c("seed", "which.lam", "folds", "measure")] <- NULL}
  if(!is.null(lags)){output$call$rule <- NULL} else {output$call$lags <- NULL}
  if(is.null(adaGam)){output$call$adaGam <- output$call$adaMethod <- NULL}
  if(is.null(group)){output$call$group <- output$call$grPenalty <- NULL}
  if(is.null(moderators)){output$call$moderators <- NULL}
  if("AIC" %in% lambda){output$call$gamma <- NULL}
  if(verbose){pb <- txtProgressBar(min = 0, max = ncol(data), style = 3)}
  mods <- lapply(1:ncol(data), function(x){
    nodes <- nodeMods(data = data, y = x, type = type, lambda = lambda, which.lam = which.lam, 
                      folds = folds, lags = lags, gamma = gamma, seed = seed, alpha = alpha,
                      moderators = moderators, scale = scale, measure = measure, group = group,
                      adaMethod = adaMethod, adaGam = adaGam, grPenalty = grPenalty, std = std)
    if(verbose){setTxtProgressBar(pb, x)}
    return(nodes)
  })
  if(verbose){close(pb)}
  adjMats <- combineMods(data = data, mods = mods, type = type, scale = scale, 
                         threshold = threshold, moderators = moderators, 
                         rule = rule, lags = lags)
  if(!is.null(lags) & length(lags) > 1){
    if(!is.null(moderators)){
      output$adjMats <- lapply(1:length(lags), function(z) adjMats[[z]][[1]]$adjMat)
      output$edgeColors <- lapply(1:length(lags), function(z) adjMats[[z]][[1]]$colMat)
      output$interactions <- lapply(1:length(lags), function(z) adjMats[[z]][[2]])
      n <- 4
    } else {
      output$adjMats <- lapply(1:length(lags), function(z) adjMats[[z]]$adjMat)
      output$edgeColors <- lapply(1:length(lags), function(z) adjMats[[z]]$colMat)
      n <- 3
    }
    for(i in 2:n){names(output[[i]]) <- paste0("lag", lags)}
  } else {
    if(!is.null(moderators)){
      output$adjMat <- adjMats[[1]]$adjMat
      output$edgeColors <- adjMats[[1]]$colMat
      if(is.null(lags)){
        output$nodewise <- list(adjNW = adjMats[[1]]$nodewise, edgeColsNW = adjMats[[1]]$nwCols)
        absMeans <- do.call(rbind, lapply(lapply(adjMats[[2]], abs), mean))
        trueMeans <- do.call(rbind, lapply(adjMats[[2]], mean))
        output$interactions <- list(absoluteMeans = absMeans, trueMeans = trueMeans, coefs = adjMats[[2]])
      } else {
        output$interactions <- adjMats[[2]]
      }
    } else {
      output$adjMat <- adjMats$adjMat
      output$edgeColors <- adjMats$colMat
      output$nodewise <- list(adjNW = adjMats$nodewise, edgeColsNW = adjMats$nwCols)
    }
  }
  output$mods <- lapply(mods, function(z) z[-which(names(z) == "fitobj")])
  output$fitobj <- lapply(mods, function(z) z[["fitobj"]])
  if(all(lambda != "CV")){output$fitobj <- lapply(output$fitobj, function(z) z$fit)}
  if(saveData == TRUE){output$data <- data}
  if(verbose){print(Sys.time() - t1)}
  class(output) <- c('list', ifelse(is.null(lags) | identical(as.numeric(lags), 0),'ggm', 'SURnet'))
  output
}

##### plotNet: Plot model results
plotNet <- function(object, which.net = 'temporal', threshold = FALSE, layout = 'spring', 
                    predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                    scale = FALSE, lag = NULL, con = 'adjR2', cat = 'nCC', 
                    plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                    binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    if(all(adjMat %in% 0:1)){colMat[colMat != 'darkgrey'] <- 'darkgrey'}
    dimnames(colMat) <- dimnames(adjMat)
    colMat
  }
  if(is.logical(which.net)){threshold <- TRUE; which.net <- 'temporal'}
  if(any(class(object) %in% c('qgraph', 'bootnetResult', 'bootnet'))){
    if(is(object, 'bootnet')){object <- object$sample}
    if(is(object, 'bootnetResult')){
      if(object$default == 'graphicalVAR'){
        xx <- unname(which(sapply(c('t', 'c', 'b', 'p'), function(z) startsWith(which.net, z))))
        if(length(xx) == 0){stop('Must specify which.net')}
        which.net <- switch(xx, 'beta', 'kappa', 'beta', match.arg(toupper(which.net), c('PCC', 'PDC')))
        object <- object$results[[match(which.net, names(object$results))]]
        if(which.net == 'beta'){object <- object[, -1]}
      }
    }
    object <- getWmat(object)
  }
  if(is(object, 'bootNet') | (isTRUE(attr(object, 'resample')) & 'fit0' %in% names(object))){
    object <- object$fit0
  }
  if(is(object, 'mgmSim')){
    names(object)[grep('^trueNet$|^b1$', names(object))] <- 'adjMat'
    dimnames(object$adjMat) <- rep(list(1:ncol(object$adjMat)), 2)
    object$edgeColors <- getEdgeColors(object$adjMat)
    if('b2' %in% names(object)){object$modEdges <- 1 * (object$b2 != 0) + 1}
  }
  if(is(object, 'mgm')){object <- net(object, which.net, threshold)}
  atts <- names(attributes(object))
  if(!is(object, 'list') | 'simMLgvar' %in% atts){
    predict <- NULL; nodewise <- FALSE
    if(is(object, "mlVARsim") | 'simMLgvar' %in% atts){
      if('mm' %in% names(object) & which.net %in% c('temporal', 'beta')){
        nodewise <- TRUE
        modEdges <- t(1 * (object$mm$mb2 != 0) + 1)
      }
      object <- t(net(object, which.net))
    }
    stopifnot(ncol(object) == nrow(object))
    if(isTRUE(names)){names <- colnames(object)}
    object <- list(adjMat = object, edgeColors = getEdgeColors(object))
    if(nodewise){
      nodewise <- FALSE
      object$modEdges <- modEdges
    }
  }
  if(any(c("SURnet", "mlGVAR", "lmerVAR") %in% atts)){
    if(is.numeric(which.net)){which.net <- c("t", "c", "b", "pdc")[which.net]}
    which.net <- match.arg(tolower(which.net), c(
      "temporal", "contemporaneous", "between", "beta", "pdc", "pcc"))
    which.net <- switch(which.net, pcc = "contemporaneous", beta = "temporal", which.net)
    if(isTRUE(attr(object, "mlGVAR"))){
      object <- object[[switch(which.net, between = "betweenNet", "fixedNets")]]
    }
    if("SURnet" %in% names(object)){
      if(!is.null(mselect) & length(object$call$moderators) > 1 & isTRUE(object$call$exogenous)){
        if(isTRUE(mselect)){mselect <- object$call$moderators[1]}
        adjm <- netInts(fit = object, n = 'temporal', threshold = TRUE, 
                        avg = FALSE, rule = 'none', empty = FALSE, 
                        mselect = mselect)
        object$SURnet$temporal$modEdges <- abs(sign(adjm)) + 1
      }
      object <- object$SURnet
    }
    if(!"adjMat" %in% names(object)){
      object[c("adjMat", "edgeColors")] <- eval(parse(text = paste0("object$", switch(
        which.net, pdc = "temporal$PDC", which.net))))[c("adjMat", "edgeColors")]
      if(!startsWith(which.net, "c") & "modEdges" %in% names(object$temporal)){
        object$modEdges <- object$temporal$modEdges
      }
    }
  }
  obmods <- object$call$moderators
  exo <- ifelse('exogenous' %in% names(object$call), object$call$exogenous, TRUE)
  if(isTRUE(attr(object, 'ggm')) & ifelse(
    all(!sapply(list(obmods, mselect), is.null)), 
    length(obmods) > 1 & exo, FALSE)){
    if(isTRUE(mselect)){mselect <- obmods[1]}
    adjm <- netInts(fit = object, n = 'between', threshold = TRUE, 
                    avg = !nodewise, rule = rule, empty = FALSE, 
                    mselect = mselect)
    if(nodewise){
      object$nodewise$modEdgesNW <- abs(sign(adjm)) + 1
      diag(object$nodewise$modEdgesNW) <- 0
    } else {
      object$modEdges <- abs(sign(adjm)) + 1
      diag(object$modEdges) <- 0
    }
  }
  if(threshold != FALSE){
    if(!is.numeric(threshold)){threshold <- .05}
    if(mnet & "mnet" %in% names(object)){
      mn <- object$call$moderators
      adj1 <- net(object, "beta", threshold, rule)
      if(isTRUE(attr(object, "SURnet"))){
        rn <- gsub("[.]y$", "", rownames(object$mnet$adjMat))
        ind <- which(rn == mn)
        object$mnet$adjMat[-ind, -ind] <- adj1
        object$mnet$adjMat[-ind, ind] <- object$mnet$adjMat[-ind, ind] * ifelse(
          object$temporal$coefs$pvals[, mn] <= threshold, 1, 0)
      } else {
        ind <- nrow(object$mnet$adjMat)
        object$mnet$adjMat[-ind, -ind] <- adj1
        mps <- object$mods0$Bm
        mps2 <- matrix(0, nrow = nrow(adj1), ncol = 4)
        mps2[match(rownames(mps), rownames(adj1)), ] <- mps
        object$mnet$adjMat[ind, -ind] <- object$mnet$adjMat[ind, -ind] * ifelse(mps2[, 4] <= threshold, 1, 0)
        #object$mnet$adjMat[ind, -ind] <- object$mnet$adjMat[ind, -ind] * ifelse(object$mods0$Bm[, 4] <= threshold, 1, 0)
      }
    } else {
      if("ggm" %in% atts & isTRUE(nodewise)){
        object$nodewise$adjNW <- net(object, "beta", threshold, rule, TRUE)
      } else {
        object$adjMat <- net(object, which.net, threshold, rule)
      }
    }
  }
  if(isTRUE(names)){
    names <- gsub("[.]y$", "", rownames(object$adjMat))
    if(length(names) == 0){names <- 1:nrow(object$adjMat)}
  } else if(is.null(names) | all(names == FALSE)){
    names <- 1:nrow(object$adjMat)
  }
  names <- names[1:nrow(object$adjMat)]
  if(is.null(lag) & "adjMats" %in% names(object)){
    stop("More than one lag modeled; need to specify which to plot")
  } else if(!"adjMats" %in% names(object)){
    if(any(grepl("lag", colnames(object$adjMat)))){lag <- 1}
  }
  if(!is.null(predict) & !mnet & !identical(predict, FALSE)){
    con <- match.arg(con, choices = c("R2", "adjR2", "MSE", "RMSE"))
    cat <- match.arg(cat, choices = c("nCC", "CC", "CCmarg"))
    type <- object$call$type
    if("ggm" %in% names(attributes(object))){
      type <- unname(sapply(object$fitobj, attr, "family"))
      if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
      if("binomial" %in% type){type[type == "binomial"] <- "c"}
    }
    tt <- length(type)
    if(class(predict) == "list"){
      if(isTRUE(attr(predict, "mlGVAR"))){
        predict <- predict[[switch(which.net, between = "betweenNet", "fixedNets")]]
      }
      predict <- list(predictNet(predict, all = FALSE, scale = scale), 
                      predictNet(object, all = FALSE, scale = scale))
      stopifnot(length(predict) == 2)
      pie <- list(); pieColor <- list()
      for(i in 1:tt){
        if(type[i] == "g"){
          pie[[i]] <- c(predict[[1]][i,con][[1]], 
                        unlist(predict[[2]][i,con]) - unlist(predict[[1]][i,con]))
          pieColor[[i]] <- c("lightblue", "lightblue4")
          if(pie[[i]][1] < 0 & pie[[i]][2] < 0){
            pie[[i]][1] <- abs(pie[[i]][1])
            pie[[i]][2] <- abs(pie[[i]][2])
            pieColor[[i]] <- c("tomato", "tomato4")
          }
          if(pie[[i]][2] < 0){
            pie[[i]][1] <- pie[[i]][1] + pie[[i]][2]
            pie[[i]][2] <- abs(pie[[i]][2])
            pieColor[[i]][2] <- "tomato"
          }
          if(pie[[i]][1] < 0){
            pie[[i]][1] <- abs(pie[[i]][1])
            pie[[i]][2] <- pie[[i]][2] - pie[[i]][1]
            pieColor[[i]][1] <- "tomato"
          }
        }
        if(type[i] == "c"){
          pie[[i]] <- c(predict[[1]][i,cat][[1]], 
                        unlist(predict[[2]][i,cat]) - unlist(predict[[1]][i,cat]))
          pieColor[[i]] <- c("peachpuff", "peachpuff4")
          if(pie[[i]][2] < 0){
            pie[[i]][1] <- pie[[i]][1] + pie[[i]][2]
            pie[[i]][2] <- abs(pie[[i]][2])
            pieColor[[i]][2] <- "tomato"
          }
        }
      }
    } else {
      predict <- predictNet(object, all = FALSE, scale = scale)
      pie <- c(); pieColor <- c()
      for(i in 1:tt){
        if(type[i] == "g"){
          pie[i] <- predict[i,con]
          pieColor[i] <- "lightblue"
        }
        if(type[i] == "c"){
          pie[i] <- predict[i,cat]
          pieColor[i] <- "peachpuff"
        }
      }
      if(any(pie < 0)){
        pieColor[pie < 0] <- "tomato"
        pie[pie < 0] <- abs(pie[pie < 0][[1]])
      }
    }
  } else {
    pie <- NULL; pieColor <- NULL
  }
  if(!is.null(lag)){
    if(lag != 1 | lag == 1 & "adjMats" %in% names(object)){
      if(lag > length(object$adjMats)){
        lag <- which(sapply(lapply(lapply(object$adjMats, colnames), 
                                   function(z) grep(paste0("lag", lag), z)), length) != 0)
      }
      qgraph(input = t(object$adjMats[[lag]]), layout = layout, labels = names, 
             edge.color = t(object$edgeColors[[lag]]), edge.labels = elabs,
             edge.label.cex = elsize, DoNotPlot = !plot, pie = pie, 
             pieColor = pieColor, ...)
    } else if(!mnet){
      if("modEdges" %in% names(object) & mlty){lty <- t(object$modEdges)} else {lty <- 1}
      qgraph(input = t(object$adjMat), layout = layout, edge.color = t(object$edgeColors),
             labels = names, DoNotPlot = !plot, pie = pie, pieColor = pieColor, 
             edge.labels = elabs, lty = lty, edge.label.cex = elsize, ...)
    } else {
      if(class(names) %in% c("numeric", "integer")){
        names <- 1:ncol(object$mnet$adjMat)
      } else {
        names <- gsub(".lag1.", "", colnames(object$mnet$adjMat))
      }
      lty <- t(object$mnet$modEdges)
      qgraph(input = t(object$mnet$adjMat), layout = layout, labels = names,
             edge.color = t(object$mnet$edgeColors), DoNotPlot = !plot, 
             pie = pie, pieColor = pieColor, shape = object$mnet$shape, 
             edge.labels = elabs, lty = lty, edge.label.cex = elsize, ...)
    }
  } else {
    if(!nodewise){
      if(!mnet){
        if(binarize){object$adjMat[object$adjMat != 0] <- 1; object$edgeColors <- NA}
        if("modEdges" %in% names(object) & mlty){lty <- object$modEdges} else {lty <- 1}
        qgraph(input = object$adjMat, layout = layout, edge.color = object$edgeColors, 
               labels = names, DoNotPlot = !plot, pie = pie, pieColor = pieColor, 
               edge.labels = elabs, lty = lty, edge.label.cex = elsize, ...)
      } else {
        pp <- ncol(object$mnet$adjMat) - 1
        if(class(names) %in% c("numeric", "integer")){
          names <- 1:(pp + 1)
        } else {
          names <- c(names, colnames(object$mnet$adjMat)[pp + 1])
        }
        lty <- object$mnet$modEdges
        qgraph(input = object$mnet$adjMat, layout = layout, labels = names,
               edge.color = object$mnet$edgeColors, DoNotPlot = !plot, pie = pie,
               pieColor = pieColor, edge.labels = elabs, edge.label.cex = elsize, 
               lty = lty, directed = object$mnet$d, shape = c(rep("circle", pp), "square"), ...)
      }
    } else {
      if("modEdgesNW" %in% names(object$nodewise) & mlty){lty <- object$nodewise$modEdgesNW} else {lty <- 1}
      qgraph(input = object$nodewise$adjNW, layout = layout, labels = names,
             edge.color = object$nodewise$edgeColsNW, DoNotPlot = !plot, 
             pie = pie, pieColor = pieColor, edge.labels = elabs, 
             edge.label.cex = elsize, lty = lty, ...)
    }
  }
}

##### plotNet2: plot the temporal and contemporaneous networks in the same window
plotNet2 <- function(object, whichNets = NULL, whichTemp = c("temporal", "PDC"),
                     titles = c("PDC ", "PCC "), ...){
  whichTemp <- match.arg(whichTemp)
  if("lmerVAR" %in% names(attributes(object))){whichTemp <- "temporal"}
  if(!is.null(whichNets)){
    tt <- ifelse(all(tolower(whichNets) == "all"), 3, 2)
    if(tt == 3){whichNets <- c(whichTemp, "contemporaneous", "between")}
    if(tt == 2){
      l <- averageLayout(plotNet(object, which.net = whichNets[1], plot = FALSE),
                         plotNet(object, which.net = whichNets[2], plot = FALSE))
    } else if(tt == 3){
      l <- averageLayout(plotNet(object, which.net = whichNets[1], plot = FALSE),
                         plotNet(object, which.net = whichNets[2], plot = FALSE),
                         plotNet(object, which.net = whichNets[3], plot = FALSE))
      if(length(titles) == 2){titles <- c("Temporal Effects", "Contemporaneous Effects", "Between-Subject Effects")}
    }
  } else {
    tt <- 2
    l <- averageLayout(plotNet(object, which.net = whichTemp, plot = FALSE), 
                       plotNet(object, which.net = "contemporaneous", plot = FALSE))
  }
  layout(t(1:tt))
  if(tt == 2){
    if(all(titles == c("PDC ", "PCC "))){
      if("lmerVAR" %in% names(attributes(object)) & all(whichNets == c("temporal", "contemporaneous"))){
        title1 <- "Temporal Effects"
        title2 <- "Contemporaneous Effects"
        titles <- c(title1, title2)
      } else {
        title1 <- ifelse(whichTemp == "temporal", "Temporal Effects", "Partial Directed Correlations")
        title2 <- "Partial Contemporaneous Correlations"
      }
    } else {
      title1 <- titles[1]
      title2 <- titles[2]
    }
  }
  if(!is.null(whichNets)){
    if(length(titles) == 2){if(all(titles == c("PDC ", "PCC "))){titles <- whichNets}}
    for(i in 1:tt){
      plotNet(object, which.net = whichNets[i], layout = l, title = titles[i], ...)
    }
  } else {
    plotNet(object, which.net = whichTemp, layout = l, title = title1, ...)
    plotNet(object, which.net = "contemporaneous", layout = l, title = title2, ...)
  }
}

##### plotNet3: Designed for mlGVAR and lmerVAR output
plotNet3 <- function(object, ..., nets = c('temporal', 'contemporaneous', 'between'),
                     titles = TRUE, l = 3, label = NULL, xpos = 0, ypos = .5){
  args0 <- list(...)
  avlay(l = l)
  args0$object <- object
  stopifnot(length(nets) == 3)
  mains <- c('Temporal network', 'Partial Contemporaneous Correlations', 'Between-subjects network')
  invisible(lapply(1:3, function(i){
    args <- replace(args0, 'which.net', list(which.net = nets[i]))
    if(isTRUE(titles)){args$title <- mains[i]}
    do.call(plotNet, args)
  }))
  if(!is.null(label)){
    plot.new()
    text(x = xpos, y = ypos, labels = label)
  }
}

##### getEdgeColors
getEdgeColors <- function(adjMat){
  obj <- sign(as.vector(adjMat))
  colMat <- rep(NA, length(obj))
  if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
  if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
  if(any(obj == -1)){colMat[obj == -1] <- "red"}
  colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
  if(all(adjMat %in% 0:1)){colMat[colMat != 'darkgrey'] <- 'darkgrey'}
  dimnames(colMat) <- dimnames(adjMat)
  colMat
}

##### predictNet: Calculate prediction error
predictNet <- function(object, data = NULL, all = FALSE, scale = FALSE){
  if("SURnet" %in% c(names(object), names(attributes(object)))){
    if("SURnet" %in% names(object)){object <- object$SURnet}
    data <- object$data
    mods <- object$mods
    type <- object$call$type
    moderators <- object$call$moderators
    if(is.character(moderators)){moderators <- which(colnames(data$X) %in% moderators)}
    lags <- object$call$lags
    predObj <- Y <- list()
    for(i in 1:length(mods)){
      predObj[[i]] <- cbind(1, data$X) %*% mods[[i]]$model
      Y[[i]] <- as.numeric(data$Y[, i])
    }
    Y <- data.frame(do.call(cbind, Y))
    yhat <- data.frame(as.matrix(do.call(cbind, predObj)))
    colnames(Y) <- colnames(data$Y)
    colnames(yhat) <- paste0(colnames(data$Y), "hat")
    data <- data$Y
    probObj <- vector("list", length = ncol(data))
  } else {
    if(is.null(data)){
      x <- try(data <- object$data)
      if(class(x) == "try-error"){stop("Must supply dataset")}
    }
    mods <- object$mods
    type <- object$call$type
    moderators <- object$call$moderators
    if(is.character(moderators)){
      if("mods0" %in% names(object)){
        moderators <- which(colnames(object$mods0$dat) %in% moderators)
      } else {
        moderators <- which(colnames(data) %in% moderators)
      }
    }
    if("lags" %in% names(object$call)){lags <- object$call$lags} else {lags <- NULL}
    if(!is.null(moderators)){if(length(moderators) == 1){if(moderators == 0){moderators <- NULL}}}
    if(is.null(colnames(data))){colnames(data) <- paste0("V", 1:ncol(data))}
    predObj <- list(); Y <- list(); probObj <- vector("list", length = ncol(data))
    if("ggm" %in% names(attributes(object))){
      mods <- lapply(mods, function(z) list(model = z$model))
      type <- unname(sapply(object$fitobj, attr, "family"))
      if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
      if("binomial" %in% type){type[type == "binomial"] <- "c"}
      for(i in 1:length(mods)){
        if(type[i] == "c"){
          dat <- setup(data = object$mods0$dat, type = type, y = i, lags = lags, scale = scale)
          if(!is.null(lags)){y <- 1} else {y <- i}
          X <- interactionMatrix(data = dat, y = y, type = type, moderators = moderators, lags = lags)
          n <- nrow(X); p <- ncol(X)
          coefs <- list(mods[[i]]$model * -1, mods[[i]]$model)
          n_cats <- length(coefs)
          potentials <- matrix(NA, nrow(data), n_cats)
          for(cat in 1:n_cats){potentials[, cat] <- exp(cbind(1, X) %*% coefs[[cat]])[, 1]}
          probs <- potentials[, 1:n_cats]/rowSums(potentials[, 1:n_cats])
          probObj[[i]] <- probs
          predObj[[i]] <- sort(unique(dat[, y]))[apply(probs, 1, which.max)]
          Y[[i]] <- as.numeric(dat[, y])
        } else {
          predObj[[i]] <- predict(object = object$fitobj[[i]])
          Y[[i]] <- as.numeric(data[, i])
        }
      }
    } else {
      for(i in 1:length(mods)){
        dat <- setup(data = data, type = type, y = i, lags = lags, scale = scale)
        if(!is.null(lags)){y <- 1} else {y <- i}
        X <- interactionMatrix(data = dat, y = y, type = type, moderators = moderators, lags = lags)
        n <- nrow(X); p <- ncol(X)
        if(type[i] == "c"){
          coefs <- mods[[i]]$model
          n_cats <- length(coefs)
          potentials <- matrix(NA, n, n_cats)
          for(cat in 1:n_cats){potentials[,cat] <- exp(cbind(1, X) %*% coefs[[cat]])[,1]}
          probs <- potentials[,1:n_cats]/rowSums(potentials[,1:n_cats])
          probObj[[i]] <- probs
          predObj[[i]] <- sort(unique(dat[,y]))[apply(probs, 1, which.max)]
        } else {
          predObj[[i]] <- cbind(1, X) %*% mods[[i]]$model
        }
        Y[[i]] <- as.numeric(dat[,y])
      }
    }
    Y <- data.frame(do.call(cbind, Y))
    yhat <- data.frame(as.matrix(do.call(cbind, predObj)))
    colnames(Y) <- paste0(colnames(data), ".y")
    colnames(yhat) <- paste0(colnames(data), ".yhat")
    names(probObj) <- paste0(colnames(data), ".probs")
  }
  calcError <- function(y, yhat, type, mod){
    if(type == "g"){
      R2 <- 1 - sum((y - yhat)^2)/sum((y - mean(y))^2)
      adjR2 <- 1 - (1 - R2) * ((length(y) - 1)/(length(y) - sum(mod != 0)))
      MSE <- sum((y - yhat)^2)/(length(y) - sum(mod != 0))
      RMSE <- sqrt(MSE)
      predError <- list(R2 = R2, adjR2 = adjR2, MSE = MSE, RMSE = RMSE)
    }
    if(type == "c"){
      CC <- sum((y == yhat)/length(y))
      nCC <- (CC - max(table(y)/sum(table(y))))/(1 - max(table(y)/sum(table(y))))
      CCmarg <- (nCC - CC)/(nCC - 1)
      predError <- list(CC = CC, nCC = nCC, CCmarg = CCmarg)
    }
    lapply(predError, round, 3)
  }
  errors <- lapply(1:ncol(data), function(z) calcError(Y[,z], yhat[,z], type[z], mods[[z]]$model))
  if(all(c("g", "c") %in% type)){
    errors <- lapply(errors, function(z) do.call(cbind, z))
    for(i in 1:length(errors)){
      if("R2" %in% colnames(errors[[i]])){
        errors[[i]] <- c(errors[[i]], rep(NA, 3))
      } else {
        errors[[i]] <- c(rep(NA, 4), errors[[i]])
      }
    }
    errors <- data.frame(Variable = colnames(data), do.call(rbind, errors))
    colnames(errors)[2:8] <- c("R2", "adjR2", "MSE", "RMSE", "CC", "nCC", "CCmarg")
  } else {
    errors <- data.frame(Variable = colnames(data), do.call(rbind, errors))
  }
  if(all == TRUE){
    output <- list(Y = Y, preds = yhat, probs = probObj, errors = errors)
    if(all(type == "g")){output$probs <- NULL}
  } else {
    output <- data.frame(Variable = errors[, 1], apply(errors[, -1], 2, unlist))
  }
  return(output)
}

##### predictDelta: Calculate change in prediction for across nested models
predictDelta <- function(mod1, mod0, scale = FALSE){
  if("SURnet" %in% names(mod0)){
    stopifnot(ncol(mod0$SURnet$data$Y) == ncol(mod1$SURnet$data$Y))
    type <- rep("g", ncol(mod0$SURnet$data$Y))
    r0 <- attr(mod0, 'rank')
    r1 <- attr(mod1, 'rank')
    if(r1 < r0){mod00 <- mod1; mod1 <- mod0; mod0 <- mod00}
  } else {
    type <- mod1$call$type
    if("ggm" %in% names(attributes(mod1))){
      if(is.list(type)){
        type <- unname(sapply(mod1$fitobj, attr, "family"))
        if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
        if("binomial" %in% type){type[type == "binomial"] <- "c"}
      } else {
        type <- ifelse(length(type) == 1, ifelse(
          type %in% c("g", "gaussian"), rep("g", ncol(data)), rep("c", ncol(data))), type)
      }
    } else {
      stopifnot(mod1$call$type == mod0$call$type)
    }
  }
  errors <- list(predictNet(mod0, all = FALSE, scale = scale), predictNet(mod1, all = FALSE, scale = scale))
  p <- ifelse(all(type == "c"), 4, ifelse(all(type == "g"), 5, 8))
  delta <- data.frame(matrix(NA, nrow = length(type), ncol = p))
  for(i in 1:length(type)){
    if(type[i] == "g"){
      delta[i,2] <- unlist(errors[[2]][i, "R2"]) - unlist(errors[[1]][i, "R2"])
      delta[i,3] <- unlist(errors[[2]][i, "adjR2"]) - unlist(errors[[1]][i, "adjR2"])
      delta[i,4] <- unlist(errors[[2]][i, "MSE"]) - unlist(errors[[1]][i, "MSE"])
      delta[i,5] <- unlist(errors[[2]][i, "RMSE"]) - unlist(errors[[1]][i, "RMSE"])
    }
    if(type[i] == "c"){
      delta[i,(p-2)] <- unlist(errors[[2]][i, "CC"]) - unlist(errors[[1]][i, "CC"])
      delta[i,(p-1)] <- unlist(errors[[2]][i, "nCC"]) - unlist(errors[[1]][i, "nCC"])
      delta[i,p] <- unlist(errors[[2]][i, "CCmarg"]) - unlist(errors[[1]][i, "CCmarg"])
    }
  }
  delta[,1] <- errors[[1]][,1]
  colnames(delta) <- colnames(errors[[1]])
  delta
}

##### avlay: wrapper for 'averageLayout'
avlay <- function(..., which.net = 'temporal', threshold = FALSE, 
                  collapse = FALSE, net = NULL, l = NULL, args = NULL){
  if(is.null(l)){
    x <- list(...)
    stopifnot(length(x) > 0)
    if(length(x) == 1){x <- x[[1]]} else if(collapse){x <- appd(x)}
    args0 <- list(object = NA, which.net = which.net, 
                  threshold = threshold, plot = FALSE)
    if(!is.null(net)){args0$which.net <- net}
    args <- append(args[setdiff(names(args), names(args0))], args0)
    averageLayout(lapply(x, function(z) do.call(
      plotNet, replace(args, 'object', list(object = z)))))
  } else if(identical(l, 1) | identical(l, 2) | identical(l, 3) | identical(l, 4)){
      layout(switch(l, t(1:2), t(matrix(1:4, 2, 2)),
                    matrix(c(1, 1, 2, 2, 4, 3, 3, 4), nrow = 2, ncol = 4, byrow = T),
                    matrix(c(4, 1, 1, 4, 2, 2, 3, 3), nrow = 2, ncol = 4, byrow = T)))
  }
}


### ======================================================================== ###
### ======================================================================== ###
##### SURsampler: Sample data for fitting system of lagged SURs
SURsampler <- function(B = NULL, S, n, seed = NULL, beta, beta2 = NULL,
                       full = TRUE, mu = 0, cholesky = FALSE, 
                       time = FALSE, allDat = TRUE){
  p <- ncol(S)
  if(!is.null(B)){
    if(class(B) == "list"){B <- do.call(rbind, B)}
    if(missing(beta)){
      if(ncol(B) == (p + 1)){
        beta <- B
      } else {
        beta <- B[, 1:(p + 1)]
        beta2 <- B[, -c(1:(p + 1)), drop = FALSE]
      }
    }
  }
  if(class(beta) == "list"){beta <- do.call(rbind, beta)}
  if(dim(beta)[2] != dim(beta)[1] + 1){
    if(dim(beta)[2] == dim(beta)[1]){
      beta <- cbind(0, beta)
    } else {stop("Invalid dimensions for beta matrix")}
  }
  if(!is.null(beta2)){
    if(class(beta2) == "list"){beta2 <- do.call(rbind, beta2)}
    stopifnot(dim(beta2)[2] == dim(beta2)[1] - 1)
  }
  if(length(mu) == 1){mu <- rep(mu, p)}
  if(!is.null(seed)){set.seed(seed)}
  tx <- Sys.time()
  if(!cholesky){
    R <- MASS::mvrnorm(n = n, mu = mu, Sigma = S)
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

##### SURfit: fit SUR model with or without constraints
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

##### SURll: log-likelihood of SUR model with LRT to compare models ***FITOBJ***
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


### ======================================================================== ###
### ======================================================================== ###
##### detrender: detrend variables based on linear model
detrender <- function(data, timevar = NULL, vars = NULL, 
                      rmTimevar = TRUE, verbose = TRUE){
  if(!is(data, 'data.frame')){data <- data.frame(data)}
  data_detrend <- data
  if(is.null(timevar)){
    timevar <- 'time'
    data_detrend$time <- 1:nrow(data_detrend)
  }
  if(is.null(vars)){vars <- setdiff(colnames(data_detrend), timevar)}
  vars <- setdiff(vars, timevar)
  for(i in seq_along(vars)){
    ff <- as.formula(paste0(vars[[i]], ' ~ ', timevar))
    fit <- lm(ff, data = data_detrend)
    if(anova(fit)$P[1] < .05){
      if(verbose){message(paste0('Detrending variable ', i))}
      data_detrend[[vars[i]]][!is.na(data_detrend[[vars[i]]])] <- residuals(fit)
    }
  }
  if(rmTimevar){data_detrend <- data_detrend[, setdiff(colnames(data_detrend), timevar)]}
  return(data_detrend)
}

##### getConsec
getConsec <- function(data, beepno = NULL, dayno = NULL, makeAtt = TRUE){
  stopifnot(!is.null(beepno) & !is.null(dayno))
  stopifnot(is.data.frame(data))
  beepday <- list(beepno, dayno)
  stopifnot(sum(sapply(c(1, nrow(data)), function(i){
    all(sapply(beepday, length) == i)})) == 1)
  if(all(sapply(beepday, length) == 1)){
    if(is.character(beepno)){beepno <- which(colnames(data) == beepno)}
    if(is.character(dayno)){dayno <- which(colnames(data) == dayno)}
    data0 <- data[, -c(beepno, dayno)]
    beepday <- as.list(data[, c(beepno, dayno)])
    data <- data0
  }
  consec <- mgm:::beepday2consec(beepvar = beepday[[1]], dayvar = beepday[[2]])
  out <- mgm:::lagData(data = data, lags = 1, consec = consec)[[3]][-1]
  if(makeAtt){
    attr(data, 'samp_ind') <- which(out)
    out <- data
  }
  return(out)
}

##### lagMat: create lagged matrices for fitting SUR models
lagMat <- function(data, type = "g", m = NULL, covariates = NULL, center = TRUE,
                   scale = FALSE, exogenous = TRUE, lags = 1, 
                   consec = NULL, checkType = FALSE){
  if(lags != 1){stop("Only lag = 1 currently supported")}
  if("samp_ind" %in% names(attributes(data))){samp_ind <- attr(data, "samp_ind")}
  data <- data.frame(data)
  vs <- ynames <- colnames(data)
  ynames <- paste0(ynames, ".y")
  binary <- c()
  if(class(type) == "list" | ifelse(length(type) != ncol(data), TRUE, FALSE)){
    if(is.list(type) & "moderators" %in% names(attributes(type))){
      vm <- which(vs %in% attr(type, "moderators"))
    }
    type <- rep("gaussian", ncol(data))
    binary <- unname(which(apply(data, 2, function(z) dim(table(z)) <= 2)))
    if(length(binary) > 0){type[binary] <- "binomial"}
    if(all(type == "binary")){stop("Can't have binary outcomes for lagged models")}
  } else if(any(!type %in% c("g", "gaussian"))){
    binary <- which(!type %in% c("g", "gaussian"))
  }
  names(type) <- ynames
  if(center | scale){
    if(length(binary) > 0){
      data[, -binary] <- apply(data[, -binary], 2, scale, center, scale)
    } else {
      data <- apply(data, 2, scale, center, scale)
    }
  }
  X <- data[-nrow(data), ]
  Y <- data[-1, ]
  colnames(Y) <- ynames
  makeMods <- function(X, m){
    vs <- colnames(X)
    mods <- list()
    if(length(vs) == 2){m <- m[1]}
    for(i in seq_along(m)){
      mform <- as.formula(paste0("~ . * ", vs[m[i]]))
      mods[[i]] <- model.matrix(mform, data.frame(X))
      mods[[i]] <- mods[[i]][, grepl(":", colnames(mods[[i]]))]
    }
    mods <- do.call(cbind, mods)
    if(any(duplicated(colnames(mods)))){mods <- mods[, !duplicated(colnames(mods))]}
    if(ncol(mods) == 1){colnames(mods) <- paste0(colnames(X), collapse = ":")}
    return(mods)
  }
  if(!is.null(covariates)){
    stopifnot(class(covariates) %in% c("numeric", "integer"))
    Y <- as.matrix(Y[, -covariates])
    if(!is.null(m)){
      if(exogenous & length(union(m, covariates)) < length(ynames)){
        ynames <- ynames[-union(m, covariates)]
      } else {
        ynames <- ynames[-covariates]
      }
      covx <- X[, covariates]
      covnames <- vs[covariates]
      X <- X[, -covariates]
      m <- which(colnames(X) %in% vs[m])
      mods <- makeMods(X, m)
      mnames <- colnames(mods)
      xnames <- colnames(X)
      X <- cbind(X, covx, mods)
      colnames(X) <- c(xnames, covnames, mnames)
      if(exogenous & length(m) < ncol(Y)){Y <- Y[, -m]}
    }
  } else if(!is.null(m)){
    if(exogenous & length(m) < ncol(Y)){
      Y <- Y[, -m]
      ynames <- ynames[-m]
    }
    mods <- makeMods(X, m)
    X <- as.matrix(cbind(X, mods))
  }
  if(class(Y) != "matrix"){
    Y <- as.matrix(Y, ncol = 1)
    colnames(Y) <- ynames
  }
  if(is.null(m) & exists("vm", inherits = FALSE) & exogenous){Y <- Y[, -vm]}
  if(exists("samp_ind", inherits = FALSE)){
    Y <- Y[samp_ind, ]
    X <- X[samp_ind, ]
  }
  rownames(Y) <- rownames(X) <- NULL
  #full <- cbind.data.frame(Y, X)
  out <- list(Y = Y, X = X) #, full = full) # FULLFIX
  if("binomial" %in% type[match(colnames(out$Y), names(type))] & checkType){
    stop("Can't have binary outcomes for lagged models")
  }
  if(!is.null(consec)){
    if(exists('samp_ind', inherits = FALSE)){consec <- consec[samp_ind]}
    out$Y <- out$Y[consec, ]
    out$X <- out$X[consec, ]
  }
  return(out)
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

##### net: get adjacency matrices from fit object(s)
net <- function(fit, n = "beta", threshold = FALSE, rule = "OR", 
                binary = FALSE, nodewise = FALSE, d = 14, r = NULL){
  if(inherits(fit, c('splitNets', 'try-error'))){return(NULL)}
  if(is(fit, 'bootNet') | (isTRUE(attr(fit, 'resample')) & 'fit0' %in% names(fit))){fit <- fit$fit0}
  if(is(fit, 'mgmSim')){return(fit[[grep('^trueNet$|^b1$|^adjMat$', names(fit))]])}
  if(is(fit, "matrix")){return(fit)}
  if(is(fit, "mgm")){
    out <- fit$pairwise$wadj * replace(fit$pairwise$signs, is.na(fit$pairwise$signs), 0)
    if(!is.null(r)){out <- out[-r, -r]}
    return(out)
  }
  if(isTRUE(n)){
    if(tolower(threshold) %in% c('and', 'or')){rule <- threshold}
    n <- threshold <- "beta"
  }
  n <- match.arg(tolower(n), c(
    "beta", "contemporaneous", "between", "pdc", "pcc", "kappa", "temporal"))
  n1 <- switch(n, between = "", beta = , pdc = , temporal = "$temporal", "$contemporaneous")
  n2 <- switch(n, kappa = "kappa", pdc = "PDC$adjMat", "adjMat")
  n3 <- paste0("fit", n1, "$", n2)
  if(isTRUE(attr(fit, "mlGVAR"))){fit <- if(n == "between"){fit$betweenNet} else {fit$fixedNets}}
  if("SURnet" %in% names(fit)){fit <- fit$SURnet}
  if(isTRUE(attr(fit, "ggm"))){n3 <- ifelse(nodewise, "fit$nodewise$adjNW", "fit$adjMat")}
  if("fixedPCC" %in% names(fit)){
    threshold <- FALSE
    n <- switch(n, contemporaneous = "pcc", temporal = "beta", n)
    if(n != "between" & "fixedResults" %in% names(fit)){fit <- fit$fixedResults}
    out <- fit[[grep(n, tolower(names(fit)))[1]]]
    if(n == "beta" & ncol(out) != nrow(out)){out <- out[, -1]}
    if(n == "between"){diag(out) <- 0}
    if(n == "pdc"){out <- t(out)}
    if(!is.null(d)){out <- round(out, d)}
  } else {
    if(isTRUE(attr(fit, "lmerVAR")) & n == "between"){n3 <- paste0("fit$between$adjMat")}
    out <- eval(parse(text = n3))
  }
  if(threshold != FALSE){
    if(!is.numeric(threshold)){threshold <- .05}
    rule <- match.arg(tolower(rule), c("or", "and"))
    atts <- names(attributes(fit))
    if("SURnet" %in% atts){
      stopifnot(!n == "between")
      n4 <- ifelse(n %in% c("beta", "temporal", "pdc"), "temporal", "contemporaneous")
      n5 <- ifelse(n4 == "temporal", "fit$temporal$coefs", "fit$contemporaneous")
      pvals <- eval(parse(text = n5))$pvals
      if(ncol(pvals) != nrow(pvals)){
        if(any(grepl("lag", colnames(pvals)))){
          colnames(pvals) <- gsub("[.]lag1[.]$", "", colnames(pvals))
        }
        pvals <- pvals[, intersect(colnames(pvals), gsub("[.]y$", "", rownames(pvals)))]
      }
    } else {
      n4 <- ifelse(n %in% c("beta", "temporal", "pdc"), "Beta", 
                   ifelse(n == "between", "gammaOmega", "gammaTheta"))
      pvals <- if("ggm" %in% atts){fit$nodewise$pvalsNW} else {fit$coefs[[n4]]$Pvals}
      if("ggm" %in% atts){n4 <- ifelse(nodewise & rule == "or", "Beta", "")}
    }
    if(any(is.na(pvals))){pvals[is.na(pvals)] <- 1}
    if(n4 %in% c("temporal", "Beta")){
      out <- out * ifelse(pvals <= threshold, 1, 0)
    } else if(rule == "or"){
      out <- out * ifelse(pvals <= threshold | t(pvals) <= threshold, 1, 0)
    } else if(rule == "and"){
      out <- out * ifelse(pvals <= threshold & t(pvals) <= threshold, 1, 0)
    }
  }
  if(binary){out <- abs(sign(out))}
  if(!is.null(r)){out <- out[-r, -r]}
  return(out)
}

##### netInts: get interactions from fit objects
netInts <- function(fit, n = 'temporal', threshold = FALSE, avg = FALSE, 
                    rule = 'none', r = NULL, empty = TRUE, mselect = NULL){
  rules <- c('none', 'or', 'and')
  eout <- function(fit, empty = TRUE){
    n <- tryCatch({ncol(net(fit))}, error = function(e){TRUE})
    cond <- isTRUE(empty | isTRUE(n) | is.null(n) | is.na(n))
    return(switch(2 - cond, list(), diag(0, n)))
  }
  if(is(fit, 'splitNets')){return(fit$ints)}
  if(is(fit, 'mgm')){
    m <- ifelse(!is.null(r), r, fit$call$moderators[length(fit$call$moderators)])
    mgmInt <- function(x, m){
      ni <- out <- x$interactions$indicator
      if(length(ni) == 2){
        ni <- t(apply(ni[[2]], 1, function(z) setdiff(z, m)))
        vals <- unlist(x$interactions$weightsAgg[[2]]) * x$interactions$signs[[2]]
        out <- matrix(0, ncol(net(x)) - 1, ncol(net(x)) - 1)
        for(i in seq_along(vals)){out[ni[i, 1], ni[i, 2]] <- vals[i]}
        out <- out + t(out)
      }
      return(out)
    }
    out <- mgmInt(fit, m)
    if(!is.null(r) & !is(out, 'list')){out <- out[-r, -r]}
    if(is(out, 'list')){return(eout(fit, empty))} else {return(out)}
  }
  if(is(fit, 'mgmSim')){if('b2' %in% names(fit)){return(fit$b2)} else {return(eout(fit, empty))}}
  if(is(fit, 'bootNet') | (isTRUE(attr(fit, 'resample')) & 'fit0' %in% names(fit))){fit <- fit$fit0}
  if(isTRUE(n) & isTRUE(threshold)){avg <- TRUE}
  if(isTRUE(n)){
    if(tolower(as.character(threshold)) %in% rules){rule <- threshold}
    n <- threshold <- "temporal"
  }
  n <- match.arg(tolower(n), c("temporal", "between"))
  atts <- names(attributes(fit))
  if("mlGVAR" %in% atts){fit <- switch(2 - isTRUE(n == "temporal"), fit$fixedNets, fit$betweenNet)}
  if("SURnet" %in% names(fit)){fit <- fit$SURnet}
  if("lmerVAR" %in% atts){
    stopifnot("ints" %in% names(fit$coefs))
    out <- fit$coefs$ints$coefs
    pvals <- fit$coefs$ints$Pvals
  } else if(any(c('mlVARsim', 'simMLgvar') %in% atts)){
    stopifnot("mm" %in% names(fit))
    out <- fit$mm$mb2
  } else if(!'interactions' %in% names(fit)){
    out <- tryCatch({mmat(fit = fit, m = mselect)}, error = function(e){TRUE})
    if(isTRUE(out)){return(eout(fit, empty))}
    pvals <- out$pvals
    out <- out$betas
  } else {
    out <- fit$interactions[[1]]
    if(isTRUE(attr(fit, "ggm"))){
      pvals <- out[[2]][[2]]
      out <- out[[2]][[1]]
    } else {
      rule <- 'none'
      pvals <- fit$interactions$pvals
      if(ncol(out) != nrow(out) & !is.null(mselect)){
        if(isTRUE(mselect)){mselect <- fit$call$moderators[1]}
        m1 <- paste0(':', mselect); m2 <- paste0(mselect, ':')
        mm <- paste0(c(m1, m2), collapse = '|')
        my <- paste0(gsub('[.]y$', '', rownames(out)), collapse = '|')
        mb1 <- out[, grep(mm, colnames(out))]
        mb1 <- mb1[, grep(my, colnames(mb1))]
        mp1 <- pvals[, grep(mm, colnames(pvals))]
        mp1 <- mp1[, grep(my, colnames(mp1))]
        out <- mb1
        pvals <- mp1
      }
    }
  }
  if(threshold != FALSE & !"mlVARsim" %in% atts){
    rule <- match.arg(tolower(rule), rules)
    if(isTRUE(attr(fit, 'ggm')) & rule == 'none'){rule <- 'or'}
    if(!is.numeric(threshold)){threshold <- .05}
    if(rule == 'none' | isTRUE(attr(fit, 'SURnet'))){
      out <- out * ifelse(pvals <= threshold, 1, 0)
    } else if(rule == 'or'){
      out <- out * ifelse(pvals <= threshold | t(pvals) <= threshold, 1, 0)
    } else if(rule == 'and'){
      out <- out * ifelse(pvals <= threshold & t(pvals) <= threshold, 1, 0)
    }
  }
  if(avg){out <- (t(out) + out)/2}
  if(!is.null(r)){out <- out[-r, -r]}
  if(is(out, 'list')){return(eout(fit, empty))}
  return(out)
}

##### mmat: get interaction matrices when there are multiple exogenous moderators
mmat <- function(fit, m = NULL){
  if(isTRUE(attr(fit, 'ggm'))){
    stopifnot(!is.null(fit$call$moderators))
    allmods <- fit$call$moderators
    p <- length(fit$mods)
    vs <- varnames <- names(fit$mods)
    exind <- function(fit, mm, ex){
      lapply(lapply(lapply(fit$mods, '[[', ifelse(ex == 'b', 'model', 'pvals')), function(i){
        i[grep(mm, rownames(i)), ]
      }), function(z){
        orignames <- names(z)
        names(z) <- gsub(mm, '', names(z))
        if(any(!names(z) %in% vs)){
          orignames <- orignames[which(names(z) %in% vs)]
          z <- z[which(names(z) %in% vs)]
        }
        attr(z, 'orignames') <- orignames
        return(z)
      })
    }
    if(is.null(m)){
      m <- allmods[1]
    } else if(!tolower(m) %in% tolower(allmods)){
      if(any(grepl(tolower(m), tolower(allmods)))){
        m <- allmods[grep(tolower(m), tolower(allmods))[1]]
      } else {
        stop(paste0('No interactions with ', m))
      }
    } else if(!m %in% allmods){
      m <- allmods[which(tolower(allmods) == tolower(m))]
    }
    m1 <- paste0(':', m); m2 <- paste0(m, ':')
    mm <- paste0(c(m1, m2), collapse = '|')
    betas <- exind(fit = fit, mm = mm, ex = 'b')
    pvals <- exind(fit = fit, mm = mm, ex = 'p')
    b1 <- p1 <- structure(matrix(0, p, p), dimnames = rep(list(vs), 2))
    for(i in 1:p){
      nvi <- attr(betas[[i]], 'orignames')
      b1[i, match(names(betas[[i]]), vs)] <- betas[[i]]
      p1[i, match(names(pvals[[i]]), vs)] <- pvals[[i]]
      if(any(grepl(m2, nvi))){
        varnames[which(paste0(m2, vs) %in% nvi)] <- paste0(m2, vs)[which(paste0(m2, vs) %in% nvi)]
      }
    }
    varnames[!grepl(m2, varnames)] <- paste0(varnames[!grepl(m2, varnames)], m1)
    b1 <- t(b1)
    p1 <- t(p1); diag(p1) <- 1
    rownames(b1) <- rownames(p1) <- varnames
    out <- list(betas = b1, pvals = p1)
    attr(out, 'moderator') <- m
  } else {
    stop('Not developed yet')
  }
  return(out)
}

##### matrixDist: compute similarity between two matrices
matrixDist <- function(mat1, mat2 = NULL, ind = "correlation", directed = TRUE, 
                       similarity = TRUE, distMat = FALSE){
  if(length(mat1) == 0){return(NA)}
  if(is.null(mat2) & is(mat1, 'list')){mat2 <- mat1[[2]]; mat1 <- mat1[[1]]}
  mat1 <- as.matrix(mat1); mat2 <- as.matrix(mat2)
  if(!all(dim(mat1) == dim(mat2))){stop("Matrices must have the same dimensions")}
  d <- as.vector(mat1) - as.vector(mat2)
  if(distMat != FALSE){
    distMat <- match.arg(distMat, choices = c(0, 1, 2))
    if(distMat == 0){return(matrix(d, ncol = ncol(mat1), nrow = nrow(mat1)))}
    if(distMat == 1){return(matrix(abs(d), ncol = ncol(mat1), nrow = nrow(mat1)))}
    if(distMat == 2){return(matrix((d^2), ncol = ncol(mat1), nrow = nrow(mat1)))}
  }
  ind <- match.arg(tolower(ind), c("cosine", "mse", "rmse", "ssd", "mae", "msd", "correlation"))
  if(!directed){
    #if(ind == "cosine"){ind <- "correlation"}
    mat1 <- mat1[lower.tri(mat1)]
    mat2 <- mat2[lower.tri(mat2)]
  } else {
    mat1 <- c(mat1)
    mat2 <- c(mat2)
  }
  if(ind == "cosine"){
    #f1 <- norm(x = mat1, type = "F")
    #f2 <- norm(x = mat2, type = "F")
    #top <- sum(diag(t(mat1) %*% mat2))
    f1 <- sqrt(sum(mat1^2))
    f2 <- sqrt(sum(mat2^2))
    top <- sum(mat1 * mat2)
    return(ifelse(similarity == TRUE, top/(f1 * f2), 1 - (top/(f1 * f2))))
  }
  d <- c(mat1) - c(mat2)
  if(ind == "correlation"){
    if(all(is.na(mat1)) || all(is.na(mat2))){return(NA)}
    if(sd(mat1, na.rm = TRUE) == 0 | sd(mat2, na.rm = TRUE) == 0){return(0)}
    return(cor(mat1, mat2, use = "pairwise"))
  }
  if(ind == "ssd"){return(sum(d^2))}
  if(ind == "mse"){return(sum(d^2)/length(d))}
  if(ind == "rmse"){return(sqrt(sum(d^2)/length(d)))}
  if(ind == "mae"){return((sum(abs(d))/length(d)))}
  if(ind == "msd"){return(sum(d/length(d)))}
}


### ======================================================================== ###
### ======================================================================== ###
##### lassoSelect: performs variable selection using the LASSO (glmnet)
lassoSelect <- function(data, yvar, type = "g", criterion = "EBIC", 
                        gamma = .5, nfolds = 10, nlam = 100, alpha = 1){
  if(any(grepl(":", colnames(data)))){data <- data[,!grepl(":", colnames(data))]}
  criterion <- match.arg(criterion, c("CV", "EBIC", "BIC", "AIC"))
  y <- as.numeric(data[, 1])
  x <- data <- as.matrix(data[, -1])
  if(length(type) > 1){type <- type[yvar]}
  fam <- ifelse(type %in% c("g", "gaussian"), "gaussian", "binomial")
  if(criterion == "CV"){
    fit <- cv.glmnet(x = x, y = y, family = fam, type.measure = "deviance", 
                     nfolds = nfolds, nlambda = nlam, alpha = alpha)
  } else {
    fit <- glmnet(x = x, y = y, family = fam, alpha = alpha, nlambda = nlam)
  }
  getIndices <- function(fit, y, x, fam = "gaussian", criterion = "EBIC", 
                         gamma = .5, lam = "null", keepFit = FALSE){
    n <- length(y); p <- ncol(x)
    fam <- match.arg(fam, c("gaussian", "binomial", "multinomial"))
    lam <- match.arg(lam, c("null", "lambda.min", "lambda.1se"))
    if(criterion != "CV"){n_lambdas <- length(fit$lambda)}
    if(fam == "gaussian"){
      if(criterion != "CV"){
        beta0 <- matrix(coef(fit, s = 1)[1], ncol = 1)
        yhat <- rep(1, n) * as.vector(beta0)
        n_neighbors <- sapply(1:n_lambdas, function(z){
          colSums(as.matrix(coef(fit, s = fit$lambda[z])[-1,]) != 0)
        })
        LL_model <- sum(dnorm(y, mean = yhat, sd = sqrt(sum((y - yhat)^2)/n), log = TRUE))
      } else {
        mods <- lapply(c("lambda.min", "lambda.1se"), function(z){
          betas <- matrix(coef(fit, s = z), ncol = 1)
          yhat <- cbind(1, x) %*% as.vector(betas)
          n_neighbors <- colSums(matrix(coef(fit, s = z)[-1, ], ncol = 1) != 0)
          LL_model <- sum(dnorm(y, mean = yhat, sd = sqrt(sum((y - yhat)^2)/n), log = TRUE))
          ic_lambda <- -2 * LL_model + n_neighbors * log(n) + 2 * gamma * n_neighbors * log(p)
          return(list(betas = matrix(betas[-1, ], ncol = 1), EBIC = ic_lambda))
        })
      }
    } else if(fam %in% c("multinomial", "binomial")){
      lam <- ifelse(criterion == "CV", list(c("lambda.min", "lambda.1se")), list(1))[[1]]
      mods <- lapply(lam, function(lam0){
        cats <- unique(y)
        n_cats <- length(cats)
        m_respdum <- matrix(NA, n, n_cats)
        m_coefs <- matrix(NA, n, n_cats)
        m_LL_parts <- matrix(NA, nrow = n, ncol = n_cats + 1)
        X <- cbind(rep(1, n), x)
        for(catIter in 1:n_cats){
          m_respdum[, catIter] <- (y == cats[catIter]) * 1
          if(fam == "multinomial"){
            m_coefs[, catIter] <- X %*% matrix(coef(fit, s = lam0)[[catIter]], ncol = 1)
          } else {
            m_coefs[, catIter] <- X %*% ((matrix(coef(fit, s = lam0), ncol = 1) * ifelse(catIter == 1, -1, 1))/2)
          }
          m_LL_parts[, catIter] <- m_respdum[, catIter] * m_coefs[, catIter]
        }
        m_LL_parts[, (n_cats + 1)] <- -log(rowSums(exp(m_coefs)))
        LL_model <- sum(rowSums(m_LL_parts))
        if(lam0 == 1){
          n_lambdas <- length(fit$lambda)
          n_neighbors <- c()
          for(NN in 1:n_lambdas){
            coefs_bin <- vector("list", length = n_cats)
            for(ca in 1:n_cats){
              if(fam == "multinomial"){
                coefs_bin[[ca]] <- as.matrix(coef(fit, s = fit$lambda[NN])[[ca]][-1, ]) != 0
              } else {
                coefs_bin[[ca]] <- as.matrix(coef(fit, s = fit$lambda[NN])) != 0
              }
            }
            n_neighbors[NN] <- colSums(Reduce("+", coefs_bin) != 0)
            if(fam == "binomial"){n_neighbors[NN] <- n_neighbors[NN] - 1}
          }
          return(list(LL_model = LL_model, n_neighbors = n_neighbors))
        } else {
          coefs_bin <- vector("list", length = n_cats)
          if(fam == "multinomial"){
            for(ca in 1:n_cats){coefs_bin[[ca]] <- as.matrix(coef(fit, s = lam0)[[ca]][-1, ]) != 0}
          } else {
            for(ca in 1:n_cats){coefs_bin[[ca]] <- as.matrix(coef(fit, s = lam0)[-1, ]) != 0}
          }
          n_neighbors <- colSums(Reduce("+", coefs_bin) != 0)
          ic_lambda <- -2 * LL_model + n_neighbors * log(n) + 2 * gamma * n_neighbors * log(p)
          betas <- matrix(coef(fit, s = lam0), ncol = 1)
          return(list(betas = matrix(betas[-1, ], ncol = 1), EBIC = ic_lambda))
        }
      })
      if(criterion != "CV"){
        LL_model <- mods[[1]]$LL_model
        n_neighbors <- mods[[1]]$n_neighbors
      }
    }
    if(criterion != "CV"){
      LL_sat <- 1/2 * fit$nulldev + LL_model
      deviance <- (1 - fit$dev.ratio) * fit$nulldev
      LL_lambda_models <- -1/2 * deviance + LL_sat
      ic_lambda <- -2 * LL_lambda_models + n_neighbors * ifelse(
        criterion == "AIC", 2, log(n)) + ifelse(
          criterion == "EBIC", list(2 * gamma * n_neighbors * log(p)), list(0))[[1]]
      allCoefs <- lapply(seq_len(n_lambdas), function(z) coef(fit)[, z])
      betas <- allCoefs[[which.min(ic_lambda)]][-1]
      coefs <- Matrix(betas, sparse = TRUE)
      rownames(coefs) <- names(betas)
      fitobj <- list(fit = fit, fit0 = NA, crit = ic_lambda)
      if(keepFit){fitobj$fit0 <- glmnet(x, y, fam, lambda = fit$lambda[which.min(ic_lambda)])}
      names(fitobj)[3] <- criterion
      output <- list(mod0 = names(betas)[betas != 0], coefs = coefs, 
                     fitobj = fitobj, allCoefs = allCoefs)
      if(length(output$mod0) == 0){output$mod0 <- 1}
    } else {
      coefs <- Matrix(do.call(cbind, lapply(mods, '[[', "betas")), sparse = TRUE)
      rownames(coefs) <- colnames(x)
      colnames(coefs) <- paste0("mod", c("0", "1se"))
      fitobj <- list(fitCV = fit, fit0 = NA, fit1se = NA)
      if(keepFit){
        fitobj$fit0 <- glmnet(x, y, fam, lambda = fit$lambda.min)
        fitobj$fit1se <- glmnet(x, y, fam, lambda = fit$lambda.1se)
      }
      attr(fitobj$fit0, "EBIC") <- mods[[1]]$EBIC
      attr(fitobj$fit1se, "EBIC") <- mods[[2]]$EBIC
      output <- list(mod0 = colnames(x)[coefs[, 1] != 0], 
                     mod1se = colnames(x)[coefs[, 2] != 0],
                     coefs = coefs, fitobj = fitobj)
      if(length(output$mod0) == 0){output$mod0 <- 1}
      if(length(output$mod1se) == 0){output$mod1se <- 1}
    }
    attr(output, "family") <- fam
    return(output)
  }
  out <- getIndices(fit, y, x, fam, criterion, gamma)
  return(out)
}

##### fitHierLASSO: performs variable selection using the hierarchical LASSO
fitHierLASSO <- function(data, yvar, type = "g", m = NULL, criterion = "CV", 
                         method = "glinternet", gamma = .5, nfolds = 10, 
                         nlam = 50, lags = NULL, useSE = TRUE, diag = FALSE, 
                         outMsgs = FALSE, dmnames = NULL, verbose = TRUE){
  if(any(grepl(":", colnames(data)))){data <- data[,!grepl(":", colnames(data))]}
  method <- match.arg(tolower(method), c("hiernet", "glinternet"))
  criterion <- match.arg(criterion, c("CV", "EBIC", "BIC", "AIC"))
  y <- as.numeric(data[, 1])
  x <- data <- as.matrix(data[, -1])
  if(method == "hiernet"){
    out1 <- capture.output({fitPath <- hierNet.path(x, y, nlam = nlam, strong = TRUE, diagonal = diag)})
    out2 <- capture.output({fitCV <- hierNet.cv(fitPath, x, y, nfolds = nfolds)})
    out3 <- capture.output({fit0 <- hierNet(x, y, lam = fitCV$lamhat, strong = TRUE, diagonal = diag)})
    out4 <- capture.output({fit1se <- hierNet(x, y, lam = fitCV$lamhat.1se, strong = TRUE, diagonal = diag)})
    mod0 <- c(fit0$bp - fit0$bn, fit0$th[lower.tri(fit0$th)])
    mod1se <- c(fit1se$bp - fit1se$bn, fit1se$th[lower.tri(fit1se$th)])
    coefs <- Matrix(cbind(mod0, mod1se), sparse = TRUE)
  } else if(method == "glinternet"){
    if(length(type) > 1){type <- type[yvar]}
    fam <- ifelse(type %in% c("g", "gaussian"), "gaussian", "binomial")
    type <- rep(1, ncol(x))
    if(criterion == "CV"){
      fitCV <- tryCatch({glinternet.cv(x, y, type, nFolds = nfolds, nLambda = nlam, 
                                       interactionCandidates = m, family = fam)}, 
                        error = function(e){
                          failed <- TRUE
                          take <- 1
                          if(verbose){cat("\n")}
                          while(failed == TRUE){
                            if(take <= 5){
                              if(verbose){cat("  Failed.. trying again, take =", take, "\n")}
                              fitCV <- try(glinternet.cv(x, y, type, nFolds = nfolds, nLambda = nlam, 
                                                         interactionCandidates = m, family = fam), silent = TRUE)
                              if(class(fitCV) == "try-error"){
                                failed <- TRUE
                                take <- take + 1
                              } else {
                                failed <- FALSE
                              }
                            } else if(take <= 10){
                              if(verbose){cat("  Failed.. trying nlam = 20, take =", take, "\n")}
                              fitCV <- try(glinternet.cv(x, y, type, nFolds = nfolds, nLambda = 20, 
                                                         interactionCandidates = m, family = fam), silent = TRUE)
                              if(class(fitCV) == "try-error"){
                                failed <- TRUE
                                take <- take + 1
                              } else {
                                failed <- FALSE
                              }
                            } else {
                              if(verbose){cat("  Failed.. trying nlam = 20 & nFolds = 3, take =", take, "\n")}
                              fitCV <- try(glinternet.cv(x, y, type, nFolds = 3, nLambda = 20, 
                                                         interactionCandidates = m, family = fam), silent = TRUE)
                              if(class(fitCV) == "try-error"){
                                failed <- TRUE
                                take <- take + 1
                                if(take == 20){break}
                              } else {
                                failed <- FALSE
                              }
                            }
                          }
                          fitCV
                        })
      if(useSE == TRUE){
        lamlist <- fitCV$lambda
        errm <- fitCV$cvErr
        errse <- fitCV$cvErrStd <- fitCV$cvErrStd/sqrt(nfolds)
        o <- which.min(errm)
        lamhat <- lamlist[o]
        oo <- errm <= errm[o] + errse[o]
        fitCV$lambdaHat1Std <- lamlist[oo & lamlist >= lamhat][1]
      }
      which.lam0 <- which(fitCV$lambda == fitCV$lambdaHat)
      while(is.null(fitCV$glinternetFit$activeSet[[which.lam0]])){
        if(verbose){cat("\n  Mod0 empty.. choosing new lambda\n")}
        which.lam0 <- which.lam0 + 1
      }
      fitCV$lambdaHat <- fitCV$lambda[which.lam0]
      which.lam1se <- which(fitCV$lambda == fitCV$lambdaHat1Std)
      while(is.null(fitCV$glinternetFit$activeSet[[which.lam1se]])){
        if(verbose){cat("\n  Mod1SE empty.. choosing new lambda\n")}
        which.lam1se <- which.lam1se + 1
      }
      fitCV$lambdaHat1Std <- fitCV$lambda[which.lam1se]
      fit0 <- glinternet(x, y, type, lambda = fitCV$lambdaHat, 
                         interactionCandidates = m, family = fam)
      fit1se <- glinternet(x, y, type, lambda = fitCV$lambdaHat1Std, 
                           interactionCandidates = m, family = fam)
      attributes(fit1se)$useSE <- attributes(fitCV)$useSE <- useSE
      mod0 <- coef(fit0)[[2]]
      mod1se <- coef(fit1se)[[2]]
    } else {
      fit <- glinternet(x, y, type, interactionCandidates = m, 
                        family = fam, nLambda = nlam)
      coefs <- coef(fit)[-1]
      mains <- 1:ncol(x)
      ints <- t(combn(mains, 2))
      ints2 <- as.numeric(apply(ints, 1, paste, collapse = ""))
      if(is.null(lags) & !is.null(m)){
        vs <- colnames(x)
        vs1 <- vs[m]
        if(length(vs1) > 1){vs1 <- paste0("(", paste(vs1, collapse = " + "), ")")}
        vs2 <- as.formula(paste0("~ . * ", vs1))
        dmnames <- colnames(model.matrix(vs2, data.frame(x)))[-1]
      }
      allCoefs <- lapply(coefs, function(z){
        zmain1 <- z$mainEffects$cont
        zmain2 <- z$mainEffectsCoef$cont
        if(length(zmain1) != 0){
          if(any(!mains %in% zmain1)){
            zmiss1 <- mains[!mains %in% zmain1]
            zcoefs1 <- c(zmain2, rep(0, length(zmiss1)))[order(c(zmain1, zmiss1))]
          } else {
            zcoefs1 <- zmain2[order(zmain1)]
          }
        } else {
          zcoefs1 <- rep(0, length(mains))
        }
        zint1 <- z$interactions$contcont
        zint2 <- z$interactionsCoef$contcont
        if(length(zint1) != 0){
          zints1 <- as.numeric(apply(zint1, 1, paste, collapse = ""))
          if(nrow(ints) != nrow(zint1)){
            zcoefs2 <- rep(0, nrow(ints))
            zcoefs2[which(ints2 %in% zints1)] <- zint2
          } else {
            zcoefs2 <- zint2[match(zints1, ints2)]
          }
        } else {
          zcoefs2 <- rep(0, nrow(ints))
        }
        betas <- unlist(c(zcoefs1, zcoefs2))
        names(betas) <- c(colnames(x), apply(combn(colnames(x), 2), 2, paste, collapse = ":"))
        if(!is.null(m)){
          betas <- betas[which(names(betas) %in% dmnames)]
          #if(is.null(lags)){
          #  x2 <- c(colnames(x), paste0(colnames(x)[-m], ":", colnames(x)[m]))
          #  betas <- betas[which(names(betas) %in% x2)]
          #} else {
          #  betas <- betas[which(names(betas) %in% dmnames)]
          #}
        }
        return(betas)
      })
      n_neighbors <- sapply(allCoefs, function(z) sum(z != 0))
      LL_models <- sapply(1:length(allCoefs), function(z){
        s2 <- sum((y - fit$fitted[, z + 1])^2)/length(y)
        sum(dnorm(y, mean = fit$fitted[, z + 1], sd = sqrt(s2), log = TRUE))
      })
      p <- length(mains) + nrow(ints)
      if(!is.null(m)){
        if(all(m == 0)){
          p <- length(mains)
        } else {
          p <- length(c(colnames(x), paste0(colnames(x)[-m], ":", colnames(x)[m])))
        }
      }
      ic_lambda <- -2 * LL_models + n_neighbors * ifelse(
        criterion == "AIC", 2, log(nrow(x))) + ifelse(
          criterion == "EBIC", list(2 * gamma * n_neighbors * log(p)), list(0))[[1]]
      betas <- allCoefs[[which.min(ic_lambda)]]
      lambda_min <- fit$lambda[which.min(ic_lambda) + 1]
      coefs <- Matrix(betas, sparse = TRUE)
      rownames(coefs) <- names(betas)
      fitobj <- list(fit = fit, fit0 = NA, crit = ic_lambda)
      if(ifelse(!is.null(m), ifelse(all(m == 0), FALSE, TRUE), TRUE) & method == "hiernet"){
        fitobj$fit0  <- glinternet(x, y, type, interactionCandidates = m, 
                                   lambda = lambda_min, family = fam)
      }
      names(fitobj)[3] <- criterion
      output <- list(mod0 = names(betas)[betas != 0], coefs = coefs, 
                     fitobj = fitobj, allCoefs = allCoefs)
      attr(output, "family") <- ifelse(fam == "gaussian", "g", "c")
      return(output)
    }
    mains <- 1:ncol(x)
    ints <- t(combn(mains, 2))
    ints2 <- as.numeric(apply(ints, 1, paste, collapse = ""))
    if(length(mod0$mainEffects$cont) != 0){
      if(any(!mains %in% mod0$mainEffects$cont)){
        mod0miss1 <- mains[!mains %in% mod0$mainEffects$cont]
        mod0coefs1 <- c(mod0$mainEffectsCoef$cont, 
                        rep(0, length(mod0miss1)))[order(c(mod0$mainEffects$cont, mod0miss1))]
      } else {
        mod0coefs1 <- mod0$mainEffectsCoef$cont[order(mod0$mainEffects$cont)]
      }
    } else {
      mod0coefs1 <- rep(0, length(mains))
    }
    if(length(mod1se$mainEffects$cont) != 0){
      if(any(!mains %in% mod1se$mainEffects$cont)){
        mod1semiss1 <- mains[!mains %in% mod1se$mainEffects$cont]
        mod1secoefs1 <- c(mod1se$mainEffectsCoef$cont, 
                          rep(0, length(mod1semiss1)))[order(c(mod1se$mainEffects$cont, mod1semiss1))]
      } else {
        mod1secoefs1 <- mod1se$mainEffectsCoef$cont[order(mod1se$mainEffects$cont)]
      }
    } else {
      mod1secoefs1 <- rep(0, length(mains))
    }
    if(length(mod0$interactions$contcont) != 0){
      mod0ints1 <- as.numeric(apply(mod0$interactions$contcont, 1, paste, collapse = ""))
      if(nrow(ints) != nrow(mod0$interactions$contcont)){
        mod0coefs2 <- rep(0, nrow(ints))
        mod0coefs2[which(ints2 %in% mod0ints1)] <- mod0$interactionsCoef$contcont
      } else {
        mod0coefs2 <- mod0$interactionsCoef$contcont[match(mod0ints1, ints2)]
      }
    } else {
      mod0coefs2 <- rep(0, nrow(ints))
    }
    if(length(mod1se$interactions$contcont) != 0){
      mod1seints1 <- as.numeric(apply(mod1se$interactions$contcont, 1, paste, collapse = ""))
      if(nrow(ints) != nrow(mod1se$interactions$contcont)){
        mod1secoefs2 <- rep(0, nrow(ints))
        mod1secoefs2[which(ints2 %in% mod1seints1)] <- mod1se$interactionsCoef$contcont
      } else {
        mod1secoefs2 <- mod1se$interactionsCoef$contcont[match(mod1seints1, ints2)]
      }
    } else {
      mod1secoefs2 <- rep(0, nrow(ints))
    }
    mod0 <- unlist(c(mod0coefs1, mod0coefs2))
    mod1se <- unlist(c(mod1secoefs1, mod1secoefs2))
    coefs <- Matrix(cbind(mod0, mod1se), sparse = TRUE)
  }
  allNames <- c(colnames(data), apply(combn(colnames(data), 2), 2, paste, collapse = ":"))
  rownames(coefs) <- allNames
  if(outMsgs == TRUE & method == "hiernet"){
    output <- list(mod0 = allNames[mod0 != 0], mod1se = allNames[mod1se != 0], coefs = coefs, 
                   fitobj = list(fitCV = fitCV, fit0 = fit0, fit1se = fit1se),
                   outMsgs = list(outPath = out1, outCV = out2, out0 = out3, out1se = out4))
  } else {
    output <- list(mod0 = allNames[mod0 != 0], mod1se = allNames[mod1se != 0], coefs = coefs, 
                   fitobj = list(fitCV = fitCV, fit0 = fit0, fit1se = fit1se))
  }
  attr(output, "family") <- ifelse(method == "glinternet", ifelse(fam == "gaussian", "g", "c"), "g")
  return(output)
}

##### varSelect: returns output for a given variable selection procedure
varSelect <- function(data, m = NULL, criterion = "AIC", method = "glmnet",
                      lags = NULL, exogenous = TRUE, type = "g", center = TRUE,
                      scale = FALSE, gamma = .5, nfolds = 10, varSeed = NULL, 
                      useSE = TRUE, nlam = NULL, covs = NULL, verbose = TRUE, 
                      beepno = NULL, dayno = NULL){
  dat <- data
  ALL <- FALSE
  dmnames <- NULL # VARSELECT START
  mall <- which(sapply(c(0, 'all'), identical, as.character(m)))
  if(any(mall)){
    m <- switch(mall, NULL, 1:ncol(dat))
    if(mall == 2){ALL <- TRUE}
  }
  if(is.null(lags)){
    if(length(m) >= ncol(dat) - 1){exogenous <- FALSE}
    if(!exogenous){ALL <- TRUE}
  } else if(any(!sapply(c(beepno, dayno), is.null))){
    stopifnot(!is.null(beepno) & !is.null(dayno))
    dat <- getConsec(data = dat, beepno = beepno, dayno = dayno)
  }
  if(!is(dat, 'list')){
    dat <- structure(data.frame(dat), samp_ind = attr(data, "samp_ind"))
    if(is.null(m) | ALL | !is.null(lags)){
      if(!is.null(lags)){
        vs <- colnames(dat)
        if(!is.null(m) & !ALL){mname <- vs[m]}
        if(!is.null(covs)){
          dat <- dat[, -covs]
          if(!is.null(m) & !ALL){m <- which(colnames(dat) %in% mname)}
        }
        dmnames <- colnames(lagMat(data = dat, m = m)$X)
        dat <- lagMat(data = dat, type = type, center = center, scale = scale)
        dat$full <- cbind.data.frame(dat$X, dat$Y) # FULLFIX
        if(ALL | (!is.null(m) & all(m == 0))){exogenous <- FALSE; m <- 1:ncol(dat$Y)}
        #if(ALL | is.null(m)){exogenous <- FALSE; m <- 1:ncol(dat$Y)}
        if(!is.null(m)){
          if(exogenous & length(m) < ncol(dat$Y)){
            dat$Y <- dat$Y[, -m]
            dat$full <- dat$full[, -m]
            if(class(dat$Y) != "matrix"){
              dat$Y <- as.matrix(dat$Y, ncol = 1)
              colnames(dat$Y) <- colnames(dat$full)[1]
            }
          } else if(length(m) == ncol(dat$Y)){
            mname <- gsub("[.]y$", "", colnames(dat$Y))
            exogenous <- FALSE
            ALL <- TRUE
          }
        }
      } else {
        mname <- colnames(dat)[m]
        dat <- list(Y = dat, X = dat)
      }
    } else if(length(m) >= ncol(dat) - 1){
      mname <- colnames(dat)[m]
      exogenous <- FALSE; ALL <- TRUE
    } else if(class(m) == "list"){
      dat <- list(Y = dat, X = data.frame(dat, m))
      m <- ncol(dat$Y)
    } else if(class(m) %in% c("numeric", "integer")){
      dn <- colnames(dat)
      dat <- list(Y = dat[, -m], X = data.frame(dat[, -m], dat[, m]))
      colnames(dat$X) <- c(dn[-m], dn[m])
      m <- ncol(dat$Y)
    }
  }
  if(!is.null(m)){
    if(class(m) == "list"){
      stopifnot(!ALL)
      mname <- names(m)
      m <- ncol(dat$Y)
      if(m < (ncol(dat$X) - 1)){
        covariates <- TRUE
        dn <- colnames(dat$X)
        m <- which(dn == mname)
        dat$X <- data.frame(dat$X[, -m], dat$X[, m])
        colnames(dat$X) <- c(dn[-m], dn[m])
        m <- ncol(dat$X) - 1
      }
    }
    if(all(m != 0)){m0 <- m}
    if(!method %in% c('hiernet', 'glinternet')){method <- 'glinternet'}
  }
  if(!is.null(varSeed)){set.seed(varSeed)}
  p <- ncol(dat$Y); data <- dat$X; Y <- dat$Y
  method <- match.arg(tolower(method), c(
    "hiernet", "glinternet", "subset", "backward", 
    "forward", "seqrep", "glmnet", "lasso"))
  criterion <- toupper(match.arg(tolower(criterion), c(
    "cv", "aic", "bic", "ebic", "cp", "rss", "adjr2", "rsq", "r2")))
  ### VARSELECT
  if(method %in% c("hiernet", "glinternet")){
    if(is.null(nlam)){nlam <- ifelse(method == "hiernet", 20, 50)}
    hierMods <- list(); t <- c()
    for(i in 1:p){
      if(all(colnames(dat$Y) %in% colnames(dat$X))){
        data <- cbind(y = Y[, i], dat$X[, -i])
      } else {
        data <- cbind(y = Y[, i], dat$X)
      }
      if(verbose == TRUE){cat("Fitting model ", i, "/", p, "...", sep = "")}
      if(verbose == "pbar"){if(i == 1){pb <- txtProgressBar(min = 0, max = p, style = 2, width = 43)}}
      if(ALL & !is.null(m)){
        if(all(m != 0)){
          m <- switch(2 - (i %in% m), NULL, which(colnames(data)[-1] %in% mname))
        }
      }
      tx <- Sys.time()
      hierMods[[i]] <- fitHierLASSO(data = data, yvar = i, type = type, m = m,
                                    nlam = nlam, nfolds = nfolds, gamma = gamma,
                                    method = method, useSE = useSE, lags = lags,
                                    criterion = criterion, dmnames = dmnames, 
                                    verbose = ifelse(is.logical(verbose), verbose, FALSE))
      t[i] <- tx <- Sys.time() - tx
      names(t)[i] <- attr(tx, "units")
      if(ALL & exists("m0", inherits = FALSE)){m <- m0}
      if(verbose == TRUE){
        cat("  Complete! ", "(", round(t[i], 2), " ", attr(tx, "units"), ")\n", sep = "")}
      if(verbose == "pbar"){setTxtProgressBar(pb, i); if(i == p){close(pb)}}
    }
    if(verbose == TRUE){
      if(length(unique(names(t))) == 1){
        cat("####### Total time:", round(sum(t), 2), names(t)[1], "\n\n")
      } else {
        tt <- t
        if(length(unique(names(tt))) == 2){
          if(all(sort(unique(names(tt))) == c("mins", "secs"))){
            tt[names(tt) == "mins"] <- 60 * tt[names(tt) == "mins"]
            cat("####### Total time:", round(sum(tt)/60, 2), "mins\n\n")
          } else if(all(sort(unique(names(tt))) == c("hours", "mins"))){
            tt[names(tt) == "hours"] <- 60 * tt[names(tt) == "hours"]
            cat("####### Total time:", round(sum(tt)/60, 2), "hours\n\n")
          }
        } else if(all(sort(unique(names(tt))) == c("hours", "mins", "secs"))){
          tt[names(tt) == "hours"] <- 360 * tt[names(tt) == "hours"]
          tt[names(tt) == "mins"] <- 60 * tt[names(tt) == "mins"]
          cat("####### Total time:", round(sum(tt)/360, 2), "hours\n\n")
        }
      }
    }
    names(hierMods) <- colnames(Y)
    attributes(hierMods)$method <- method
    attributes(hierMods)$criterion <- criterion
    if(criterion == "EBIC"){attributes(hierMods)$gamma <- gamma}
    attributes(hierMods)$time <- t
    if(exists("covariates", inherits = FALSE)){attributes(hierMods)$covariates <- TRUE}
    if(!is.null(covs)){attributes(hierMods)$covs <- covs}
    if(!is.null(lags)){
      attributes(hierMods)$moderators <- mname
      attributes(hierMods)$exogenous <- exogenous
    }
    return(hierMods)
  }
  ### LASSO SELECTION
  if(method %in% c("lasso", "glmnet")){
    lassoMods <- list(); t <- c()
    if(is.null(nlam)){nlam <- 100}
    for(i in 1:p){
      if(all(colnames(dat$Y) %in% colnames(dat$X))){
        data <- cbind(y = Y[, i], dat$X[, -i])
      } else {
        data <- cbind(y = Y[, i], dat$X)
      }
      if(verbose != FALSE){if(i == 1){pb <- txtProgressBar(min = 0, max = p, style = 2, width = 43)}}
      tx <- Sys.time()
      lassoMods[[i]] <- lassoSelect(data = data, yvar = i, criterion = criterion, 
                                    type = type, gamma = gamma, nfolds = nfolds, 
                                    nlam = nlam)
      t[i] <- tx <- Sys.time() - tx
      names(t)[i] <- attr(tx, "units")
      if(verbose != FALSE){setTxtProgressBar(pb, i); if(i == p){close(pb)}}
    }
    names(lassoMods) <- colnames(Y)
    attributes(lassoMods)[c("method", "criterion", "time")] <- list("glmnet", criterion, t)
    if(criterion %in% c("EBIC", "CV")){attr(lassoMods, "gamma") <- gamma}
    if(exists("covariates", inherits = FALSE)){attributes(lassoMods)$covariates <- TRUE}
    if(!is.null(covs)){attributes(lassoMods)$covs <- covs}
    return(lassoMods)
  }
  if(method %in% c("subset", "backward", "forward", "seqrep")){
    if(criterion == "CV"){criterion <- "Cp"}
    if(tolower(criterion) %in% c("aic", "ebic")){criterion <- "bic"}
    ind <- match.arg(tolower(criterion), c("cp", "bic", "adjr2", "rsq", "r2", "rss"))
    if(method == "subset"){method <- "exhaustive"}
    if(ind == "r2"){ind <- "rsq"}
    best <- ifelse(ind %in% c("cp", "bic", "rss"), which.min, which.max)
    regMods <- bestMods <- list()
    for(i in 1:p){
      if(all(colnames(dat$Y) %in% colnames(dat$X))){data <- dat$X[,-i]}
      regMods[[i]] <- summary(regsubsets(data, Y[,i], nvmax = ncol(data), method = method))
      bestMods[[i]] <- ifelse(regMods[[i]]$which, 1, 0)[best(regMods[[i]][[ind]]), -1]
    }
    bestMods <- lapply(bestMods, function(z) names(z)[z == 1])
    bestMods <- lapply(1:p, function(z){
      bestOut <- list(mod0 = bestMods[[z]], fitobj = regMods[[z]])
      attr(bestOut, "family") <- "g"
      return(bestOut)
    })
    names(bestMods) <- colnames(Y)
    attributes(bestMods)$method <- "regsubsets"
    attributes(bestMods)$criterion <- ind
    return(bestMods)
  }
}

##### resample: bootstrapping or multi-sample splits for variable selection
resample <- function(data, m = NULL, niter = 10, sampMethod = "bootstrap", criterion = "AIC",
                     method = "glmnet", rule = "OR", gamma = .5, nfolds = 10, 
                     nlam = 50, which.lam = "min", threshold = FALSE, bonf = FALSE, 
                     alpha = .05, exogenous = TRUE, split = .5, center = TRUE, 
                     scale = FALSE, varSeed = NULL, seed = NULL, verbose = TRUE, 
                     lags = NULL, binary = NULL, type = 'g', saveMods = TRUE,
                     saveData = FALSE, saveVars = FALSE, fitit = TRUE, 
                     nCores = 1, cluster = 'mclapply', block = FALSE,
                     beepno = NULL, dayno = NULL, ...){
  # The 'split' argument, for sampMethod %in% c('split', 'stability') 
  # Is the proportion of the data to be included in the training set
  t1 <- Sys.time() # RESAMPLE START
  if(is.null(lags) & !identical(block, FALSE)){
    message('Assuming lags = 1 since block resampling was requested')
    lags <- 1
  }
  if(!is.null(lags)){
    lags <- switch(2 - identical(as.numeric(lags), 0), NULL, 1)
    if(any(!sapply(c(beepno, dayno), is.null))){
      stopifnot(!is.null(beepno) & !is.null(dayno))
      data <- getConsec(data = data, beepno = beepno, dayno = dayno)
    }
  }
  consec <- switch(2 - (!is.null(lags) & 'samp_ind' %in% names(attributes(data))), 
                   attr(data, 'samp_ind'), NULL)
  N <- ifelse(is.null(consec), nrow(data) - ifelse(is.null(lags), 0, lags), length(consec))
  data <- data.frame(data)
  if(!is.null(m)){if(all(m == 0)){m <- NULL}}
  method <- ifelse(!is.null(m), 'glinternet', ifelse(
    !method %in% c('glmnet', 'subset'), 'glmnet', method))
  if(isTRUE(block)){block <- floor(3.15 * N^(1/3))}
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
                 varSeed = varSeed, type = type, lags = lags, block = block)
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
  sampler2 <- function(N, size = N, block = FALSE, method = 'bootstrap', consec = NULL){
    allinds <- switch(2 - is.null(consec), seq_len(N), consec)
    if(identical(block, FALSE)){
      sample(allinds, size, replace = (method == 'bootstrap'))
    } else if(method == 'bootstrap'){
      stopifnot(block < N & block > 1)
      nblox <- 1
      while(N > block * nblox){nblox <- nblox + 1}
      possible <- seq_len(N - block)
      starts <- replicate(nblox, sample(possible, 1))
      sampinds <- c(sapply(starts, function(z) allinds[z:(z + block - 1)]))
      if(length(sampinds) > N){sampinds <- sampinds[1:N]}
      return(sampinds)
    } else {
      stop('Multi-split block resampling not yet implemented')
    }
  }
  sampInd <- samps <- vars <- vars1 <- fits <- train <- test <- list()
  if((exogenous | is.null(m)) & sampMethod != 'bootstrap'){
    mno <- as.numeric(!is.null(m) & !identical(m, 0))
    lno <- as.numeric(!is.null(lags) & !identical(lags, 0))
    minsize <- sampleSize(p = ncol(data) - mno, m = mno,
                          lags = lno, print = FALSE)
    minn <- ifelse(split > .5, N - floor(N * split), floor(N * split))
    if(minn < minsize & split != .5){
      split <- .5
      minn <- ifelse(split > .5, N - floor(N * split), floor(N * split))
      if(minn > minsize){warning('Split size set to .5 to produce large enough subsamples')}
    }
    if(minn < minsize){stop('Sample-split size is smaller than required')}
  }
  if(nCores > 1 | isTRUE(nCores)){
    if(isTRUE(nCores)){nCores <- parallel::detectCores()}
    if(grepl('Windows', sessionInfo()$running)){cluster <- 'SOCK'}
    if(sampMethod != 'bootstrap'){
      if(split <= 0 | split >= 1){stop('split size must be between 0 and 1')}
      n <- floor(N * split)
    } else {
      n <- N
    }
    for(i in 1:niter){ #### SAMPLING 1
      set.seed(seeds[i])
      sampInd[[i]] <- sampler2(N = N, size = n, block = block, method = sampMethod, consec = consec)
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
    } ##### SAMPLING 1 END
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
      for(i in 1:niter){ ##### SAMPLING 2
        set.seed(seeds[i]) # RESAMPLE SPLIT
        sampInd[[i]] <- sampler2(N = N, size = n, block = block, method = sampMethod, consec = consec)
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
        if(verbose){ ##### SAMPLING 2 END
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
      for(i in 1:niter){ ##### SAMPLING 3
        set.seed(seeds[i]) # RESAMPLE BOOT
        sampInd[[i]] <- sampler2(N = N, block = block, consec = consec) # Currently only for sampMethod == 'bootstrap'
        if(is.null(lags)){
          samps[[i]] <- data[sampInd[[i]], ]
        } else {
          samps[[i]] <- data
          attr(samps[[i]], "samp_ind") <- sampInd[[i]]
        }
        if(verbose == TRUE){ ##### SAMPLING 3 END
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
  ### RESAMPLE VARMODS
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
    attr(out, 'resample') <- TRUE
    class(out) <- c('list', 'resample')
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
  sci <- sapply(adjCIs0, function(z) !all(z$select_ci == z$select))
  if(any(sci)){
    message(paste0('CIs different than p-values for:\n', paste(
      names(sci)[sci], collapse = '\n')))
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
  class(out) <- c('list', 'resample')
  out
}

##### bootCover: coverage probabilities based on bootstrapped resamples
bootCover <- function(x, true = NULL, selectedOnly = FALSE){
  if(is(x, 'resample')){
    niter <- x$call$niter
    y1 <- lapply(x$samples$coefs, '[[', 'lower')
    y2 <- lapply(x$samples$coefs, '[[', 'upper')
  } else if(is(x, 'bootNet')){
    stopifnot('bootFits' %in% names(x))
    stopifnot('fitobj' %in% names(x$bootFits$boot1))
    if(is.null(true)){true <- x$fit0}
    niter <- length(x$bootFits)
    x <- lapply(x$bootFits, '[[', 'fitobj')
    y <- lapply(names(x$boot1), function(z) lapply(x, '[[', z))
    y <- setNames(lapply(names(x$boot1), function(z){
      z0 <- lapply(x, '[[', z)
      z1 <- do.call(rbind, lapply(z0, function(i) setNames(i$lower, i$predictor)))
      z2 <- do.call(rbind, lapply(z0, function(i) setNames(i$upper, i$predictor)))
      return(list(lower = z1, upper = z2))
    }), names(x$boot1))
    y1 <- lapply(y, '[[', 'lower')
    y2 <- lapply(y, '[[', 'upper')
  }
  stopifnot(!is.null(true))
  if(is(true, 'mgmSim')){
    t1 <- net(true)
    v <- paste0(ifelse('X1' %in% colnames(true$data), 'X', 'V'), 1:ncol(t1))
    dimnames(t1) <- rep(list(v), 2)
    if('b2' %in% names(true)){
      t1 <- rbind(t1, true$m1)
      rownames(t1)[nrow(t1)] <- 'M'
      t2 <- netInts(true)
      dimnames(t2) <- list(paste0(v, ':M'), v)
      t1 <- rbind(t1, t2)
    }
    true <- t1
  } else if(is(true, 'ggm')){
    t1 <- net(true, nodewise = TRUE)
    if('interactions' %in% names(true)){
      t1 <- rbind(t1, true$mnet$adjMat[nrow(true$mnet$adjMat), -ncol(true$mnet$adjMat), drop = FALSE])
      t1 <- rbind(t1, netInts(true))
    }
    true <- t1
  }
  out <- matrix(0, ncol = ncol(true), nrow = nrow(true))
  dimnames(out) <- dimnames(true)
  for(i in 1:ncol(true)){
    lo <- y1[[i]]
    hi <- y2[[i]]
    real <- true[-i, i]
    k <- which(names(real) %in% colnames(lo))
    for(j in seq_along(k)){
      lower <- lo[, j]
      upper <- hi[, j]
      if(selectedOnly){
        lower <- na.omit(lower)
        upper <- na.omit(upper)
      } else {
        lower[is.na(lower)] <- 0
        upper[is.na(upper)] <- 0
      }
      ind <- ifelse(k[j] < i, k[j], k[j] + 1)
      out[ind, i] <- sum(lower <= real[k[j]] & real[k[j]] <= upper)/length(lower)
    }
  }
  out
}


### ======================================================================== ###
### ======================================================================== ###
##### plotCoefs: plot coefficients from SURfit with confidence intervals
plotCoefs <- function(fit, true = FALSE, alpha = .05, plot = TRUE, col = "blue", 
                      flip = TRUE, data = NULL, select = TRUE, size = 1, 
                      labels = TRUE, title = NULL){
  if(isTRUE(attr(fit, 'mlGVAR')) & 'varMods' %in% names(fit)){
    if(class(true) %in% c('logical', 'character', 'numeric')){
      if(is.numeric(true)){true <- as.logical(true - 1)}
      true <- ifelse(is.logical(true), c('fixed', 'between')[true + 1], match.arg(
        true, c('fixed', 'between', 'temporal', 'contemporaneous', 'ggm')))
    } else stop('Must provide resample output to include true network')
    true <- switch(true, ggm = , between = 'between', 'fixed')
    fit <- fit$varMods[[grep(true, names(fit$varMods))]]
    if(is.null(title) & !identical(title, FALSE)){
      title <- switch(true, fixed = 'Temporal Network', 'Between-subjects network')
    }
    true <- FALSE
  }
  if(is(fit, "systemfit")){
    x <- as.matrix(confint(fit, level = (1 - alpha)))
    x <- cbind(x, coef(fit))
    x <- data.frame(x)
    colnames(x) <- c("lower", "upper", "b")
    x <- x[!grepl("Intercept", rownames(x)),]
    x$predictor <- rownames(x)
    rownames(x) <- 1:nrow(x)
    x$Y <- NA
    for(i in 1:length(fit$eq)){
      x[grepl(paste0("eq", i), x$predictor), "Y"] <- as.character(fit$eq[[i]]$terms[[2]])
    }
    x$Y <- as.factor(x$Y)
    x$predictor <- as.factor(qdapRegex::rm_between(x$predictor, "eq", "_"))
    x <- data.frame(Y = x$Y, predictor = x$predictor, lower = x$lower, b = x$b, upper = x$upper)
    if(!is.null(true)){
      dat <- tryCatch({eval(fit$call$data)}, error = function(e){
        if(!is.null(data)){data} else {stop("Need to supply data to index names of predictors")}
      })
      if(class(dat) == "list"){dat <- cbind.data.frame(dat$X, dat$Y)} # FULLFIX
      nY <- rep(levels(x$Y), each = (ncol(true) - 1))
      nP <- rep(colnames(dat)[!colnames(dat) %in% levels(x$Y)], length(levels(x$Y)))
      trueDat <- data.frame(nY, nP, B = as.vector(t(true[,-1])))
    }
    dat <- list()
    for(i in 1:length(fit$eq)){dat[[i]] <- x[grepl(paste0("^X", i, ".y$"), x$Y),]}
    dat <- do.call(rbind, dat)
  } else {
    if(all(c("mod0", "mod1se") %in% names(fit))){
      stop("Must index either 'mod0' or 'mod1se' to indicate which to use")
    } else if(any(grepl("CIs|freqs", names(fit)))){
      if(is.null(true) & !any(grepl("freqs", names(fit)))){
        true <- if("adjCIs" %in% names(fit)){fit$fits$fit1} else {fit$fit1}} 
      fit <- ifelse("adjCIs" %in% names(fit), list(fit$adjCIs), 
                    ifelse("freqs" %in% names(fit), list(fit$freqs), 
                           list(fit$fitCIs)))[[1]]
    }
    if(is.logical(true)){true <- NULL}
    ynames <- names(fit)
    x <- do.call(rbind, fit)
    if("select" %in% colnames(fit[[1]])){
      x <- cbind.data.frame(Y = rep(ynames, sapply(fit, nrow)), x)
    } else {
      x <- cbind.data.frame(Y = rep(ynames, each = length(levels(x$predictor))), x)
    }
    if(!is.null(true)){
      if(is.list(true)){
        true2 <- unlist(lapply(lapply(true$mods, '[[', "model"), function(z) z[-1, ]))
        if(length(as.vector(true2)) != nrow(x)){
          true3 <- names(true2)
          true4 <- rep(NA, nrow(x))
          xx <- unname(apply(x[, 1:2], 1, paste0, collapse = "."))
          true4[match(true3, xx)] <- true2
          true2 <- true4
        }
        trueDat <- data.frame(nY = x$Y, nP = as.character(x$predictor), B = as.vector(unname(true2)))
      } else if(!is.logical(true)){
        if(ncol(true) == (2 * nrow(true))){true <- true[, -1]}
        trueDat <- data.frame(nY = x$Y, nP = as.character(x$predictor), B = as.vector(t(true)))
      }
    }
    if(any(rowSums(is.na(x)) != 0)){x <- x[rowSums(is.na(x)) == 0, ]}
    x$Y <- factor(x$Y)
    x$predictor <- factor(x$predictor)
    dat <- x
  }
  rownames(dat) <- 1:nrow(dat)
  dat$ord <- paste(dat$Y, dat$predictor, sep = "_")
  if(!is.null(true)){
    trueDat$ord <- paste(trueDat$nY, trueDat$nP, sep = "_")
    dat$B <- trueDat[trueDat$ord %in% dat$ord, "B"]
  }
  if(select != FALSE & "select" %in% names(dat)){
    if(!isTRUE(select)){dat$select <- !(dat$lower <= 0 & dat$upper >= 0)}
    dat <- dat[dat$select, ]
  }
  if(plot == TRUE){
    plotCI <- function(dat, xlabs, true, flip = TRUE, size = 1, labels = TRUE, title = NULL){
      invisible(suppressMessages(require(ggplot2)))
      if(is.logical(labels)){if(!labels){xlabs <- NULL}} else {xlabs <- labels}
      p <- ggplot(dat, aes(x = reorder(ord, b), y = b, group = Y)) +
        facet_wrap(Y ~., scales = "free") + geom_point(size = size, na.rm = TRUE) + 
        geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
        geom_hline(aes(yintercept = 0), linetype = "dashed") +
        scale_x_discrete(labels = xlabs) + xlab("Predictor") + 
        ylab(expression(hat(beta))) + theme_bw()
      if(!is.null(true)){
        p <- p + geom_point(data = dat, aes(x = reorder(ord, B), y = B, group = Y),
                            col = col, alpha = .85, shape = 8, na.rm = TRUE, size = size)
      }
      if(!is.null(title)){p <- p + ggtitle(title)}
      if(flip == TRUE){return(p + coord_flip())} else {return(p)}
    }
    plotCI(dat, xlabs = function(x){sub("[^_]*_", "", x)}, true = true, 
           flip = flip, size = size, labels = labels, title = title)
  } else {
    return(dat)
  }
}

##### plotPvals: plot the ECDF of pvalues from resampler
plotPvals <- function(obj, outcome = 1, predictor = 1, alpha = .05){
  stopifnot("adjCIs" %in% names(obj))
  pvals <- lapply(obj$samples$coefs, '[[', "P")
  if(is.character(outcome)){outcome <- which(names(pvals) == outcome)}
  if(is.character(predictor)){predictor <- which(colnames(pvals[[outcome]]) == predictor)}
  ps <- pvals[[outcome]][, predictor]
  n <- length(ps)
  gammas <- seq(ceiling(alpha * n)/n, 1 - 1/n, by = 1/n)
  h <- hist(ps, breaks = n, plot = FALSE)
  h$density <- h$density/sum(h$density)
  fdr_line <- function(p, gammas, alpha = .05){pmax(.05, (1 - (log(min(gammas)))/alpha) * p)}
  plot(h, freq = FALSE, main = "", ylim = c(0, 1), xlab = "Adjusted P-values")
  lines(sort(ps), fdr_line(sort(ps), gammas, alpha), lty = 2, lwd = 2)
  lines(sort(ps), seq(.05, 1, length = n), lwd = 2)
}

##### modSelect: creates necessary input for fitNetwork when selecting variables
modSelect <- function(obj, data = NULL, fit = FALSE, select = "select", 
                      thresh = NULL, ascall = TRUE, type = "gaussian", ...){
  if(is.null(data)){ascall <- FALSE}
  cis <- c("adjCIs", "freqs")[which(c("adjCIs", "freqs") %in% names(obj))]
  if(length(cis) != 0){obj <- obj[[cis]]}
  allNames <- lapply(lapply(obj, '[[', "predictor"), as.character)
  vs <- names(allNames) <- names(obj)
  if(is.numeric(select)){
    for(i in seq_along(obj)){obj[[i]]$select <- obj[[i]]$freq >= select}
    select <- "select"
  } else if(!is.null(thresh)){
    if(grepl('select|ci', select)){select <- 'freq'}
    select <- match.arg(select, c("freq", "split1", "split2", "Pvalue"))
    for(i in seq_along(obj)){
      if(select != "Pvalue"){
        obj[[i]]$select <- obj[[i]][[select]] >= thresh
      } else {
        obj[[i]]$select <- obj[[i]][[select]] <= thresh
      }
    }
    select <- "select"
  } else if(select == "ci"){
    select <- "select_ci"
  }
  if(select == 'select_ci'){
    sel <- sapply(obj, function(z) all(z$select_ci == z$select))
    if(all(sel)){message('CIs lead to same decisions as p-values')}
  }
  varMods <- lapply(seq_along(obj), function(z){
    list(mod0 = as.character(obj[[z]][which(obj[[z]][[select]]), "predictor"]))
  })
  if(any(sapply(lapply(varMods, '[[', "mod0"), length) == 0)){
    v0 <- which(sapply(lapply(varMods, '[[', "mod0"), length) == 0)
    for(jj in seq_along(v0)){varMods[[v0[jj]]]$mod0 <- "1"}
  }
  attributes(varMods) <- attributes(obj)
  type <- ifelse("type" %in% names(attributes(obj)), 
                 list(attr(obj, "type")), list(type))[[1]]
  for(i in seq_along(obj)){
    if(any(grepl(":", varMods[[i]]$mod0))){
      ints <- varMods[[i]]$mod0[grepl(":", varMods[[i]]$mod0)]
      mains <- unique(unlist(strsplit(ints, ":")))
      mains <- mains[order(match(mains, vs))]
      varMods[[i]]$mod0 <- union(varMods[[i]]$mod0, union(mains, ints))
      varMods[[i]]$mod0 <- varMods[[i]]$mod0[match(allNames[[i]], varMods[[i]]$mod0)]
      varMods[[i]]$mod0 <- varMods[[i]]$mod0[!is.na(varMods[[i]]$mod0)]
    }
    attr(varMods[[i]], "family") <- ifelse(
      length(type) == 1, ifelse(type %in% c("g", "gaussian"), "g", "c"), 
      ifelse(type[i] %in% c("g", "gaussian"), "g", "c")
    )
  }
  if(!is.null(data)){
    if(!is(data, "data.frame")){
      data <- data.frame(data, check.names = FALSE)
    }
    args <- list(...)
    atts <- attributes(obj)
    atts[c("names", "criterion", "method")] <- NULL
    atts$data <- data
    atts$type <- varMods
    if("lags" %in% names(attributes(atts$type))){
      if(attr(atts$type, "lags") == 1 & "moderators" %in% names(attributes(atts$type))){
        m <- attr(atts$type, "moderators")
        attr(atts$type, "moderators") <- colnames(data)[m]
        attr(varMods, "moderators") <- colnames(data)[m]
      }
    }
    if("err" %in% names(atts)){atts$err <- NULL}
    dupArgs <- intersect(names(args), names(atts))
    newArgs <- setdiff(names(args), names(atts))
    attr(varMods, "call") <- atts <- append(
      replace(atts, dupArgs, args[dupArgs]), args[newArgs])
    if(fit){return(do.call(fitNetwork, atts))}
  }
  if(ascall){return(attr(varMods, "call"))} else {return(varMods)}
}

##### getFitCIs: provides model coefficients with CIs
getFitCIs <- function(fit, allNames = NULL, alpha = .05){
  if("SURnet" %in% names(fit)){
    if(!'SURfit' %in% names(fit)){stop('Requires SURfit')}
    fit <- append(fit$SURnet, list(fitobj = fit$SURfit$eq))
  }
  if(is.null(allNames)){
    allNames <- lapply(fit$mods, function(z) rownames(z$model)[-1])
  } else if(all(c("mods", "call") %in% names(allNames))){
    allNames <- lapply(allNames$mods, function(z) rownames(z$model)[-1])
  }
  fitobj <- fit$fitobj
  fitCoefs <- suppressWarnings(suppressMessages(
    lapply(seq_along(fitobj), function(z){
      coefs1 <- summary(fitobj[[z]])$coefficients[-1, , drop = FALSE]
      if(nrow(coefs1) == 0){return(NA)}
      coefs2 <- confint(fitobj[[z]], level = 1 - alpha)[-1, , drop = FALSE]
      coefs <- cbind.data.frame(coefs1, coefs2)
      colnames(coefs) <- c("b", "se", "t", "P", "lower", "upper")
      coefs <- coefs[match(allNames[[z]], rownames(coefs)), ]
      coefs <- data.frame(predictor = factor(allNames[[z]]), lower = coefs$lower,
                          b = coefs$b, upper = coefs$upper, Pvalue = coefs$P,
                          select = ifelse(is.na(coefs$b), FALSE, TRUE))
      rownames(coefs) <- 1:nrow(coefs)
      if(any(is.na(coefs))){coefs[is.na(coefs)] <- NA}
      return(coefs)
    })))
  names(fitCoefs) <- names(fitobj)
  return(fitCoefs)
}

##### adjustedCI: uses the union-bound approach for obtaining adjCIs
adjustedCI <- function(lci, rci, centers, ses, df.res, gamma.min, ci.level, var, 
                       multi.corr = FALSE, verbose = FALSE, s0 = list(s0 = NA)){
  findInsideGamma <- function(low, high, ci.info, verbose){
    range.length <- 10
    test.range <- seq(low, high, length.out = range.length)
    cover <- mapply(doesItCoverGamma, beta.j = test.range, ci.info = list(ci.info = ci.info))
    while(!any(cover)){
      range.length <- 10 * range.length
      test.range <- seq(low, high, length.out = range.length)
      cover <- mapply(doesItCoverGamma, beta.j = test.range, ci.info = list(ci.info = ci.info))
      if(range.length > 10^3){
        message("FOUND NO INSIDE POINT")
        message("number of splits")
        message(length(ci.info$centers))
        message("centers")
        message(ci.info$centers)
        message("ses")
        message(ci.info$ses)
        message("df.res")
        message(ci.info$df.res)
        return(NULL)
        #stop("couldn't find an inside point between low and high. The confidence interval doesn't exist!")
      }
    }
    if(verbose){cat("Found an inside point at granularity of", range.length, "\n")}
    min(test.range[cover])
  }
  doesItCoverGamma <- function(beta.j, ci.info){
    if(missing(ci.info)){stop("ci.info is missing to the function doesItCoverGamma")}
    centers <- ci.info$centers
    ci.lengths <- ci.info$ci.lengths
    no.inf.ci <- ci.info$no.inf.ci
    ses <- ci.info$ses
    df.res <- ci.info$df.res
    gamma.min <- ci.info$gamma.min
    multi.corr <- ci.info$multi.corr
    s0 <- ci.info$s0
    alpha <- 1 - ci.info$ci.level
    pval.rank <- rank(-abs(beta.j - centers)/(ci.lengths/2))
    nsplit <- length(pval.rank) + no.inf.ci
    gamma.b <- pval.rank/nsplit
    if(multi.corr){
      if(any(is.na(s0))){stop("need s0 information to be able to create multiple testing corrected pvalues")}
      level <- (1 - alpha * gamma.b/(1 - log(gamma.min) * s0))
    } else {
      level <- (1 - alpha * gamma.b/(1 - log(gamma.min)))
    }
    a <- (1 - level)/2
    a <- 1 - a
    if(all(gamma.b <= gamma.min)){
      TRUE
    } else {
      fac <- qt(a, df.res)
      nlci <- centers - fac * ses
      nrci <- centers + fac * ses
      all(nlci[gamma.b > gamma.min] <= beta.j) && all(nrci[gamma.b > gamma.min] >= beta.j)
    }
  }
  bisectionBounds <- function(shouldcover, shouldnotcover, ci.info, verbose){
    reset.shouldnotcover <- FALSE
    if(doesItCoverGamma(beta.j = shouldnotcover, ci.info = ci.info)){
      reset.shouldnotcover <- TRUE
      if(verbose){cat("finding a new shouldnotcover bound\n")}
      while(doesItCoverGamma(beta.j = shouldnotcover, ci.info = ci.info)){
        old <- shouldnotcover
        shouldnotcover <- shouldnotcover + (shouldnotcover - shouldcover)
        shouldcover <- old
      }
      if(verbose){
        cat("new\n")
        cat("shouldnotcover", shouldnotcover, "\n")
        cat("shouldcover", shouldcover, "\n")
      }
    }
    if(!doesItCoverGamma(beta.j = shouldcover, ci.info = ci.info)){
      if(reset.shouldnotcover){stop("Problem: we first reset shouldnotcover and are now resetting shouldcover, this is not supposed to happen")}
      if(verbose){cat("finding a new shouldcover bound\n")}
      while(!doesItCoverGamma(beta.j = shouldcover, ci.info = ci.info)){
        old <- shouldcover
        shouldcover <- shouldcover + (shouldcover - shouldnotcover)
        shouldnotcover <- old
      }
      if(verbose){
        cat("new\n")
        cat("shouldnotcover", shouldnotcover, "\n")
        cat("shouldcover", shouldcover, "\n")
      }
    }
    return(list(shouldcover = shouldcover, shouldnotcover = shouldnotcover))
  }
  bisectionCoverage <- function(outer, inner, ci.info, verbose, eps.bound = 10^(-7)){
    checkBisBounds(shouldcover = inner, shouldnotcover = outer, ci.info = ci.info, verbose = verbose)
    eps <- 1
    while(eps > eps.bound){
      middle <- (outer + inner)/2
      if(doesItCoverGamma(beta.j = middle, ci.info = ci.info)){
        inner <- middle
      } else {
        outer <- middle
      }
      eps <- abs(inner - outer)
    }
    solution <- (inner + outer)/2
    if(verbose){cat("finished bisection...eps is", eps, "\n")}
    return(solution)
  }
  checkBisBounds <- function(shouldcover, shouldnotcover, ci.info, verbose){
    if(doesItCoverGamma(beta.j = shouldnotcover, ci.info = ci.info)){
      stop("shouldnotcover bound is covered! we need to decrease it even more! (PLZ implement)")
    } else if(verbose){cat("shouldnotcover bound is not covered, this is good")}
    if(doesItCoverGamma(beta.j = shouldcover, ci.info = ci.info)){
      if(verbose){cat("shouldcover is covered!, It is a good covered bound")}
    } else {
      stop("shouldcover is a bad covered bound, it is not covered!")
    }
  }
  inf.ci <- is.infinite(lci) | is.infinite(rci)
  no.inf.ci <- sum(inf.ci)
  if(verbose){cat("number of Inf ci:", no.inf.ci, "\n")}
  if((no.inf.ci == length(lci)) || (no.inf.ci >= (1 - gamma.min) * length(lci))){return(c(-Inf, Inf))}
  lci <- lci[!inf.ci]
  rci <- rci[!inf.ci]
  centers <- centers[!inf.ci]
  df.res <- df.res[!inf.ci]
  ses <- ses[!inf.ci]
  s0 <- s0[!inf.ci]
  ci.lengths <- rci - lci
  ci.info <- list(lci = lci, rci = rci, centers = centers, ci.lengths = ci.lengths, 
                  no.inf.ci = no.inf.ci, ses = ses, s0 = s0, df.res = df.res, 
                  gamma.min = gamma.min, multi.corr = multi.corr, ci.level = ci.level)
  inner <- findInsideGamma(low = min(centers), high = max(centers), ci.info = ci.info, verbose = verbose)
  if(is.null(inner)){
    alpha <- 1 - ci.level
    beta.j <- seq(min(centers), max(centers), length = length(centers))
    pval.rank <- rank(-abs(beta.j - centers)/(ci.lengths/2))
    nsplit <- length(pval.rank) + no.inf.ci
    gamma.b <- pval.rank/nsplit
    level <- (1 - alpha * gamma.b/(1 - log(gamma.min)))
    a <- 1 - ((1 - level)/2)
    fac <- qt(a, df.res)
    nlci <- centers - fac * ses
    nrci <- centers + fac * ses
    return(c(0, 0))
  }
  outer <- min(lci)
  new.bounds <- bisectionBounds(shouldcover = inner, shouldnotcover = outer, ci.info = ci.info, verbose = verbose)
  inner <- new.bounds$shouldcover
  outer <- new.bounds$shouldnotcover
  l.bound <- bisectionCoverage(outer = outer, inner = inner, ci.info = ci.info, verbose = verbose)
  if(verbose){cat("lower bound ci aggregated is", l.bound, "\n")}
  outer <- max(rci)
  new.bounds <- bisectionBounds(shouldcover = inner, shouldnotcover = outer, ci.info = ci.info, verbose = verbose)
  inner <- new.bounds$shouldcover
  outer <- new.bounds$shouldnotcover
  u.bound <- bisectionCoverage(inner = inner, outer = outer, ci.info = ci.info, verbose = verbose)
  if(verbose){cat("upper bound ci aggregated is", u.bound, "\n")}
  return(c(l.bound, u.bound))
}


### ======================================================================== ###
### ======================================================================== ###
##### stability: get empirical selection probabilities from resampled data
stability <- function(obj){
  objcoefs <- obj$freqs
  p <- length(objcoefs)
  vs <- names(objcoefs)
  allNames <- lapply(lapply(objcoefs, '[[', "predictor"), as.character)
  k <- sapply(allNames, length)
  niter <- obj$call$niter
  nlam <- obj$call$nlam - 1
  splits <- crits <- list()
  for(i in 1:2){
    s0 <- paste0("split", i)
    splits[[i]] <- lapply(1:p, function(z){
      pIts <- lapply(1:niter, function(it){
        lapply(1:k[z], function(j){
          sapply(obj[[s0]][[it]]$varMods[[z]]$allCoefs, '[[', j)
        })
      })
      pOut <- lapply(1:k[z], function(zz){
        pCoefs <- do.call(rbind, lapply(pIts, '[[', zz))
        return(colMeans(abs(sign(pCoefs))))
      })
      pOut2 <- lapply(1:k[z], function(zz){
        pCoefs2 <- do.call(rbind, lapply(pIts, '[[', zz))
        return(abs(sign(pCoefs2)))
      })
      freqOut <- do.call(cbind.data.frame, pOut)
      names(freqOut) <- allNames[[z]]
      return(list(freqs = freqOut, signs = pOut2))
    })
    crits[[i]] <- lapply(1:p, function(z){
      sapply(1:niter, function(zz){
        obj[[s0]][[zz]]$varMods[[z]]$fitobj[[3]]
      })
    })
    names(crits[[i]]) <- names(splits[[i]]) <- vs
  }
  simulVars <- lapply(1:p, function(z){
    simul1 <- sapply(1:k[z], function(zz){
      colMeans(splits[[1]][[z]]$signs[[zz]] * splits[[2]][[z]]$signs[[zz]])
    })
    simul1 <- data.frame(simul1)
    colnames(simul1) <- allNames[[z]]
    return(simul1)
  })
  names(simulVars) <- vs
  for(i in 1:2){for(j in 1:p){splits[[i]][[j]] <- splits[[i]][[j]]$freqs}}
  output <- list(split1 = append(splits[[1]], list(crits = crits[[1]])), 
                 split2 = append(splits[[2]], list(crits = crits[[2]])), 
                 simult = simulVars)
  return(output)
}

##### plotStability: plot stability paths for a given variable
plotStability <- function(obj, pp = 1, s = 3, thresh = .5, color = "black"){
  if(obj$call$criterion == "CV"){stop("Not possible when criterion = CV")}
  objcall <- obj$call
  obj <- obj$stability
  node <- names(obj[[s]])[pp]
  p <- length(obj[[s]]) - ifelse(s == 3, 0, ifelse(!objcall$exogenous, 2, 1))
  stopifnot(pp <= p + 1)
  n <- nrow(obj[[s]][[pp]])
  k <- ncol(obj[[s]][[pp]])
  plot(0, type = "n", ylim = c(0, 1), xlim = c(1, n + 1), axes = F,
       main = node, xlab = expression(lambda), ylab = "Selection Probability") 
  axis(1); axis(2)
  if(any(objcall$moderators == 0)){
    colors <- rep("black", k)
    lty <- rep(1, k)
  } else {
    colors <- c(rep("black", p), rep("grey50", k - p))
    lty <- c(rep(1, p), rep(2, k - p))
  }
  if(color == "terrain"){colors <- terrain.colors(k)}
  if(!is.null(thresh)){
    if(thresh != 0){
      lines(1:n, rep(thresh, n), lwd = 2, lty = ifelse(color == "terrain", 2, 1),
            col = ifelse(color == "terrain", "black", rgb(1, 0, 0, .2)))
    }
  }
  #if(!is.null(thresh)){if(thresh != 0){abline(h = thresh, col = rgb(1, 0, 0, .2), lwd = 2)}}
  for(i in 1:k){
    lines(obj[[s]][[pp]][, i], lwd = 2, col = colors[i], lty = lty[i])
  }
  if(color == "terrain"){ 
    legend(n + 3, 1.1, xpd = T, legend = colnames(obj[[s]][[pp]]), 
           col = colors, lty = lty, lwd = 2, bty = "n")
  }
}


### ======================================================================== ###
### ======================================================================== ###
if(exists("clear")){
  if(isTRUE(clear)){
    message("Clearing .GlobalEnv")
    rm(clear)
    .modnets <- new.env()
    .modnets <- globalenv()
    suppressMessages(attach(.modnets))
    rm(list = setdiff(ls(), ifelse(exists('keep'), 'keep', '')))
  }
}
message("This is modnets 1.0.0")
