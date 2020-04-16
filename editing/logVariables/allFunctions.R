if(!exists("pkgs")){pkgs <- c("glmnet", "mgm", "qgraph", "grpregOverlap", "systemfit", "leaps", "hierNet", "glinternet")}
if(!"none" %in% pkgs){invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))}
if(exists("clear")){if(clear == TRUE){rm(list = ls()); clear <- TRUE}} else {rm(pkgs)}
if(".modnets" %in% search()){detach(".modnets")}
setwd('~/C: Desktop/COMPS/METHODS/CODE/modnets')
files <- paste0('./', c('ggm', 'testing', 'centrality', 'sim', 'mlGVAR', 'simGVAR', 'penalized', 'power', 'plots'),'.R')
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
                "m3", "m5", "mod", "PTSD", "fried", "wichers", "new", "dep1", 
                "dep2", "bfi1", "bfi2", "bfi3", "obama", "bfiDat", "constantini")
  if(d){return(data.frame(sort(whatData)))}
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
                       getLL = TRUE, saveMods = TRUE, binarize = FALSE, ...){
  t1 <- Sys.time() # START
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
  ### SETUP
  if(is.null(lambda)){
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
      if(!saveMods){output$fitobj <- output$mods0$models <- NULL}
      if(verbose){print(Sys.time() - t1)}
      return(output)
    } else {
      ### TEMPORAL
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
                    covs = covariates, sur = std, maxiter = maxiter)
      dat <- lagMat(data = data, type = type, m = moderators, covariates = covariates, 
                    center = center, scale = scale, exogenous = exogenous)
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
      if(!saveMods){output$SURfit <- NULL}
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
                    binarize = FALSE, ...){
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
    if("SURnet" %in% names(object)){object <- object$SURnet}
    if(!"adjMat" %in% names(object)){
      object[c("adjMat", "edgeColors")] <- eval(parse(text = paste0("object$", switch(
        which.net, pdc = "temporal$PDC", which.net))))[c("adjMat", "edgeColors")]
      if(!startsWith(which.net, "c") & "modEdges" %in% names(object$temporal)){
        object$modEdges <- object$temporal$modEdges
      }
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
        object$mnet$adjMat[ind, -ind] <- object$mnet$adjMat[ind, -ind] * ifelse(
          object$mods0$Bm[, 4] <= threshold, 1, 0)
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
      if("modEdges" %in% names(object)){lty <- t(object$modEdges)} else {lty <- 1}
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
        if("modEdges" %in% names(object)){lty <- object$modEdges} else {lty <- 1}
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
      if("modEdgesNW" %in% names(object$nodewise)){lty <- object$nodewise$modEdgesNW} else {lty <- 1}
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
avlay <- function(..., net = 'temporal', collapse = FALSE, args = NULL){
  x <- list(...)
  stopifnot(length(x) > 0)
  if(length(x) == 1){x <- x[[1]]} else if(collapse){x <- appd(x)}
  args0 <- list(object = NA, which.net = net, plot = FALSE)
  args <- append(args[setdiff(names(args), names(args0))], args0)
  averageLayout(lapply(x, function(z) do.call(
    plotNet, replace(args, 'object', list(object = z)))))
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
                   covs = NULL, sur = TRUE, ...){
  if(!is(varMods, 'list')){type <- varMods; varMods <- NULL}
  eqs <- surEqs(data = data, varMods = varMods, mod = match.arg(mod, c("min", "1se")), 
                m = m, exogenous = exogenous, covs = covs)
  dat <- lagMat(data = data, type = type, m = m, covariates = covs, center = center, 
                scale = scale, exogenous = exogenous)
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
##### lagMat: create lagged matrices for fitting SUR models
lagMat <- function(data, type = "g", m = NULL, covariates = NULL, center = TRUE,
                   scale = FALSE, exogenous = TRUE, lags = 1, checkType = FALSE){
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
  if(is(fit, 'bootNet') | (isTRUE(attr(fit, 'resample')) & 'fit0' %in% names(fit))){fit <- fit$fit0}
  if(is(fit, 'mgmSim')){fit <- fit[[grep('^trueNet$|^b1$', names(fit))]]}
  if(is(fit, "matrix")){return(fit)}
  if(is(fit, "mgm")){
    return(fit$pairwise$wadj * replace(
      fit$pairwise$signs, is.na(fit$pairwise$signs), 0))
  }
  if(isTRUE(n)){n <- threshold <- "beta"}
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
                    rule = 'none', r = NULL){
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
    if(is(out, 'list')){return(list())} else {return(out)}
  }
  if(is(fit, 'mgmSim')){if('b2' %in% names(fit)){return(fit$b2)} else {return(list())}}
  if(is(fit, 'bootNet') | (isTRUE(attr(fit, 'resample')) & 'fit0' %in% names(fit))){fit <- fit$fit0}
  if(isTRUE(n) & isTRUE(threshold)){avg <- TRUE}
  if(isTRUE(n)){
    if(is.character(threshold)){
      if(tolower(threshold) %in% c('and', 'or')){
        rule <- threshold
      }
    }
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
  } else {
    if(!'interactions' %in% names(fit)){return(list())}
    out <- fit$interactions[[1]]
    if(isTRUE(attr(fit, "ggm"))){
      pvals <- out[[2]][[2]]
      out <- out[[2]][[1]]
    } else {
      pvals <- fit$interactions$pvals
    }
  }
  if(threshold != FALSE & !"mlVARsim" %in% atts){
    rule <- match.arg(tolower(rule), c('none', 'or', 'and'))
    if(!is.numeric(threshold)){threshold <- .05}
    if(rule == 'none'){
      out <- out * ifelse(pvals <= threshold, 1, 0)
    } else if(rule == 'or'){
      out <- out * ifelse(pvals <= threshold | t(pvals) <= threshold, 1, 0)
    } else if(rule == 'and'){
      out <- out * ifelse(pvals <= threshold & t(pvals) <= threshold, 1, 0)
    }
  }
  if(avg){out <- (t(out) + out)/2}
  if(!is.null(r)){out <- out[-r, -r]}
  if(is(out, 'list')){return(list())}
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
                         outMsgs = FALSE, dmnames = NULL){
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
                          cat("\n")
                          while(failed == TRUE){
                            if(take <= 5){
                              cat("  Failed.. trying again, take =", take, "\n")
                              fitCV <- try(glinternet.cv(x, y, type, nFolds = nfolds, nLambda = nlam, 
                                                         interactionCandidates = m, family = fam), silent = TRUE)
                              if(class(fitCV) == "try-error"){
                                failed <- TRUE
                                take <- take + 1
                              } else {
                                failed <- FALSE
                              }
                            } else if(take <= 10){
                              cat("  Failed.. trying nlam = 20, take =", take, "\n")
                              fitCV <- try(glinternet.cv(x, y, type, nFolds = nfolds, nLambda = 20, 
                                                         interactionCandidates = m, family = fam), silent = TRUE)
                              if(class(fitCV) == "try-error"){
                                failed <- TRUE
                                take <- take + 1
                              } else {
                                failed <- FALSE
                              }
                            } else {
                              cat("  Failed.. trying nlam = 20 & nFolds = 3, take =", take, "\n")
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
        cat("\n  Mod0 empty.. choosing new lambda\n")
        which.lam0 <- which.lam0 + 1
      }
      fitCV$lambdaHat <- fitCV$lambda[which.lam0]
      which.lam1se <- which(fitCV$lambda == fitCV$lambdaHat1Std)
      while(is.null(fitCV$glinternetFit$activeSet[[which.lam1se]])){
        cat("\n  Mod1SE empty.. choosing new lambda\n")
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
                      useSE = TRUE, nlam = NULL, covs = NULL, verbose = TRUE){
  dat <- data
  ALL <- FALSE
  dmnames <- NULL # VARSELECT START
  if(is.null(lags)){
    if(!is.null(m)){if(length(m) >= ncol(dat) - 1){exogenous <- FALSE}}
    if(!exogenous){ALL <- TRUE}
  }
  if(ifelse(is.null(m), FALSE, ifelse(
    all(m == 0) | all(m == "all"), TRUE, FALSE))){
    m <- 0
    #if(!is.null(lags)){m <- NULL} else if(m == "all"){m <- 0}
    ALL <- TRUE
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
                                    criterion = criterion, dmnames = dmnames)
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
resample <- function(data, m = NULL, niter = 10, sampMethod = "split", criterion = "AIC",
                     method = "glmnet", rule = "OR", gamma = .5, nfolds = 10, 
                     nlam = 50, which.lam = "min", threshold = FALSE, bonf = FALSE, 
                     alpha = .05, exogenous = TRUE, split = .5, center = TRUE, 
                     scale = FALSE, varSeed = NULL, seed = NULL, verbose = TRUE, 
                     lags = NULL, binary = NULL, type = 'g', saveMods = TRUE,
                     saveData = FALSE, saveVars = FALSE, fitit = TRUE, 
                     nCores = 1, cluster = 'mclapply', ...){
  t1 <- Sys.time() # RESAMPLE START
  if(!is.null(lags)){lags <- switch(2 - identical(as.numeric(lags), 0), NULL, 1)}
  if(!is.null(m)){if(all(m == 0)){m <- NULL}}
  method <- ifelse(!is.null(m), 'glinternet', ifelse(
    !method %in% c('glmnet', 'subset'), 'glmnet', method))
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


### ======================================================================== ###
### ============================= BOOTSTRAPPING ============================ ###
### ======================================================================== ###
##### bootNet: essentially 'bootnet', but with fitNetwork for all estimation
bootNet <- function(data, m = NULL, nboots = 10, lags = NULL, caseDrop = FALSE, rule = 'OR',
                    ci = .95, caseMin = .05, caseMax = .75, caseN = 10, threshold = FALSE,
                    fits = NULL, type = 'g', saveMods = TRUE, verbose = TRUE,
                    size = NULL, nCores = 1, cluster = 'mclapply', ...){
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
                  threshold = sampThresh, verbose = FALSE, saveMods = FALSE)
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
          parallel::clusterExport(cl, c('fitNetwork', 'Matrix', 'net', 'netInts', obj1), envir = environment())
        }
      } else {
        cl <- nCores
      }
      if(verbose){
        pbapply::pboptions(type = 'timer', char = '-')
        fits <- suppressWarnings(structure(pbapply::pblapply(seq_len(nboots), function(z){
          args0$data <- switch(lags + 1, data[inds[[z]], ], structure(data, samp_ind = inds[[z]]))
          fit <- tryCatch({do.call(fitNetwork, args0)}, error = function(e){
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
  strFUN <- switch(2 - net3, function(x){diag(x) <- 0; rowSums(x)}, function(x){diag(x) <- 0; colSums(x)})
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
                     textPos = 'value', ...){
  args <- tryCatch({list(...)}, error = function(e){list()})
  fit0 <- switch(2 - ('fit0' %in% names(obj)), obj$fit0, list())
  if(isTRUE(attr(obj, 'resample'))){
    if(!'data' %in% names(obj)){return(plotCoefs(fit = obj, ...))}
    if(length(fit0) == 0){
      fit0 <- switch(2 - ('fit0' %in% names(args)), args$fit0, TRUE)
    }
    obj <- do.call(bootNet, append(list(data = obj, fit0 = fit0), replace(args, 'fit0', NULL)))
  }
  if(!identical(threshold, FALSE) & 'bootFits' %in% names(obj)){
    dat <- paste0('obj$fit0$', ifelse('temporal' %in% names(obj), 'SURnet$data', 'data'))
    if(isTRUE(attr(obj$bootFits, 'resample'))){attributes(obj$bootFits)$resample <- NULL}
    obj <- bootNet(data = eval(parse(text = dat)), fits = obj$bootFits, 
                   threshold = threshold, fit0 = obj$fit0)
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
    if("pairwise" %in% names(obj)){
      boots <- obj$boots
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
  if(!lags & !identical(text, FALSE) & type == 'edges'){
    boots <- boots[[switch(2 - interactions, 'boot_int_edges', 'boot_edges')]]
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
                      alloutput = TRUE, ...){
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
                              AND = (rule == 'and'), gamma = gamma))
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
intsPlot <- function(out, y, moderator = NULL, nsims = 500, alpha = .05){
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


### ======================================================================== ###
### ========================== SIMNETS FUNCTIONS =========================== ###
### ======================================================================== ###
##### SimNets: simulate multilevel GVAR data to evaluate parameter recovery
SimNets <- function(N, nsims = 10, ntime = 50, p = 3, m = "random", FUN = "mlGVAR",
                    inds = c("cosine", "mse", "mae"), saveMods = FALSE, 
                    betweenArgs = list(moderators = NULL), d = 14,
                    verbose = TRUE, ...){
  t1 <- Sys.time()
  A0 <- list(...)
  A1 <- append(list(nTime = ntime, nPerson = N, nNode = p, m = m), A0)
  A1 <- A1[intersect(formalArgs("mlGVARsim"), names(A1))]
  A2 <- append(list(data = NULL, m = p + 1, betweenArgs = betweenArgs, verbose = FALSE), A0)
  FUN <- match.arg(FUN, c("mlGVAR", "lmerVAR"))
  FUNargs <- setdiff(formalArgs(FUN), "...")
  if("selectFUN" %in% names(A2) & FUN == "mlGVAR"){
    if(isTRUE(A2$selectFUN)){A2$selectFUN <- "varSelect"}
    FUNargs <- setdiff(union(formalArgs(switch(match.arg(
      A2$selectFUN[1], c("varSelect", "resample", "stability", "bootstrap")), 
      varSelect = "varSelect", "resample")), FUNargs), "...")
  }
  saveVars <- isTRUE("sampMethod" %in% FUNargs)
  A2 <- A2[intersect(FUNargs, names(A2))]
  if(is.null(m)){A2$m <- NULL}
  nit <- paste0("iter", 1:nsims)
  allNets <- c("temporal", "contemporaneous", "between")
  output <- setNames(lapply(seq_along(N), function(z){
    A1$nPerson <- N[z]
    if(verbose){
      cat("\n----------------------\n SAMPLE SIZE: N = ", N[z], "\n----------------------\n", sep = "")
      message("Simulating datasets")
      pb1 <- txtProgressBar(max = nsims, style = 3)  
    }
    X <- setNames(lapply(seq_len(nsims), function(i){
      x <- do.call(mlGVARsim, A1)
      nx <- setNames(lapply(allNets, function(ii){net(x, ii, d = d)}), allNets)
      if(!is.null(m)){nx$interactions <- netInts(x)}
      if(verbose){setTxtProgressBar(pb1, i)}
      return(list(x = x, nets = list(nets0 = nx)))
    }), nit)
    if(verbose){
      close(pb1)
      message("Fitting models")
      pb2 <- txtProgressBar(max = nsims, style = 3)
    }
    Y <- setNames(lapply(seq_len(nsims), function(i){
      A2$data <- switch(2 - is.null(m), X[[i]]$x$data, X[[i]]$x$mm$data)
      x <- do.call(match.fun(FUN), A2)
      nx0 <- setNames(lapply(allNets, function(ii) net(x, ii)), allNets)
      nx1 <- setNames(lapply(allNets, function(ii) net(x, ii, T)), allNets)
      nx2 <- setNames(lapply(allNets, function(ii) net(x, ii, T, "and")), allNets)
      if(!is.null(m)){
        nx0$interactions <- netInts(x)
        nx1$interactions <- nx2$interactions <- netInts(x, T)
      }
      if(verbose){setTxtProgressBar(pb2, i)}
      return(list(x = x, nets = list(nets0 = nx0, nets1 = nx1, nets2 = nx2)))
    }), nit)
    if(verbose){close(pb2)}
    out <- setNames(lapply(seq_len(nsims), function(n1){
      setNames(lapply(seq_along(Y[[n1]]$nets), function(n2){
        indout1 <- do.call(rbind, lapply(seq_along(Y[[n1]]$nets[[n2]]), function(n3){
          sapply(inds, function(n4){matrixDist(
            Y[[n1]]$nets[[n2]][[n3]], X[[n1]]$nets$nets0[[n3]], 
            n4, isTRUE(n3 %in% c(1, 4)))})
        }))
        indout2 <- do.call(rbind, lapply(seq_along(Y[[n1]]$nets[[n2]]), function(n3){
          performance(Y[[n1]]$nets[[n2]][[n3]], X[[n1]]$nets$nets0[[n3]])
        }))
        netind <- names(Y[[n1]]$nets[[n2]])
        indout <- cbind.data.frame(net = netind, iter = rep(n1, length(netind)), indout1, indout2)
        return(indout)
      }), paste0("nets", 0:2))
    }), nit)
    out0 <- do.call(rbind, lapply(out, '[[', 1))
    out1 <- do.call(rbind, lapply(out, '[[', 2))
    out2 <- do.call(rbind, lapply(out, '[[', 3))
    rownames(out0) <- rownames(out1) <- rownames(out2) <- 1:nrow(out0)
    results <- list(nets0 = out0, nets1 = out1, nets2 = out2)
    if(saveVars){Y1 <- lapply(Y, function(z) z$x$varMods)}
    if(!saveMods){
      X <- lapply(X, '[[', "nets")
      Y <- lapply(Y, '[[', "nets")
    }
    output <- list(sims = X, fits = Y, results = results)
    if(saveVars){output$varMods <- Y1}
    attr(output, "allNets") <- allNets
    return(output)
  }), paste0("N", N))
  outcall <- as.list(match.call())[-1]
  if(any(sapply(outcall, class) %in% c("symbol", "name"))){
    oc <- which(sapply(outcall, class) %in% c("symbol", "name"))
    for(j in seq_along(oc)){
      outcall[[oc[j]]] <- as.character(outcall[[oc[j]]])
      if(outcall[[oc[j]]] %in% c("T", "F")){
        outcall[[oc[j]]] <- as.logical(outcall[[oc[j]]])
      }
    }
  }
  if(any(sapply(outcall, class) == "call")){
    oc <- which(sapply(outcall, class) == "call")
    for(j in seq_along(oc)){outcall[[oc[j]]] <- eval(outcall[[oc[j]]])}
  }
  outcall <- append(outcall, list(
    allNets = switch(2 - is.null(m), allNets, c(allNets, "interactions"))))
  output <- append(list(call = outcall), output)
  attr(output, "time") <- Sys.time() - t1
  if(verbose){cat("\n"); print(Sys.time() - t1)}
  return(output)
}

##### simNet: simulate network structures and data
simNet <- function(N = 100, p = 5, m = FALSE, m2 = .1, b1 = NULL, b2 = NULL, 
                   sparsity = .5, intercepts = NULL, nIter = 250, msym = FALSE, 
                   onlyDat = FALSE, pbar = TRUE, div = 1000, gibbs = TRUE,
                   ordinal = FALSE, nLevels = 5, mord = FALSE, time = TRUE,
                   mbinary = FALSE, minOrd = 3, m1 = NULL, m1_range = NULL,
                   m2_range = c(.1, .3), modType = 'none', lags = NULL, V = 2, 
                   skewErr = FALSE, onlyNets = FALSE, netArgs = NULL,
                   nCores = 1, cluster = 'SOCK', ...){
  t1 <- Sys.time()
  args <- tryCatch({list(...)}, error = function(e){list()})
  if(!is.null(netArgs)){list2env(netArgs[setdiff(names(netArgs), 'data')], envir = environment())}
  if(!is.null(lags) & !identical(as.numeric(lags), 0)){
    args1 <- append(list(nTime = N, nPerson = 1, nNode = p, lag = 1, GGMsparsity = sparsity), as.list(match.call())[-1])
    if(!identical(m, FALSE)){args1 <- replace(args1, 'm', ifelse(mbinary, 'binary', ifelse(mord, 'ordinal', m)))}
    if('nPerson' %in% names(args)){args1 <- replace(args1, 'nPerson', args$nPerson)}
    args1 <- append(args1, args[setdiff(names(args), names(args1))])
    out <- do.call(mlGVARsim, args1[intersect(names(args1), formalArgs('mlGVARsim'))])
    return(out)
  }
  skewErr <- ifelse(is.numeric(skewErr), skewErr, ifelse(!identical(skewErr, FALSE), 3, FALSE))
  gibbsFun <- switch(2 - identical(skewErr, FALSE), function(x){rnorm(1, x, 1)}, 
                     function(x){sn::rsn(1, x, alpha = skewErr)[1]})
  if(!is.null(m1_range)){if(!is.numeric(m1_range) | length(m1_range) != 2){m1_range <- NULL}}
  if(!is.numeric(m2_range) | length(m2_range) != 2){m2_range <- c(.1, .3)}
  if(is.numeric(modType)){modType <- switch(modType + 1, 'none', 'full', 'partial', 'full2', 'partial2')}
  modType <- match.arg(tolower(modType), c('none', 'full', 'partial', 'full2', 'partial2'))
  if(grepl('2', modType) & is.null(m1)){m1 <- .7}
  if(is.null(intercepts)){
    intercepts <- rep(0, p)
  } else if('skewed' %in% intercepts){
    if(length(intercepts) == 1 | length(intercepts) > 2){intercepts <- '3'}
    intercepts <- replicate(p, sn::rsn(1, alpha = as.numeric(setdiff(intercepts, 'skewed'))))
  } else if(length(intercepts) != p){
    intercepts <- rnorm(p)
  }
  mfun <- switch(
    2 - mbinary, 
    function(size = 1){sample(x = 0:1, size = size)}, 
    function(mean = 0){rnorm(n = 1, mean = mean)}
  )
  if(isTRUE(m)){m <- ifelse(mbinary, 1, mfun())}
  if('skewed' %in% m){
    if(length(m) == 1 | length(m) > 2){m <- '3'}
    mskew <- as.numeric(setdiff(m, 'skewed'))
    mfun <- function(mean = 0){sn::rsn(n = 1, xi = mean, alpha = mskew)[1]}
    m <- mfun()
  }
  m <- ifelse(identical(m, FALSE), 0, m)
  if(is.null(b1)){
    b1 <- simPcor(p, sparsity, finalRange = m1_range)
  }
  if(is.null(b2)){
    mat <- b2 <- diag(0, p)
    pn <- (p * (p - 1))/2
    while(all(b2 == 0)){
      if(m2 >= 0 & m2 < 1){
        b2 <- sample(0:1, pn, TRUE, prob = c(1 - m2, m2))
      } else if(m2 >= 1){
        b2 <- numeric(pn)
        b2[sample(1:pn, ifelse(m2 > pn, pn, round(m2)))] <- 1
      }
      b2 <- b2 * sample(c(-1, 1), pn, TRUE, prob = c(.5, .5))
      mat[upper.tri(mat)] <- b2 * runif(pn, min(m2_range), max(m2_range))
      b2 <- as.matrix(Matrix::forceSymmetric(mat))
      if(msym){b2 <- simPcor(p, sparsity = 1 - m2, finalRange = m2_range)}
    }
    if(all(m == 0)){b2 <- diag(0, p)}
  }
  if(modType != 'none' & !all(b2 == 0)){
    while(ifelse(grepl('full', modType), !all(b1[b2 != 0] == 0), any(b1[b2 != 0] == 0))){
      b1 <- simPcor(p, sparsity, finalRange = m1_range)
    }
  }
  diag(b1) <- diag(b2) <- 0
  if(is.null(m1) | identical(m1, 0) | all(m == 0)){
    m1 <- rep(0, p)
  } else {
    if(isTRUE(m1)){m1 <- .7}
    if(length(m1) == 1){
      if(is.null(m1_range)){m1_range <- c(.1, .4)}
      m10 <- runif(p, min(m1_range), max(m1_range))
      if(m1 >= 0 & m1 < 1){
        m10 <- m10 * sample(0:1, p, TRUE, prob = c(1 - m1, m1))
      } else if(m1 >= 1){
        m10[sample(1:p, ifelse(m1 > p, 0, round(p - m1)))] <- 0
      }
      if(grepl('2', modType) & !all(b2 == 0)){
        mm01 <- apply(b2, 1, function(z) any(z != 0))
        while(ifelse(grepl('full', modType), !all(m10[mm01] == 0), any(m10[mm01] == 0))){
          m10 <- runif(p, min(m1_range), max(m1_range))
          if(m1 >= 0 & m1 < 1){
            m10 <- m10 * sample(0:1, p, TRUE, prob = c(1 - m1, m1))
          } else if(m1 >= 1){
            m10[sample(1:p, ifelse(m1 > p, 0, round(p - m1)))] <- 0
          }
        }
      }
      m1 <- m10
    }
  }
  if(onlyNets & !onlyDat){
    return(list(b1 = b1, b2 = b2, intercepts = intercepts, m = m, m1 = m1))
  }
  if(nCores > 1 | isTRUE(nCores)){
    if(isTRUE(nCores)){nCores <- parallel::detectCores()}
    if(grepl('Windows', sessionInfo()$running)){cluster <- 'SOCK'}
    if(tolower(cluster) != 'mclapply'){
      cluster <- match.arg(toupper(cluster), c('SOCK', 'FORK'))
      cl <- parallel::makeCluster(nCores, type = cluster)
    } else {
      cl <- nCores
    }
    sampFun <- function(case, nIter, p, m, b1, b2, intercepts, div, V, m1, skewErr){
      sampling <- matrix(NA, nrow = nIter, ncol = p)
      sampling[1, ] <- rnorm(p)
      m0 <- c(ifelse(m == 0, 0, mfun(m)), numeric(nIter - 1))
      for(iter in 2:nIter){
        m0[iter] <- ifelse(m == 0, 0, mfun(m))
        vv <- switch(V, 1:p, sample(1:p, p, replace = FALSE))
        for(v in vv){
          v_mu <- 0
          v_ps <- which(b1[v, ] != 0)
          if(length(v_ps) > 0){
            for(vp in v_ps){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, sampling[v_it, vp] * b1[v, vp])
            }
          }
          if(m1[v] != 0){v_mu <- c(v_mu, m0[iter] * m1[v])}
          v_ps2 <- which(b2[v, ] != 0)
          if(length(v_ps2) > 0){
            for(vp in v_ps2){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, ((sampling[v_it, vp] * m0[iter]) * b2[v, vp]))
            }
          }
          v_mu <- intercepts[v] + sum(v_mu)
          sampling[iter, v] <- gibbsFun(v_mu)
          if(any(abs(sampling[!is.na(sampling)]) > div)){stop('Sampler diverged')}
        }
      }
      out <- list(data = sampling[nIter, ], M = m0[nIter])
      return(out)
    }
    if(cluster == 'SOCK'){
      parallel::clusterExport(cl, c('mfun', 'sampFun', 'gibbsFun'), envir = environment())
    }
    out <- tryCatch({
      if(pbar){
        pbapply::pboptions(type = 'timer', char = '-')
        pbapply::pblapply(1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1, 
                          b2 = b2, intercept = intercepts, div = div, V = V, 
                          m1 = m1, skewErr = skewErr, cl = cl)
      } else if(tolower(cluster) != 'mclapply'){
        parallel::parLapply(
          cl, 1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1, b2 = b2, 
          intercepts = intercepts, div = div, V = V, m1 = m1, skewErr = skewErr)
      } else {
        parallel::mclapply(
          1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1, b2 = b2, 
          intercepts = intercepts, div = div, V = V, m1 = m1, 
          skewErr = skewErr, mc.cores = nCores)
      }}, 
      error = function(e){TRUE}
    )
    if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
    if(isTRUE(out)){stop('Sampler diverged')}
    data <- do.call(rbind, lapply(out, '[[', 'data'))
    M <- unlist(lapply(out, '[[', 'M'))
    rm(cl)
  } else if(isTRUE(gibbs) | !all(m == 0)){
    data <- matrix(NA, nrow = N, ncol = p); M <- c()
    if(pbar){pb <- txtProgressBar(max = N, style = 3)}
    for(case in 1:N){
      sampling <- matrix(NA, nrow = nIter, ncol = p)
      sampling[1, ] <- rnorm(p)
      m0 <- c(ifelse(m == 0, 0, mfun(m)), numeric(nIter - 1))
      for(iter in 2:nIter){
        m0[iter] <- ifelse(m == 0, 0, mfun(m))
        vv <- switch(V, 1:p, sample(1:p, p, replace = FALSE))
        for(v in vv){
          v_mu <- 0
          v_ps <- which(b1[v, ] != 0)
          if(length(v_ps) > 0){
            for(vp in v_ps){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, sampling[v_it, vp] * b1[v, vp])
            }
          }
          if(m1[v] != 0){v_mu <- c(v_mu, m0[iter] * m1[v])}
          v_ps2 <- which(b2[v, ] != 0)
          if(length(v_ps2) > 0){
            for(vp in v_ps2){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, ((sampling[v_it, vp] * m0[iter]) * b2[v, vp]))
            }
          }
          v_mu <- intercepts[v] + sum(v_mu)
          sampling[iter, v] <- gibbsFun(v_mu)
          if(any(abs(sampling[!is.na(sampling)]) > div)){stop('Sampler diverged')}
        }
      }
      data[case, ] <- sampling[nIter, ]
      M[case] <- m0[nIter]
      if(pbar){
        setTxtProgressBar(pb, case)
        if(case == N){close(pb)}
      }
    }
  } else {
    Sigma <- cov2cor(solve(diag(p) - b1))
    if(identical(skewErr, FALSE)){
      data <- mvtnorm::rmvnorm(n = N, mean = intercepts, sigma = Sigma)
    } else {
      data <- sn::rmsn(n = N, xi = intercepts, Omega = Sigma, alpha = rep(skewErr, p))
    }
  }
  out <- list(b1 = b1, b2 = b2, intercepts = intercepts)
  if(ordinal){
    for(i in 1:ncol(data)){
      ord <- c()
      while(length(unique(ord)) < minOrd){
        ord <- as.numeric(cut(data[, i], sort(c(-Inf, rnorm(nLevels - 1), Inf))))
      }
      data[, i] <- ord 
    }
  }
  if(all(m == 0)){
    out$b2 <- NULL
  } else {
    if(mord & !mbinary){
      if(minOrd < 2){minOrd <- 2}
      while(length(unique(mord)) < minOrd){
        mord <- as.numeric(cut(M, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
      }
      M <- mord
    }
    data <- cbind(data, M)
  }
  out <- append(list(data = data.frame(data)), out)
  if(onlyDat){out <- data.frame(data)}
  if(!all(m == 0)){
    out <- append(out, list(m = m, m1 = m1))
    attributes(out)[c('m2', 'modType')] <- list(m2, modType)
  }
  class(out) <- c(ifelse(onlyDat, 'data.frame', 'list'), 'mgmSim')
  attr(out, 'time') <- t2 <- Sys.time() - t1
  if(time){print(Sys.time() - t1)}
  return(out)
}

##### simNet2: wrapper for simNet to prevent failures
simNet2 <- function(..., nets = NULL, maxiter = 5, pbar = TRUE, time = TRUE){
  args <- args2 <- tryCatch({list(...)}, error = function(e){list()})
  args$onlyNets <- TRUE
  args$pbar <- args2$pbar <- pbar
  args$time <- args2$time <- time
  args2$netArgs <- switch(2 - is.null(nets), do.call(simNet, args), nets)
  tt <- t <- out <- 0
  while(identical(out, 0)){
    out <- tryCatch({do.call(simNet, args2)}, error = function(e){0})
    t <- t + 1
    if(t == maxiter){
      if(!is.null(nets)){stop('Need to generate new network matrices')}
      tt <- tt + 1; t <- 0
      args2$netArgs <- do.call(simNet, args)
      if(tt == maxiter){stop('Parameters may be intractable')}
    }
  }
  return(out)
}

##### summarize: aggregate simulation output
summarize <- function(obj, ind = "mean", ...){
  N <- names(obj)[grepl("N", names(obj))]
  fun <- match.fun(match.arg(ind, c("mean", "median", "sd", "var", "quantile")))
  output <- setNames(lapply(N, function(n) setNames(lapply(obj[[n]]$results, function(z){
    z <- do.call(rbind, lapply(obj$call$allNets, function(i){
      i <- subset(z, net == i)[, -(1:2)]
      apply(i, 2, fun, na.rm = TRUE, ...)
    }))
    rownames(z) <- obj$call$allNets
    return(z)
  }), names(obj[[n]]$results))), N)
  if(length(output) == 1){output <- output[[1]]}
  return(output)
}

##### zeros: return all networks or get the number of zeros for each parameter
zeros <- function(obj, type = "sims", net = "temporal", ind = 2, n = 1, all = FALSE){
  if(is.numeric(type)){type <- c("sims", "fits")[type]}
  if(!isTRUE(tryCatch({match.arg(
    type, c("temporal", "contemporaneous", "between", "interactions"))}, 
    error = function(e){TRUE}))){
    net <- type
    type <- "sims"
  }
  type <- match.arg(type, c("sims", "fits"))
  if(type == "sims"){ind <- 1}
  if(any(grepl("N", names(obj)))){
    allNets <- obj$call$allNets
    if(is.numeric(n)){
      if(n < length(obj)){
        obj <- obj[[n + 1]]
      } else {
        obj <- obj[[grep(paste0("^", n), gsub("N", "", names(obj)))]]
      }
    } else {
      obj <- obj[[grep(n, names(obj))]]
    }
  } else {
    allNets <- attr(obj, "allNets")
  }
  if(is.numeric(net)){net <- allNets[net]}
  net <- match.arg(net, allNets)
  out0 <- do.call(abind::abind, c(lapply(obj[[type]], function(z){
    z[[ind]][[net]]}), along = 3))
  out1 <- t(sapply(seq_len(dim(out0)[1]), function(i){
    sapply(seq_len(dim(out0)[2]), function(j){
      sum(out0[i, j, ] == 0)})}))
  if(all){return(out0)} else {return(out1)}
}

##### sumPlot
sumPlot <- function(obj, ind = "cosine", nets = 1, n = 1, ylim = c(0, 1), ...){
  if(any(grepl("N", names(obj)))){
    if(length(obj) == 2){
      obj <- obj[[2]]
    } else {
      if(is.numeric(n)){
        if(n < length(obj)){
          obj <- obj[[n + 1]]
        } else {
          obj <- obj[[grep(paste0("^", n), gsub("N", "", names(obj)))]]
        }
      } else {
        obj <- obj[[grep(n, names(obj))]]
      }
    }
  }
  obj <- obj$results
  if(is.numeric(nets)){nets <- paste0("nets", 0:2)[nets]}
  nets <- match.arg(nets, paste0("nets", 0:2))
  type <- as.numeric(gsub("nets", "", nets))
  type <- switch(type + 1, "", "(thresholded, OR rule)", "(thresholded, AND rule)")
  obj <- obj[[nets]]
  ind <- match.arg(ind, colnames(obj)[-(1:2)])
  ff <- as.formula(paste0(ind, " ~ net"))
  boxplot(ff, obj, main = paste(ind, type), ylim = ylim, ...)
}

##### simFits: fit a series of mlGVAR models (OLD)
simFits <- function(data, crits, m = NULL, methods = NULL, ...){
  if(any(crits == "all")){
    crits0 <- c("aic", "bic", "ebic", "cv", "cv", "cp", "adjr2")
    crits <- c(crits0, crits[-1])
  }
  crits <- match.arg(tolower(crits), c(
    "aic", "bic", "ebic", "cv", "cp", "adjr2", "r2", "rss"), 
    several.ok = TRUE)
  lam <- rep("min", length(crits))
  seeds <- rep(1, length(crits))
  if(any(crits[duplicated(crits)] == "cv")){
    lam[which(crits == "cv")[2]] <- "1se"
    seeds[which(crits == "cv")] <- 2
  }
  if(is.null(methods)){
    methods <- rep("glmnet", length(crits))
  } else if(length(methods) == 1){
    methods <- rep(methods, length(crits))
  }
  methods[grepl("p|r", crits)] <- "subset"
  x <- length(crits) + 1
  if(!is.null(m)){
    if(!is(m, "list")){m <- rep(list(m), x)}
    if(any(!sapply(m, is.null))){
      methods[which(!sapply(m, is.null)) - 1] <- "glinternet"
    }
  } else {
    m <- rep(list(NULL), x)
  }
  fits0 <- lapply(seq_len(x), function(i){
    if(i == 1){
      pb <- txtProgressBar(min = 0, max = x, style = 3)
      assign("pb", pb, envir = parent.frame(2))
      xx <- mlGVAR(data, m = m[[i]], subjectNets = FALSE, verbose = FALSE, ...)
      setTxtProgressBar(pb, i)
    } else {
      xx <- suppressWarnings(
        mlGVAR(data, m = m[[i]], subjectNets = FALSE, selectFUN = TRUE, 
               verbose = FALSE, criterion = crits[i - 1], which.lam = lam[i - 1], 
               varSeed = switch(seeds[i - 1], NULL, 666), 
               method = methods[i - 1], ...))
      setTxtProgressBar(pb, i)
    }
    return(xx)
  })
  fnames <- paste0("fit", 0:(x - 1))
  names(fits0) <- fnames
  inds <- data.frame(net = fnames, crit = c(0, crits), 
                     lam = c(0, lam), method = c(0, methods))
  if(!is.null(m)){
    m0 <- unlist(lapply(m, function(mm) paste(colnames(data)[mm], collapse = "|")))
    if(any(m0 == "")){m0[m0 == ""] <- 0}
    inds$mods <- m0
  }
  y <- as.data.frame(SURll(fits0))
  y$net <- rownames(y)
  inds0 <- merge(y, inds, by = "net", sort = FALSE)
  fits1 <- fits0[which(!duplicated(inds0[, "LL"]))]
  inds1 <- inds0[which(!duplicated(inds0[, "LL"])), ]
  out <- list(fits = fits1, inds = inds1, fits0 = fits0, inds0 = inds0)
  list2env(out, .GlobalEnv)
  return(message("\nfits, fits0, inds, and inds0 in .GlobalEnv"))
}


### ======================================================================== ###
### ============================ SIM RESULTS =============================== ###
### ======================================================================== ###
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
                        inds = "all", rule = "OR", getVals = FALSE){
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


### ======================================================================== ###
### ======================= EXTRACT CLUSTER RESULTS ======================== ###
### ======================================================================== ###
##### GET RESULTS
getData <- function(fits, data, ind = "cosine", threshold = FALSE, ints = TRUE, 
                    exclude = "PDC", NAs = FALSE, rule = "OR"){
  if(ind %in% c("perf", "performance")){
    mets <- c("sens", "spec", "prec", "acc")
    out <- lapply(seq_along(fits), function(i){
      data.frame(t(performance(fits[[i]], data[[i]], threshold, rule = rule)))})
    out2 <- setNames(lapply(1:4, function(i){
      i <- do.call(rbind, lapply(out, '[[', i))
      colnames(i) <- mets
      return(i)
    }), c("kappa", "beta", "PCC", "PDC"))
    if(!is.null(exclude)){if(all(exclude != "none")){out2 <- out2[-which(names(out2) %in% exclude)]}}
    if(ints){
      out3 <- do.call(rbind, lapply(seq_along(fits), function(i){
        int <- tryCatch({netInts(fits[[i]], threshold = threshold)}, 
                        error = function(e){diag(0, ncol(net(fits[[i]])))})
        return(performance(int, data[[i]]$mm$mb2, threshold, rule = rule))
      }))
      out2$interactions <- as.matrix(out3)
    }
    out <- do.call(rbind, lapply(out2, function(z) data.frame(value = c(z), ind = rep(mets, each = length(fits)))))
    out$net <- factor(rep(names(out2), each = nrow(out)/length(names(out2))))
    out$ind <- factor(out$ind)
    rownames(out) <- 1:nrow(out)
  } else {
    out <- t(sapply(seq_along(fits), function(i){
      compareMods(fits[[i]], data[[i]], threshold = threshold, ind = ind, rule = rule)[[1]]}))
    out <- data.frame(x = c(out), net = rep(c("kappa", "beta", "PCC", "PDC"), each = length(fits)))
    if(!is.null(exclude)){if(all(exclude != "none")){out <- out[-which(out$net %in% exclude), ]}}
    if(ints){
      out2 <- sapply(seq_along(fits), function(i){
        int <- tryCatch({netInts(fits[[i]], threshold = threshold)}, 
                        error = function(e){diag(0, ncol(net(fits[[i]])))})
        return(matrixDist(int, data[[i]]$mm$mb2, ind = ind))
      })
      out <- rbind(out, data.frame(x = out2, net = rep("interactions", length(fits))))
    }
    names(out)[1] <- "value"
  }
  if(!NAs){if(any(is.na(out))){out[is.na(out[, 1]), 1] <- 0}}
  attr(out, "type") <- ind
  return(out)
}

##### GET DESCRIPTIVES
getMeans <- function(x, se = TRUE){
  nn <- as.character(unique(x$net))
  if("ind" %in% names(x)){
    met <- as.character(unique(x$ind))
    out1 <- do.call(rbind, lapply(seq_along(nn), function(i){
      sapply(seq_along(met), function(j){
        mean(subset(x, net == nn[i] & ind == met[j])[, 1], na.rm = TRUE)})}))
    dimnames(out1) <- list(nn, paste0(met, ".M"))
    out2 <- do.call(rbind, lapply(seq_along(nn), function(i){
      sapply(seq_along(met), function(j){
        j <- subset(x, net == nn[i] & ind == met[j])
        if(any(is.na(j[, 1]))){j <- j[!is.na(j[, 1]), ]}
        jj <- sd(j[, 1], na.rm = TRUE)
        if(se){jj <- jj/sqrt(nrow(j))}
        return(jj)
      })}))
    dimnames(out2) <- list(nn, paste0(met, ifelse(se, ".SE", ".SD")))
    out <- cbind(out1, out2)
  } else {
    ms <- sapply(seq_along(nn), function(i) mean(subset(x, net == nn[i])[, 1], na.rm = TRUE))
    ses <- sapply(seq_along(nn), function(i){
      i <- subset(x, net == nn[i])
      if(any(is.na(i[, 1]))){i <- i[!is.na(i[, 1]), ]}
      ii <- sd(i[, 1], na.rm = TRUE)
      if(se){ii <- ii/sqrt(nrow(i))}
      return(ii)
    })
    out <- rbind(ms, ses)
    dimnames(out) <- list(c("mean", ifelse(se, "SE", "SD")), nn)
    if("type" %in% names(attributes(x))){
      rownames(out) <- paste0(attr(x, "type"), c(".M", ifelse(se, ".SE", ".SD")))
    }
  }
  return(out)
}

##### PLOT RESULTS
plotData <- function(x, breaks = NULL){
  require(ggplot2)
  if("ind" %in% names(x)){
    g <- ggplot(x, aes(x = ind, y = value)) + facet_grid(~ net) + xlab("Metric")
  } else {
    g <- ggplot(x, aes(x = net, y = value)) + xlab("Network")
  }
  g <- g + geom_boxplot(outlier.size = .3) + theme_bw()
  if(!is.null(attr(x, "type"))){
    tt <- attr(x, "type")
    if(tt == "cor"){tt <- "Correlation"}
    if(tt == "perf"){tt <- "Performance"}
    if(tt == "cosine"){tt <- "Cosine"}
    if(tt == "mse"){tt <- "MSE"}
    if(tt == "mae"){tt <- "MAE"}
    g <- g + ylab(tt)
  }
  if(!is.null(breaks)){g <- g + scale_y_continuous(breaks = breaks)}
  return(g)
}

##### Extension of 'getData' for multiple folders
makeData <- function(x, ind = 'cosine', exclude = 'PDC', useFits = TRUE, 
                     ints = TRUE, alt = NULL, rule = 'OR'){
  if(!useFits){ints <- FALSE}
  if(!is.null(alt)){
    if(!grepl("^/", alt)){alt <- paste0("/", alt)}
    if(!grepl("[.]RDS$", alt)){alt <- paste0(alt, ".RDS")}
    useFits <- FALSE
  }
  dirs <- dir()[grep(paste0("^", x), dir())]
  stopifnot(length(dirs) > 0)
  dirs <- dirs[order(as.numeric(gsub(paste0(letters, collapse = "|"), "", dirs)))]
  pb <- txtProgressBar(max = length(dirs), style = 3)
  X <- setNames(lapply(seq_along(dirs), function(z){
    intsZ <- ints
    data <- readRDS(paste0(dirs[z], "/data.RDS"))
    if(intsZ & !"mm" %in% unique(lapply(data, names))[[1]]){intsZ <- FALSE}
    fits <- readRDS(paste0(dirs[z], ifelse(
      useFits, "/fits.RDS", ifelse(!is.null(alt), alt, "/vars.RDS"))))
    errs <- which(sapply(fits, length) == 0)
    if(any(errs)){fits <- fits[-errs]; data <- data[-errs]}
    if(length(fits) > 0){
      x1 <- getData(fits, data, ind, exclude = exclude, ints = intsZ, rule = rule)
      x2 <- getData(fits, data, ind, T, exclude = exclude, ints = intsZ, rule = rule)
      x3 <- getData(fits, data, "perf", exclude = exclude, ints = intsZ, rule = rule)
      x4 <- getData(fits, data, "perf", T, exclude = exclude, ints = intsZ, rule = rule)
      out <- list(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
    } else {
      out <- list()
    }
    setTxtProgressBar(pb, z)
    if(z == length(dirs)){close(pb)}
    if(any(errs)){attr(out, 'errs') <- unname(errs)}
    return(out)
  }), dirs)
  if(any(sapply(X, function(z) 'errs' %in% names(attributes(z))))){
    attr(X, 'errs') <- lapply(X, attr, 'errs')
  }
  attr(X, "type") <- ind
  return(X)
}

##### Connected to 'makeData'
makePlots <- function(outs, x = 1, ms = NULL, measure = 'cosine', plot = TRUE, 
                      se = TRUE, critNames = NULL, box = FALSE, vertical = TRUE){
  require(ggplot2)
  measure <- match.arg(tolower(measure), c('cosine', 'correlation', 'mae', 'mse'))
  if(sum(grepl(measure, names(outs))) == 1){outs <- outs[[grep(measure, names(outs))]]}
  if(!is.null(ms)){
    if(isTRUE(ms)){
      ms <- c(.05, .1, .25)
      outs <- appd(outs)
      if(!x %in% 1:2){stop('Not prepared to plot this many factors at once')}
    } else if(is.numeric(ms) & length(ms) == 1){
      outs <- outs[[switch(as.character(ms), '0.05' = 1, '0.1' = 2, '0.25' = 3, ms)]]
      ms <- NULL
    }
  }
  if("type" %in% names(attributes(outs))){outs <- list(outs); one <- TRUE} else {one <- FALSE}
  if("type" %in% names(attributes(outs[[1]]))){measure <- attr(outs[[1]], "type")}
  tt <- as.numeric(gsub(paste0(letters, collapse = "|"), "", names(outs[[1]])))
  Z <- lapply(outs, function(z) lapply(z, function(zz) lapply(zz, getMeans, se = se)))
  Z <- lapply(Z, function(z){
    if(any(sapply(z, length) == 0)){
      kk <- which(sapply(z, length) == 0)
      z <- z[-kk]
    }
    zout <- do.call(rbind, lapply(z, '[[', x))
    attr(zout, 'tt') <- switch(2 - (length(z) == length(tt)), tt, tt[-kk])
    return(zout)
  })
  if(any(grepl('[.]M$', rownames(Z[[1]])))){
    Z <- lapply(Z, function(z){
      tz <- attr(z, 'tt')
      data.frame(tt = rep(tz, ncol(z)), 
                 mean = c(z[seq(1, nrow(z), by = 2), ]),
                 se = c(z[seq(2, nrow(z), by = 2), ]), 
                 Network = rep(Hmisc::capitalize(colnames(z)), each = length(tz)))
    })
    if(is.null(ms) & length(unique(sapply(Z, nrow))) == 1 & (length(unique(gsub("[0-9]", "", c(sapply(outs, names))))) > 1 | !is.null(critNames))){
      dat <- data.frame(do.call(rbind, Z), crit = rep(
        toupper(unique(gsub("[0-9]", "", c(sapply(outs, names))))), 
        each = unique(sapply(Z, nrow))))
      if(!is.null(critNames)){dat$crit <- rep(critNames, each = unique(sapply(Z, nrow)))}
    } else if(!one){
      mm <- sapply(Z, nrow)
      xx <- rep(c('Moderator', 'No Moderator'), length(mm)/2)
      #xx <- c("Moderator", "No Moderator")[order(sapply(Z, nrow), decreasing = TRUE)]
      #mm <- sapply(Z, nrow)[order(sapply(Z, nrow), decreasing = TRUE)]
      #if(length(Z) > 2 & FALSE){
      #  if(length(unique(mm)) == 1){
      #    mm <- rep(2:1, length(Z)/2)
      #  } else {
      #    mm <- sapply(Z, nrow)
      #  }
      #  #stopifnot(all(sapply(Z, nrow)[seq(1, length(Z), by = 2)] > sapply(Z, nrow)[seq(2, length(Z), by = 2)]))
      #  xx[which(mm == max(unique(mm)))] <- "Moderator"
      #  xx[which(mm == min(unique(mm)))] <- "No Moderator"
      #  mm <- sapply(Z, nrow)
      #  crit <- c(unlist(sapply(seq_along(mm), function(z) rep(xx[z], mm[z]))))
      #} else {
      #  crit <- rep(xx, mm)
      #}
      dat <- data.frame(do.call(rbind, Z), crit = rep(xx, mm))
      if(!is.null(ms)){
        Ms0 <- apply(rbind(seq(1, length(mm), by = 2), seq(2, length(mm), by = 2)), 2, function(z) sum(mm[z]))
        dat$M <- rep(ms, Ms0)
      }
    } else if(one){
      dat <- Z[[1]]
    }
    g <- ggplot(dat, aes(x = tt, y = mean, colour = Network, fill = Network)) + theme_bw()
    if(!one & is.null(ms)){if(vertical){g <- g + facet_grid(~ crit)} else {g <- g + facet_grid(crit ~ .)}}
    if(!one & !is.null(ms)){if(vertical){g <- g + facet_grid(M ~ crit)} else {g <- g + facet_grid(crit ~ M)}}
    if(!box){
      g <- g + geom_line() + geom_point(size = 1) + 
        geom_ribbon(aes(ymin = mean - se, ymax = mean + se, colour = NULL), alpha = .3)
    } else { # WORKS, NOT PROPERLY
      g <- g + geom_boxplot(aes(colour = NULL), outlier.size = 1, alpha = .7)
    }
    if(measure %in% c("cor", "correlation")){measure <- "Correlation"}
    if(measure == "cosine"){measure <- "Cosine Similarity"}
    if(measure == "mse"){measure <- "Mean Squared Error (MSE)"}
    if(measure == "mae"){measure <- "Mean Absolute Error (MAE)"}
    g <- g + ylab(measure)
  } else {
    Z <- lapply(Z, function(z){
      tz <- attr(z, 'tt')
      xx <- unique(rownames(z))
      data.frame(tt = rep(rep(tz, each = length(xx)), 4),
                 ind = rep(c("Sensitivity", "Specificity", "Precision", "Accuracy"), each = nrow(z)),
                 mean = c(z[, 1:4]), se = c(z[, 5:8]), Network = rep(Hmisc::capitalize(xx), 4))
    })
    if(length(unique(sapply(Z, nrow))) == 1 & (length(unique(gsub("[0-9]", "", c(sapply(outs, names))))) > 1 | !is.null(critNames))){
      dat <- data.frame(do.call(rbind, Z), crit = rep(
        toupper(unique(gsub("[0-9]", "", c(sapply(outs, names))))), 
        each = unique(sapply(Z, nrow))))
      if(!is.null(critNames)){dat$crit <- rep(critNames, each = unique(sapply(Z, nrow)))}
    } else if(!one){
      mm <- sapply(Z, nrow)
      xx <- rep(c('Moderator', 'No Moderator'), length(mm)/2)
      dat <- data.frame(do.call(rbind, Z), crit = rep(xx, mm))
    } else if(one){
      dat <- Z[[1]]
    }
    g <- ggplot(dat, aes(x = tt, y = mean, colour = Network, fill = Network)) + 
      theme_bw() + ylab("Average Performance")
    if(!one){
      g <- g + facet_grid(ind ~ crit)
    } else {
      if(vertical){g <- g + facet_grid(~ ind)} else {g <- g + facet_grid(ind ~ .)}
    }
    if(!box){
      g <- g + geom_line() + geom_point(size = 1) + 
        geom_ribbon(aes(ymin = mean - se, ymax = mean + se, colour = NULL), alpha = .3)
    } else { # WORKS, NOT PROPERLY
      g <- g + geom_boxplot(aes(colour = NULL), outlier.size = 1, alpha = .7)
    }
  }
  g <- g + xlab("Number of Time Points")
  if(x %in% 1:2){attr(dat, 'measure') <- measure}
  if(plot){return(g)} else {return(dat)}
}

##### appd: take a list and append all sublists; or collect certain lists
appd <- function(x){
  if(length(x) == 1){
    inds <- c('correlation', 'cosine', 'mse', 'mae')
    inds <- ifelse(is.character(x), list(inds), list(paste0(inds, x)))[[1]]
    inds <- paste(paste0('^', inds, '$'), collapse = '|')
    out <- nlist(inds, TRUE)
  } else {
    out <- list()
    for(i in seq_along(x)){out <- append(out, x[[i]])}
  }
  return(out)
}


########## GATHERING DATA
fullData <- function(x, crit = NULL, ints = TRUE, full = FALSE, rmInts = FALSE){
  fullData2 <- function(x, rmInts = FALSE){
    out <- setNames(lapply(1:4, function(z){
      zz <- do.call(rbind, lapply(x, '[[', z))
      zz$net <- factor(zz$net)
      if(is.character(zz$crit)){zz$crit <- factor(zz$crit)}
      if(z %in% 3:4){zz$ind <- factor(zz$ind)}
      return(zz)
    }), paste0('x', 1:4))
    if(all(sapply(1:4, function(z) 'interactions' %in% out[[z]]$net))){
      if(all(sapply(1:4, function(z) 'M' %in% colnames(out[[z]])))){
        if(length(unique(sapply(out, function(z) length(unique(subset(z, net == 'interactions')$M))))) == 1){rmInts <- TRUE}
      }
      if(rmInts){
        out <- lapply(out, function(z){
          zz <- subset(z, net != 'interactions')
          zz$net <- factor(zz$net)
          return(zz)
        })
      }
    }
    return(out)
  }
  if(!is.null(crit)){
    if(all(crit == 'm')){
      crit <- c(.05, .1, .25)
      out <- lapply(1:3, function(z) fullData(x[[z]], crit[z], ints = ints, full = FALSE))
      if(full){out <- fullData2(out, rmInts = rmInts)}
      return(out)
    }
  }
  out <- setNames(lapply(1:4, function(z){
    do.call(rbind, lapply(seq_along(x), function(zz){
      nn <- names(x[[zz]])
      tt <- as.numeric(gsub(paste0(letters, collapse = "|"), "", nn))
      crit <- ifelse(is.null(crit), toupper(unique(gsub("[0-9]", "", nn))), crit)
      out <- lapply(x[[zz]], '[[', z)
      out <- cbind.data.frame(do.call(rbind, out), tt = rep(tt, sapply(out, nrow)), crit = crit)
      if(ints){out$M <- as.numeric(isTRUE("interactions" %in% unique(x[[zz]][[1]][[1]]$net)))}
      rownames(out) <- 1:nrow(out)
      return(out)
    }))
  }), paste0("x", 1:4))
  if(full){out <- fullData2(out, rmInts = rmInts)}
  return(out)
}
mms <- function(x, search = FALSE, double = FALSE, ind = 'net', y = NULL, ...){
  ff <- as.formula(paste0('value ~ .', ifelse(double, '^2', '')))
  if(!is.null(y)){x <- lapply(seq_along(x), function(z) x[[z]][x[[z]][[ind]] == y, -which(colnames(x[[z]]) == ind)])}
  if(!search){
    m <- lapply(seq_along(x), function(z) lm(ff, x[[z]]))
  } else {
    m <- lapply(seq_along(x), function(z) step(lm(ff, x[[z]]), ...))
  }
  return(m)
}
s2pool <- function(x, ind = NULL, vals = NULL, d = TRUE, all = FALSE){
  if(is(x, 'lm')){x <- x$model}
  if(is.null(vals)){vals <- c(0, 1)}
  if(is.null(ind)){ind <- 'M'}
  x1 <- x[x[, ind] == vals[1], 'value']
  x2 <- x[x[, ind] == vals[2], 'value']
  m1 <- mean(x1, na.rm = T)
  m2 <- mean(x2, na.rm = T)
  s1 <- var(x1, na.rm = T)
  s2 <- var(x2, na.rm = T)
  df1 <- length(x1) - 1
  df2 <- length(x2) - 1
  df0 <- df1 + df2
  out <- spool <- ((df1/df0) * s1) + ((df2/df0) * s2)
  if(d){out <- tryCatch({abs(m1 - m2)/sqrt(out)}, error = function(e){NA})}
  if(d){if(m1 == m2){out <- NA} else {names(out) <- paste0("M", 0:1)[which.max(c(m1, m2))]}}
  if(!all){return(out)} else {return(list(out = out, means = c(m1, m2), vars = c(s1, s2), s2pool = spool))}
}
ssx <- function(x, ind = 'net', v = 'PCC'){
  if(all(v == 'all') | isTRUE(v)){v <- unique(x[[1]][, ind])}
  if(is.factor(v)){v <- as.character(v)}
  if(length(v) == 1){
    out <- lapply(x, function(z) z[z[, ind] == v, -which(colnames(z) == ind)])
  } else {
    out <- lapply(seq_along(v), function(z) ssx(x, ind = ind, v = v[z]))
    names(out) <- paste0(ifelse(is.numeric(v), ind, ''), v)
    olen <- unique(sapply(out, length))
    stopifnot(length(olen) == 1)
    stopifnot(length(unique(lapply(out, names))) == 1)
    n2 <- unique(lapply(out, names))[[1]]
    out <- lapply(1:olen, function(z) lapply(out, '[[', z))
    names(out) <- n2
  }
  return(out)
}
s2poolBy <- function(x, ind = 'net'){
  ind <- match.arg(ind, setdiff(unique(unlist(lapply(x, colnames))), 'value'))
  if(ind == 'ind'){x <- x[which(sapply(lapply(x, colnames), function(z) 'ind' %in% z))]}
  xlen <- length(x)
  inds <- unique(x[[1]][, ind])
  if(is.factor(inds)){inds <- as.character(inds)}
  x <- setNames(lapply(inds, function(z){
    ssx(x, ind = ind, v = z)
  }), paste0(ifelse(is.numeric(inds), ind, ''), inds))
  stopifnot(all(unique(sapply(x, length)) == xlen))
  xn <- names(x[[1]])
  x <- setNames(lapply(seq_len(xlen), function(z) lapply(x, '[[', z)), xn)
  out <- setNames(lapply(seq_along(x), function(z){
    z1 <- sapply(x[[z]], s2pool, d = TRUE)
    z2 <- list()
    if(!all(is.na(z1))){
      z2 <- data.frame(do.call(rbind, strsplit(names(z1), '[.]')), d = z1, stringsAsFactors = FALSE)
      z3 <- any(sapply(names(z1), function(k) sum(strsplit(k, '')[[1]] == '.')) > 1)
      if(isTRUE(z3)){z2 <- z2[, -which(apply(z2, 2, function(z) any(grepl(ind, z))))]}
      rownames(z2) <- 1:nrow(z2)
      colnames(z2) <- c(ind, 'M', 'd')
      if(ind == 'crit'){z2$crit <- as.numeric(paste0('.', z2$crit))}
      z2$M <- as.numeric(gsub('M', '', z2$M))
      if(is.numeric(inds)){z2[, ind] <- as.numeric(gsub(paste0(letters, collapse = '|'), '', z2[, ind]))}
    }
    return(z2)
  }), names(x))
  if(any(sapply(out, length) == 0)){out <- out[-which(sapply(out, length) == 0)]}
  return(out)
}
s22 <- function(x, ind1 = 'net', ind2 = 'tt', neg = FALSE){
  inds <- unique(x[[1]][, ind1])
  if(is.factor(inds)){inds <- as.character(inds)}
  X <- lapply(inds, function(z) ssx(x, ind1, z))
  XX <- lapply(X, function(z) s2poolBy(z, ind = ind2))
  poop <- lapply(seq_along(XX), function(z) lapply(XX[[z]], function(zz) cbind(zz, ind = inds[z])))
  poop <- lapply(seq_len(unique(sapply(XX, length))), function(z) lapply(poop, '[[', z))
  poop <- lapply(poop, function(z) do.call(rbind, z))
  names(poop) <- names(XX[[1]])
  if(neg & all(sapply(poop, function(z) 'M' %in% colnames(z)))){
    if(all(sapply(poop, function(z) all(z$M %in% c(0, 1))))){
      poop <- lapply(poop, function(z){
        z[z$M == 0, 'd'] <- -z[z$M == 0, 'd']
        return(z)
      })
    }
  }
  return(poop)
}

### ======================================================================== ###
### ======================================================================== ###
##### nodeMods: Conducts nodewise regression for both MGMs and mVARs
nodeMods <- function(data, y, type, moderators = NULL, lambda = "CV", lags = NULL, 
                     folds = 10, which.lam = "lambda.min", gamma = 0.25, seed = NULL, 
                     threeWay = TRUE, scale = TRUE, measure = "deviance", 
                     alpha = 1, group = NULL, std = TRUE, adaMethod = "ridge", 
                     adaGam = NULL, grPenalty = "gel"){
  yy <- y
  fam <- ifelse(length(type) == 1, ifelse(
    type %in% c("g", "gaussian"), "gaussian", "binomial"), ifelse(
      type[y] == "g", "gaussian", ifelse(type[y] == "c", "multinomial", "binomial")))
  if(length(measure) == 2){measure <- ifelse(type[y] == "g", measure[1], measure[2])}
  data <- setup(data = data, type = type, y = y, lags = lags, scale = scale)
  if(!is.null(lags)){y <- yy <- 1}
  if(is.null(moderators) & "glinternet" %in% lambda){moderators <- 1:ncol(data)}
  X <- interactionMatrix(data = data, y = y, type = type, moderators = moderators, 
                         threeWay = threeWay, lags = lags)
  y <- data[,y]
  if(!is.null(group)){
    if(is.null(moderators)){stop("Can only use group lasso if moderators are included")}
    if(fam == "multinomial" & nrow(data.frame(table(y))) > 2){
      stop("Can't use group lasso on categorical variables with more than 2 levels")
    }
    xMain <- colnames(X)[!grepl(":", colnames(X))]
    xInts <- strsplit(colnames(X)[grepl(":", colnames(X))], ":")
    group <- list()
    for(i in 1:length(xInts)){
      group[[i]] <- c(which(xMain %in% xInts[[i]]), (length(xMain) + i))
      names(group)[i] <- paste0("grp", i)
      stopifnot(paste0(colnames(X)[group[[i]]][1], ":", colnames(X)[group[[i]]][2]) == colnames(X)[group[[i]]][3])
      group[[i]] <- colnames(X)[group[[i]]]
    }
  }
  n <- nrow(X); p <- ncol(X)
  lam <- ifelse(grepl("min", which.lam), "lambda.min", "lambda.1se")
  if("glinternet" %in% lambda){
    X2 <- as.matrix(cbind(y, data[, -yy]))
    if(is.null(lags)){
      m <- ifelse(yy %in% moderators, NA, ifelse(moderators < yy, moderators, moderators - 1))
      if(is.na(m)){m <- NULL}
    } else {
      m <- moderators
    }
    if(all(lambda %in% c("glinternet", "CV"))){
      if(!is.null(seed)){set.seed(seed)}
      fit <- fitHierLASSO(data = X2, type = type, yvar = 1, nfolds = folds, m = m, useSE = TRUE)
      fit[1:2] <- NULL
      lam <- ifelse(lam == "lambda.min", 1, 2)
      betas <- fit$coefs[, lam]
      if(!is.null(m)){
        p1 <- 1:(ncol(data) - 1)
        if(is.null(lags)){m2 <- moderators} else {m2 <- moderators + 1}
        p2 <- which(names(betas) %in% names(betas)[grepl(":", names(betas)) & grepl(paste0(names(data)[m2], collapse = "|"), names(betas))])
        betas <- betas[c(p1, p2)]
      }
      n_neighbors <- sum(betas != 0)
      betas <- as.matrix(c("(Intercept)" = fit$fitobj[[lam + 1]]$betahat[[2]][1], betas), ncol = 1)
      predicted <- cbind(1, X) %*% betas[, 1]
      s2 <- sum((y - predicted)^2)/length(y)
      LL_model <- sum(dnorm(y, mean = predicted, sd = sqrt(s2), log = TRUE))
      deviance <- sum((y - predicted)^2)
      modFitIndex <- -2 * LL_model + n_neighbors * log(n) + 2 * gamma * n_neighbors * log(p)
    } else {
      X2 <- X2[, -1]
      fit <- glinternet(X = X2, Y = y, numLevels = rep(1, ncol(X2)), 
                        interactionCandidates = m, family = fam)
      coefs <- coef(fit)[-1]
      mains <- 1:ncol(X2)
      ints <- t(combn(mains, 2))
      ints2 <- as.numeric(apply(ints, 1, paste, collapse = ""))
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
        names(betas) <- c(colnames(X2), apply(combn(colnames(X2), 2), 2, paste, collapse = ":"))
        if(!is.null(m)){
          p1 <- 1:(ncol(data) - 1)
          if(is.null(lags)){m2 <- moderators} else {m2 <- moderators + 1}
          p2 <- which(names(betas) %in% names(betas)[grepl(":", names(betas)) & grepl(paste0(names(data)[m2], collapse = "|"), names(betas))])
          betas <- betas[c(p1, p2)]
        }
        return(betas)
      })
      n_neighbors <- sapply(allCoefs, function(z) sum(z != 0))
      betas <- lapply(1:length(allCoefs), function(z){
        as.matrix(c("(Intercept)" = fit$betahat[[z + 1]][1], allCoefs[[z]]), ncol = 1)})
      preds <- lapply(1:length(betas), function(z) cbind(1, X) %*% betas[[z]])
      s2 <- sapply(1:length(betas), function(z) sum((y - preds[[z]])^2)/length(y))
      if(fam == "gaussian"){
        LL_models <- sapply(1:length(preds), function(z){
          sum(dnorm(y, mean = preds[[z]], sd = sqrt(s2[z]), log = TRUE))})
        deviance <- sapply(1:length(preds), function(z) sum((y - preds[[z]][, 1])^2))
      } else {
        LL_models <- sapply(1:length(preds), function(z){
          sum(dbinom(y, 1, (exp(preds[[z]])/(1 + exp(preds[[z]]))), log = TRUE))})
        deviance <- sapply(1:length(preds), function(z) -2 * LL_models[z])
      }
      ic_lambda <- -2 * LL_models + n_neighbors * ifelse(
        "AIC" %in% lambda, 2, log(n)) + ifelse(
          "EBIC" %in% lambda, list(2 * gamma * n_neighbors * log(p)), list(0))[[1]]
      n_neighbors <- n_neighbors[which.min(ic_lambda)]
      betas <- betas[[which.min(ic_lambda)]]
      LL_model <- LL_models[which.min(ic_lambda)]
      deviance <- deviance[which.min(ic_lambda)]
      lambda_min <- fit$lambda[which.min(ic_lambda) + 1]
      modFitIndex <- ic_lambda[which.min(ic_lambda)]
      fit <- list(fit = glinternet(
        X = X2, Y = y, numLevels = rep(1, ncol(X2)), 
        interactionCandidates = m, lambda = lambda_min, family = fam))
    }
  }
  if(all(lambda == "CV")){
    if(!is.null(seed)){set.seed(seed)}
    if(is.null(group)){
      if(is.null(adaGam)){
        fit <- cv.glmnet(x = X, y = y, family = fam, type.measure = measure, 
                         nfolds = folds, alpha = alpha, standardize = std)
      } else {
        fit <- adaLasso(x = X, y = y, fam = fam, seed = seed, adaGam = adaGam, 
                        folds = folds, alpha = alpha, adaMethod = adaMethod, 
                        which.lam = lam, measure = measure, 
                        lambda = lambda, std = std)
      }
      betas <- coef(fit, s = lam)
      deviance <- (1 - fit$glmnet.fit$dev.ratio) * fit$glmnet.fit$nulldev
      deviance <- deviance[which(fit$lambda == fit[[lam]])]
    } else {
      if(!is.null(adaGam)){warning("Adaptive lasso not implemented for group lasso")}
      if(fam == "multinomial"){fam <- "binomial"}
      fit <- cv.grpregOverlap(X = X, y = y, family = fam, group = group, 
                              dfmax = length(group) * 3, nfolds = folds, 
                              alpha = alpha, penalty = grPenalty)
      if(fam == "binomial"){
        betas <- coef(fit)/2
        betas <- list(as.matrix((betas * (-1))), as.matrix(betas))
        names(betas) <- c("1", "2")
      } else {
        betas <- as.matrix(coef(fit))
      }
      deviance <- fit$fit$loss[fit$min]
    }
    if(fam == "gaussian"){
      Betas <- matrix(coef(fit, s = lam), ncol = 1)
      predicted <- cbind(1, X) %*% as.vector(Betas)
      s2 <- sum((y - predicted)^2)/length(y)
      LL_model <- sum(dnorm(y, mean = predicted, sd = sqrt(s2), log = TRUE))
      if(is.null(group)){
        n_neighbors <- colSums(matrix(coef(fit, s = lam)[-1, ], ncol = 1) != 0)
      } else {
        n_neighbors <- colSums(matrix(coef(fit)[-1], ncol = 1) != 0)
      }
    }
    if(fam == "multinomial" | fam == "binomial"){
      cats <- unique(y)
      n_cats <- length(cats)
      m_respdum <- matrix(NA, n, n_cats)
      m_coefs <- matrix(NA, n, n_cats)
      m_LL_parts <- matrix(NA, nrow = n, ncol = n_cats + 1)
      for(catIter in 1:n_cats){
        m_respdum[, catIter] <- (y == cats[catIter]) * 1
        if(is.null(group)){
          if(fam == "multinomial"){
            m_coefs[, catIter] <- cbind(rep(1, n), X) %*% matrix(coef(fit, s = lam)[[catIter]], ncol = 1)
          } else {
            m_coefs[, catIter] <- cbind(rep(1, n), X) %*% ((matrix(coef(fit, s = lam), ncol = 1) * ifelse(catIter == 1, -1, 1))/2)
          }
        } else {
          m_coefs[, catIter] <- cbind(rep(1, n), X) %*% betas[[catIter]]
        }
        m_LL_parts[, catIter] <- m_respdum[, catIter] * m_coefs[, catIter]
      }
      m_LL_parts[, (n_cats + 1)] <- -log(rowSums(exp(m_coefs)))
      LL_model <- sum(rowSums(m_LL_parts))
      coefs_bin <- vector("list", length = n_cats)
      if(is.null(group)){
        if(fam == "multinomial"){
          for(ca in 1:n_cats){coefs_bin[[ca]] <- as.matrix(coef(fit, s = lam)[[ca]][-1, ]) != 0}
        } else {
          for(ca in 1:n_cats){coefs_bin[[ca]] <- as.matrix(coef(fit, s = lam)[-1, ]) != 0}
        }
      } else {
        for(ca in 1:n_cats){coefs_bin[[ca]] <- as.matrix(betas[[ca]][-1]) != 0}
      }
      n_neighbors <- colSums(Reduce("+", coefs_bin) != 0)
    }
    modFitIndex <- -2 * LL_model + n_neighbors * log(n) + 2 * gamma * n_neighbors * log(p)
  } else if(!"glinternet" %in% lambda){
    fit <- ModFit(X = X, y = y, fam = fam, criterion = lambda, grPenalty = grPenalty,
                  gamma = gamma, alpha = alpha, group = group, adaGam = adaGam,
                  adaMethod = adaMethod, which.lam = lam, seed = seed, 
                  folds = folds, measure = measure, std = std)
    betas <- fit$model
    names(fit)[which(names(fit) == "lambda")] <- lam
    modFitIndex <- fit[[1]]
    deviance <- fit$deviance
    n_neighbors <- fit$n_neighbors
    LL_model <- fit$LL_mod
  }
  if(is.null(group) & fam == "binomial"){betas <- list(betas)}
  crits <- c("AIC", "BIC", "EBIC")
  final_lam <- ifelse("glinternet" %in% lambda, ifelse(
    any(crits %in% lambda), lambda_min, fit$fitobj$fitCV[[lam + 6]]), fit[[lam]])
  tau <- Tau(betas = betas, x = X, d = ifelse(is.null(moderators), 1, 2))
  mod <- list(ic = modFitIndex, deviance = deviance, LL_model = LL_model,
              n_neighbors = n_neighbors, lambda = final_lam, tau = tau, 
              model = betas, fitobj = fit)
  names(mod)[1] <- ifelse("glinternet" %in% lambda, ifelse(
    any(crits %in% lambda), lambda[lambda %in% crits], "EBIC"), lambda)
  if(all(lambda == "CV")){names(mod)[1] <- "EBIC"}
  mod
}

##### ModFit: Fit single regression model using model selection indices to determine lambda
ModFit <- function(X, y, fam = "gaussian", criterion = "EBIC", gamma = 0.25, 
                   alpha = 1, group = NULL, grPenalty = "cMCP", adaGam = NULL,
                   adaMethod = "ridge", which.lam = "lambda.min", folds = 10,
                   seed = NULL, measure = "deviance", scale = TRUE, std = TRUE){
  criterion <- match.arg(arg = toupper(criterion), choices = c("EBIC", "AIC", "BIC"))
  if(class(y) == "character"){
    if(scale == TRUE){X <- apply(X, 2, scale)}
    X <- data.frame(X)
    X <- model.matrix(~., X)
    Y <- X[,y]
    X <- X[,-which(colnames(X) == y)]
    y <- Y
  }
  if(colnames(X)[1] == "(Intercept)"){X <- X[,-1]}
  if(is.null(group)){
    if(is.null(adaGam)){
      fit <- glmnet(x = X, y = y, family = fam, alpha = alpha, standardize = std)
    } else {
      fit <- adaLasso(x = X, y = y, fam = fam, adaGam = adaGam, folds = folds, 
                      alpha = alpha, adaMethod = adaMethod, which.lam = which.lam, 
                      measure = measure, lambda = criterion, std = std)
    }
  } else {
    if(fam == "multinomial"){fam <- "binomial"}
    fit <- grpregOverlap(X = X, y = y, family = fam, group = group, alpha = alpha,
                         penalty = grPenalty, dfmax = length(group)*3)
    if(fam == "binomial"){
      beta_null <- (coef(fit)[,1])/2
      beta_null <- list(as.matrix((beta_null * (-1))), as.matrix(beta_null))
      names(beta_null) <- c("1", "2")
    }
  }
  n_lambdas <- length(fit$lambda)
  n <- nrow(X); p <- ncol(X)
  if(fam == "gaussian"){
    beta_null <- matrix(coef(fit, s = 1)[1], ncol = 1)
    pred_null <- rep(1, n) * as.vector(beta_null)
    s2 <- sum((y - pred_null)^2)/length(y)
    LL_null <- sum(dnorm(y, mean = pred_null, sd = sqrt(s2), log = TRUE))
  }
  if(fam == "multinomial" | fam == "binomial"){
    cats <- unique(y)
    n_cats <- length(cats)
    m_respdum <- matrix(NA, n, n_cats)
    m_coefs <- matrix(NA, n, n_cats)
    m_LL_parts <- matrix(NA, nrow = n, ncol = n_cats + 1)
    for(catIter in 1:n_cats){
      m_respdum[, catIter] <- (y == cats[catIter]) * 1
      if(is.null(group)){
        if(fam == "multinomial"){
          m_coefs[, catIter] <- cbind(rep(1, n), X) %*% matrix(coef(fit, s = 1)[[catIter]], ncol = 1)
        } else {
          m_coefs[, catIter] <- cbind(rep(1, n), X) %*% ((matrix(coef(fit, s = 1), ncol = 1) * ifelse(catIter == 1, -1, 1))/2)
        }
      } else {
        m_coefs[, catIter] <- cbind(rep(1, n), X) %*% beta_null[[catIter]]
      }
      m_LL_parts[, catIter] <- m_respdum[, catIter] * m_coefs[, catIter]
    }
    m_LL_parts[, n_cats + 1] <- -log(rowSums(exp(m_coefs)))
    LL_null <- sum(rowSums(m_LL_parts))
  }
  if(is.null(group)){
    LL_sat <- 1/2 * fit$nulldev + LL_null
    deviance <- (1 - fit$dev.ratio) * fit$nulldev
  } else {
    LL_sat <- 1/2 * fit$loss[1] + LL_null
    deviance <- fit$loss
  }
  LL_lambda_models <- -1/2 * deviance + LL_sat
  if(fam == "gaussian"){
    if(is.null(group)){
      n_neighbors <- sapply(1:n_lambdas, function(z){
        colSums(as.matrix(coef(fit, s = fit$lambda[z])[-1,]) != 0)})
    } else {
      n_neighbors <- sapply(1:n_lambdas, function(z){
        colSums(as.matrix(coef(fit, lambda = fit$lambda[z])[-1]) != 0)})
    }
  }
  if(fam == "multinomial" | fam == "binomial"){
    n_neighbors <- c()
    if(!is.null(group)){coefs_bin2 <- vector("list", length = n_lambdas)}
    for(NN in 1:n_lambdas){
      coefs_bin <- vector("list", length = n_cats)
      if(is.null(group)){
        for(ca in 1:n_cats){
          if(fam == "multinomial"){
            coefs_bin[[ca]] <- as.matrix(coef(fit, s = fit$lambda[NN])[[ca]][-1, ]) != 0
          } else {
            coefs_bin[[ca]] <- as.matrix(coef(fit, s = fit$lambda[NN])) != 0
          }
        }
      } else {
        betas <- coef(fit, lambda = fit$lambda[NN])/2
        coefs_bin2[[NN]] <- list(as.matrix(betas * (-1)), as.matrix(betas))
        for(ca in 1:n_cats){coefs_bin[[ca]] <- as.matrix(coefs_bin2[[NN]][[ca]][-1]) != 0}
      }
      n_neighbors[NN] <- colSums(Reduce("+", coefs_bin) != 0)
      if(is.null(group) & fam == "binomial"){n_neighbors[NN] <- n_neighbors[NN] - 1}
    }
  }
  ic_lambda <- -2 * LL_lambda_models + n_neighbors * ifelse(
    criterion == "AIC", 2, log(n)) + ifelse(
      criterion == "EBIC", list(2 * gamma * n_neighbors * log(p)), list(0))[[1]]
  lambda_min <- fit$lambda[which.min(ic_lambda)]
  if(is.null(group)){
    lambda_min_model <- coef(fit, s = lambda_min)
  } else {
    if(fam == "binomial"){
      lambda_min_model <- coef(fit, lambda = lambda_min)/2
      lambda_min_model <- list(as.matrix(lambda_min_model * (-1)), as.matrix(lambda_min_model))
      names(lambda_min_model) <- c("1", "2")
    } else {
      lambda_min_model <- as.matrix(coef(fit, lambda = lambda_min))
    }
  }
  output <- list(ic = min(ic_lambda), deviance = deviance[which.min(ic_lambda)],
                 n_neighbors = n_neighbors[which.min(ic_lambda)], 
                 LL_mod = LL_lambda_models[which.min(ic_lambda)],
                 lambda = lambda_min, model = lambda_min_model, fit = fit)
  names(output)[1] <- ifelse(criterion == "AIC", "AIC", ifelse(gamma == 0, "BIC", "EBIC"))
  output
}

##### adaLasso: adaptive lasso
adaLasso <- function(x, y, fam, seed = NULL, adaGam = 1, lambda = "CV", 
                     folds = 10, alpha = 1, adaMethod = "ridge", std = TRUE,
                     which.lam = "lambda.min", measure = "deviance"){
  adaMethod <- match.arg(adaMethod, choices = c("ridge", "lasso", "glm"))
  if(adaMethod == "glm"){
    if(fam == "multinomial"){
      if(nrow(data.frame(table(y))) > 2){
        stop("Can't use glm on categorical variables with more than 2 levels")
      } else {
        fam <- "binomial"
      }
    }
    x <- data.frame(y, x)
    mod <- glm(y ~., x, family = fam)
    betas <- as.matrix(coef(mod))
    w <- 1/abs(matrix(betas[,1][2:(ncol(x[,-1])+1)]))^adaGam
  } else {
    if(adaMethod == "ridge"){alpha_init <- 0}
    if(adaMethod == "lasso"){alpha_init <- 1}
    if(!is.null(seed)){set.seed(seed)}
    initial <- cv.glmnet(x = x, y = y, family = fam, alpha = alpha_init, 
                         nfolds = folds, type.measure = measure, standardize = std)
    if(fam == "multinomial"){
      w <- 1/abs(matrix(coef(initial, s = which.lam)[[1]][,1][2:(ncol(x)+1)]))^adaGam
    } else {
      w <- 1/abs(matrix(coef(initial, s = which.lam)[,1][2:(ncol(x)+1)]))^adaGam
    }
  }
  w[w[,1] == Inf] <- 999999999
  if(lambda == "CV"){
    if(!is.null(seed)){set.seed(seed)}
    fit <- cv.glmnet(x = x, y = y, family = fam, alpha = alpha, nfolds = folds, 
                     type.measure = measure, penalty.factor = w, standardize = std)
  } else {
    fit <- glmnet(x = x, y = y, family = fam, alpha = alpha, penalty.factor = w, standardize = std)
  }
  fit
}

##### Tau: Calculate threshold for nodewise regression models
Tau <- function(betas, x, d = 1){
  p <- ncol(x)
  n <- nrow(x)
  if(class(betas) == "list"){
    B <- list(); tau <- c()
    for(i in 1:length(betas)){
      B[[i]] <- as.numeric(betas[[i]])[-1]
      tau[i] <- sqrt(d) * sqrt(sum(B[[i]]^2)) * sqrt((log(p))/n)
    }
  } else {
    betas <- as.numeric(betas)[-1]
    tau <- sqrt(d) * sqrt(sum(betas^2)) * sqrt((log(p))/n)
  }
  tau
}

##### interactionMatrix: Create design matrix for all models
interactionMatrix <- function(data, y, type, moderators = NULL, 
                              lags = NULL, threeWay = TRUE){
  p <- ncol(data)
  mainMod <- paste(colnames(data)[-y], collapse = " + ")
  if(!is.null(lags)){
    if(!is.null(moderators)){
      p <- ((p - 1)/length(lags)) + 1
      moderators <- moderators + 1
      ints <- t(combn((1:p)[-y], 2))
      ints_m <- c()
      for(i in 1:nrow(ints)){ints_m[i] <- any(ints[i,] %in% moderators)}
      ints <- ints[ints_m,]
      if(length(lags) == 1){
        intMods <- paste0(colnames(data)[ints[,1]], "*", colnames(data)[ints[,2]], collapse = "+")
      } else {
        lagInts <- list(); intMods <- list()
        for(j in 1:length(lags)){
          lagInts[[j]] <- data.frame(y = data[,"y"], data[,grep(paste0("lag", lags[j]), colnames(data))])
          intMods[[j]] <- paste0(colnames(lagInts[[j]])[ints[,1]], "*", colnames(lagInts[[j]])[ints[,2]], collapse = "+")
        }
        for(i in 1:(length(lags) - 1)){intMods[[i]] <- paste0(intMods[[i]], " + ")}
        intMods <- do.call(paste0, intMods)
      }
      form <- as.formula(paste0(colnames(data)[y], " ~ ", mainMod, " + ", intMods))
    } else {
      form <- as.formula(paste0(colnames(data)[y], " ~ ", mainMod))
    }
  } else {
    if(!is.null(moderators) & (if((y %in% moderators) & threeWay == FALSE){length(moderators) > 1} else {TRUE})){
      if((y %in% moderators) & threeWay == TRUE){
        ints <- t(combn((1:p)[-y], 2))
        intMods <- paste0("V", ints[,1], ".*V", ints[,2], ".", collapse = " + ")
      } else {
        if(y %in% moderators){moderators <- moderators[-which(moderators == y)]}
        nm <- length(moderators)
        intMods <- list()
        for(i in 1:nm){
          intMods[[i]] <- paste0(colnames(data)[-c(y, moderators[i])], "*V", moderators[i], ".", collapse = " + ")
        }
        if(nm > 1){for(i in 1:(nm - 1)){intMods[[i]] <- paste0(intMods[[i]], " + ")}}
        intMods <- do.call(paste0, intMods)
      }
      form <- as.formula(paste0(colnames(data)[y], " ~ ", mainMod, " + ", intMods))
    } else {
      form <- as.formula(paste0(colnames(data)[y], " ~ ", mainMod))
    }
  }
  X <- model.matrix(form, data)[,-1]
  X
}


### ======================================================================== ###
### ======================================================================== ###
##### combineMods: Feeds model results to 'combineAll' function; creates n lists for n lags
combineMods <- function(data, mods, type, moderators = NULL, threshold = TRUE,
                        rule = "AND", lags = NULL, scale = TRUE){
  if(length(type) == 1){
    typetype <- type
    type <- rep("c", ncol(data))
    attr(type, "family") <- typetype
  }
  if(is.null(lags)){
    output <- combineAll(data = data, mods = mods, type = type, threshold = threshold,
                         moderators = moderators, rule = rule, scale = scale)
  } else {
    if(length(lags) == 1){
      output <- combineAll(data = data, mods = mods, type = type, threshold = threshold,
                           moderators = moderators, lags = lags, scale = scale)
    } else {
      n_lags <- length(lags)
      lag_mods <- list(); output <- list()
      for(i in 1:n_lags){
        lag_mods[[i]] <- lagMods(mods = mods, type = type, lag = lags[i])
        output[[i]] <- combineAll(data = data, mods = lag_mods[[i]], type = type, 
                                  threshold = threshold, moderators = moderators, 
                                  lags = lags[i], scale = scale)
      }
    }
  }
  output
}

##### combineAll: Main function used to extract and aggregate nodewise coefs for MGMs and mVARs
combineAll <- function(data, mods, type, moderators = NULL, threshold = TRUE,
                       rule = "AND", lags = NULL, scale = TRUE){
  V <- list(); taus <- list()
  if(any(type == "c")){cats <- which(type == "c")}
  for(i in 1:length(mods)){
    taus[[i]] <- ifelse(threshold, mods[[i]]$tau, 0)
    if(type[i] == "g"){
      V[[i]] <- as.numeric(mods[[i]]$model)[-1]
      names(V[[i]]) <- rownames(mods[[i]]$model)[-1]
    } else {
      V[[i]] <- vector("list", length = length(mods[[i]]$model))
      for(j in 1:length(mods[[i]]$model)){
        V[[i]][[j]] <- as.numeric(mods[[i]]$model[[j]])[-1]
        names(V[[i]][[j]]) <- rownames(mods[[i]]$model[[j]])[-1]
      }
    }
  }
  names(V) <- paste0("V", 1:length(V))
  type2 <- sapply(V, class)
  if(any(type2 == "list")){for(i in which(type2 == "list")){V[[i]] <- do.call(rbind, V[[i]])}}
  for(i in which(type2 == "numeric")){V[[i]] <- t(as.matrix(V[[i]]))}
  Vtmp <- list(); combIn <- list(); V1 <- V; Vmods <- list()
  for(i in 1:length(V)){
    for(j in 1:nrow(V[[i]])){V[[i]][j, abs(V[[i]][j,]) < taus[[i]][j]] <- 0}
    V1[[i]] <- V[[i]]
    if(!is.null(moderators)){
      V[[i]] <- V1[[i]][,-grep(":", colnames(V1[[i]]))]
      Vmods[[i]] <- V1[[i]][,grep(":", colnames(V1[[i]]))]
      if(type2[i] == "numeric"){
        V[[i]] <- t(as.matrix(V[[i]]))
        Vmods[[i]] <- t(as.matrix(Vmods[[i]]))
      }
    }
    Vtmp[[i]] <- strsplit(ifelse(is.null(colnames(V[[i]])), list(names(V[[i]])), 
                                 list(colnames(V[[i]])))[[1]], "[.]")
    Vtmp[[i]] <- unlist(Vtmp[[i]])[grep("V", unlist(Vtmp[[i]]))]
    Vtmp[[i]] <- unique(Vtmp[[i]][duplicated(Vtmp[[i]])])
    if(length(Vtmp[[i]]) > 0){
      combIn[[i]] <- vector("list", length = length(Vtmp[[i]]))
      for(m in 1:length(Vtmp[[i]])){
        if(class(V[[i]]) == "numeric"){
          combIn[[i]][[m]] <- abs(V[[i]][grep(Vtmp[[i]][m], names(V[[i]]))])
          V[[i]] <- V[[i]][-grep(Vtmp[[i]][m], names(V[[i]]))[-1]]
        }
        if(class(V[[i]]) == "matrix"){
          combIn[[i]][[m]] <- abs(V[[i]][,grep(Vtmp[[i]][m], colnames(V[[i]]))])
          V[[i]] <- V[[i]][,-grep(Vtmp[[i]][m], colnames(V[[i]]))[-1]]
        }
        if(class(V[[i]]) == "numeric"){
          combIn[[i]][[m]] <- mean(combIn[[i]][[m]])
          V[[i]][grep(Vtmp[[i]][m], names(V[[i]]))[1]] <- combIn[[i]][[m]]
        }
        if(class(V[[i]]) == "matrix"){
          combIn[[i]][[m]] <- mean(rowMeans(combIn[[i]][[m]]))
          V[[i]][,grep(Vtmp[[i]][m], colnames(V[[i]]))[1]] <- combIn[[i]][[m]]
        }
      }
    }
  }
  V2 <- V
  if(any(sapply(V, class) == "matrix")){
    if("c" %in% type){V[which(type == "c")] <- lapply(which(type == "c"), function(z){
      if("family" %in% names(attributes(type))){
        V[[z]][ifelse(attr(type, "family") == "binomial", 1, 2), ]
      } else {
        abs(V[[z]])}
    })}
    V[which(sapply(V, class) == "matrix")] <- lapply(V[which(sapply(V, class) == "matrix")], colMeans)
    names(V) <- names(V2)
  }
  if(is.null(lags)){
    VN <- strsplit(names(unlist(V[which(sapply(V, length) != 0)])), "[.]")
    VN2 <- lapply(1:length(VN), function(z) VN[[z]][1:2])
    matches <- unique(lapply(1:length(VN2), function(x){
      which(sapply(1:length(VN2), function(a) setequal(VN2[[x]], VN2[[a]])) == TRUE)
    }))
    matches1 <- lapply(1:length(matches), function(x) unlist(V)[matches[[x]]])
    for(i in 1:length(matches1)){
      pds <- sum(strsplit(gsub(".$", "", names(matches1[[i]][1])), "")[[1]] == ".")
      names(matches1)[i] <- gsub(paste0(paste(rep(".", pds), collapse = ""), "$"), "", names(matches1[[i]][1]))
    }
    if(any(grepl("[.]$", names(matches1)))){
      names(matches1)[grepl("[.]$", names(matches1))] <- gsub(".$", "", names(matches1)[grepl("[.]$", names(matches1))])
    }
    n <- combn(paste0("V", 1:length(type)), 2)
    n1 <- rbind(as.numeric(gsub("V", "", n[1,])), as.numeric(gsub("V", "", n[2,])))
    n2 <- lapply(1:ncol(n1), function(z) type[n1[,z]])
    choice <- sapply(1:length(n2), function(z) ifelse("c" %in% n2[[z]], "A", "B"))
    if("family" %in% names(attributes(type))){choice <- rep("B", dim(n)[2])}
    for(i in which(choice == "A")){matches1[[i]] <- abs(matches1[[i]])}
    rule <- match.arg(arg = toupper(rule), choices = c("AND", "OR"))
    if(rule == "AND"){
      final <- lapply(1:length(matches1), function(z) mean(matches1[[z]]) * (!(0 %in% matches1[[z]])))
      names(final) <- names(matches1)
    }
    if(rule == "OR"){final <- lapply(matches1, mean)}
    finalCol <- vector("list", length = length(final))
    for(i in 1:length(final)){
      finalCol[[i]] <- ifelse(choice[i] == "A" | final[[i]] == 0, "darkgrey", 
                              ifelse(final[[i]] > 0, "darkgreen", "red"))
    }
    adjMat <- colMat <- matrix(0, ncol = length(V), nrow = length(V))
    edges <- lapply(strsplit(gsub("V", "", names(final)), "[.]"), as.numeric)
    for(i in 1:length(final)){
      adjMat[edges[[i]][2], edges[[i]][1]] <- final[[i]]
      adjMat[edges[[i]][1], edges[[i]][2]] <- final[[i]]
      colMat[edges[[i]][2], edges[[i]][1]] <- finalCol[[i]]
      colMat[edges[[i]][1], edges[[i]][2]] <- finalCol[[i]]
    }
    diag(colMat) <- "darkgrey"
    nwMat <- rbind(c(rep(0, length(V))), do.call(cbind, V))
    diag(nwMat) <- 0
    for(nw in 2:length(V)){nwMat[(1:(nw-1)),nw] <- V[[nw]][1:(nw-1)]}
    rownames(nwMat) <- colnames(nwMat) <- NULL
    nwCol <- c()
    for(i in 1:length(nwMat)){
      nwCol[i] <- ifelse(nwMat[i] == 0, "darkgrey", ifelse(nwMat[i] > 0, "darkgreen", "red"))
    }
    nwCol <- matrix(nwCol, ncol = ncol(nwMat), nrow = nrow(nwMat), byrow = F)
    if(any(type == "c")){
      if(length(which(type == "c")) == 1){
        level <- nrow(data.frame(table(data[,which(type == "c")])))
      } else {
        level <- sapply(apply(data[,which(type == "c")], 2, function(z) data.frame(table(z))), nrow)
      }
      for(i in which(type == "c")){nwCol[,i] <- "darkgrey"}
      if(any(level > 2)){for(i in which(level > 2)){nwCol[i,] <- "darkgrey"}}
    }
    output <- list(adjMat = adjMat, colMat = colMat, nodewise = nwMat, nwCols = nwCol)
  } else {
    V <- do.call(rbind, V)
    if(any(type == "c")){for(i in cats){V[,i] <- abs(V[,i]); V[i,] <- abs(V[i,])}}
    Vnames <- colnames(V)
    if(any(!grepl("[.]$", Vnames))){colnames(V)[!grepl("[.]$", Vnames)] <- gsub(".$", "", colnames(V)[!grepl("[.]$", Vnames)])}
    adjMat <- colMat <- V
    if(any(adjMat == 0)){colMat[adjMat == 0] <- "darkgrey"}
    if(any(adjMat < 0)){colMat[adjMat < 0] <- "red"}
    if(any(adjMat > 0)){colMat[adjMat > 0] <- "darkgreen"}
    if(any(type == "c")){colMat[,cats] <- "darkgrey"; colMat[cats,] <- "darkgrey"}
    output <- list(adjMat = adjMat, colMat = colMat)
  }
  if(!is.null(moderators)){
    if(!is.null(lags)){
      Vlag <- list(); lagNames <- list()
      for(l1 in 1:length(Vmods)){
        Vlag[[l1]] <- strsplit(colnames(Vmods[[l1]]), ":")
        for(l2 in 1:ncol(Vmods[[l1]])){
          for(l3 in 1:2){
            Vlag[[l1]][[l2]][l3] <- gsub(paste0("lag", lags, "."), "", Vlag[[l1]][[l2]][l3])
          }
        }
        Vlag[[l1]] <- unlist(lapply(Vlag[[l1]], paste0, collapse = ":"))
        lagNames[[l1]] <- colnames(Vmods[[l1]])
        colnames(Vmods[[l1]]) <- Vlag[[l1]]
      }
      out3 <- combineLag3way(V = Vmods, cats = ifelse(any(type == "c"), TRUE, FALSE), lags = lags)
    } else {
      if(any(type == "c")){Vmods <- combineMulti3way(data = data, V = Vmods, type = type, scale = scale)}
      if(sum(sapply(Vmods, length)) == 3){
        out3 <- list(unlist(Vmods))
        Vname1 <- strsplit(sapply(Vmods, colnames), ":")
        names(out3) <- paste0(union(Vname1[[1]], Vname1[[2]]), collapse = ":")
      } else {
        out3 <- combine3way(V = Vmods, rule = rule)
      }
    }
    output <- list(output, out3)
  }
  output
}

##### combineLag3way: Get moderator coefs for mVARs; combine within-model coefs for categorical vars
combineLag3way <- function(V, cats, lags){
  if(cats == TRUE){
    Vn <- list()
    for(i in 1:length(V)){
      Vn[[i]] <- strsplit(colnames(V[[i]]), ":")
      for(j in 1:length(Vn[[i]])){
        for(k in 1:2){
          if(!grepl("[.]$", Vn[[i]][[j]][k])){
            Vn[[i]][[j]][k] <- gsub(".$", "", Vn[[i]][[j]][k])
          }
        }
      }
      Vn[[i]] <- unlist(lapply(Vn[[i]], paste0, collapse = ":"))
      colnames(V[[i]]) <- Vn[[i]]
    }
    Vn2 <- lapply(V, colnames)
    if(sum(unlist(lapply(Vn2, duplicated))) != 0){
      Vn3 <- list(); Vn4 <- list()
      for(i in 1:length(Vn2)){
        Vn3[[i]] <- unique(Vn2[[i]][duplicated(Vn2[[i]])])
        Vn4[[i]] <- vector("list", length = length(Vn3[[i]]))
        for(j in 1:length(Vn3[[i]])){
          Vn4[[i]][[j]] <- V[[i]][,grep(Vn3[[i]][j], colnames(V[[i]]))]
          if(class(Vn4[[i]][[j]]) == "numeric"){Vn4[[i]][[j]] <- t(as.matrix(Vn4[[i]][[j]], nrow = 1))}
          Vn4[[i]][[j]] <- as.matrix(rowMeans(abs(Vn4[[i]][[j]])), ncol = 1)
          colnames(Vn4[[i]][[j]]) <- Vn3[[i]][j]
        }
        Vn4[[i]] <- do.call(cbind, Vn4[[i]])
        V[[i]] <- V[[i]][,-grep(paste0(Vn3[[i]], collapse = "|"), colnames(V[[i]]))]
        if(class(V[[i]]) == "numeric"){V[[i]] <- t(as.matrix(V[[i]], nrow = 1))}
        V[[i]] <- cbind(V[[i]], Vn4[[i]])
      }
    }
    for(i in 1:length(V)){
      if(nrow(V[[i]]) != 1){V[[i]] <- t(as.matrix(colMeans(abs(V[[i]])), nrow = 1))}
    }
  }
  V <- lapply(V, function(z) z[,which(z != 0)])
  names(V) <- paste0("V", 1:length(V), ".")
  V <- V[which(sapply(V, length) != 0)]
  V2 <- lapply(V, names)
  V3 <- list()
  if(length(V) != 0){
    for(i in 1:length(V2)){
      V2[[i]] <- strsplit(V2[[i]], ":")
      V3[[i]] <- vector("list", length = length(V2[[i]]))
      for(j in 1:length(V2[[i]])){
        for(k in 1:2){
          V2[[i]][[j]][k] <- gsub("[.]$", paste0(".lag", lags, "."), V2[[i]][[j]][k])
        }
        V3[[i]][[j]] <- paste0(V2[[i]][[j]], collapse = ":")
      }
      V3[[i]] <- unlist(V3[[i]])
      names(V[[i]]) <- V3[[i]]
    }
  }
  V
}

##### combineMulti3way: Aggregate within-model MGM moderator coefs for categorical vars
combineMulti3way <- function(data, V, type, scale = TRUE){
  data <- setup(data = data, type = type, scale = scale)
  if(length(which(type == "c")) == 1){
    level <- nrow(data.frame(table(data[,which(type == "c")])))
  } else {
    level <- sapply(apply(data[,which(type == "c")], 2, function(z) data.frame(table(z))), nrow)
  }
  Vn <- lapply(V, colnames)
  mainNames <- paste0("V", 1:ncol(data), ".")
  if(any(level > 2)){
    mainNames <- mainNames[which(type == "c")[which(level > 2)]]
    Vn2 <- list()
    for(i in 1:length(Vn)){
      if(length(mainNames) > 1){
        Vn2[[i]] <- vector("list", length = length(mainNames))
        for(j in 1:length(mainNames)){Vn2[[i]][[j]] <- Vn[[i]][grep(mainNames[j], Vn[[i]])]}
      } else {
        Vn2[[i]] <- Vn[[i]][grep(mainNames, Vn[[i]])]
      }
    }
    multiLev <- level[which(level > 2)] - 1
    VN2 <- Vn2; VNEW <- list(); MAIN <- mainNames; MULTILEV <- multiLev
    for(M in 1:length(MAIN)){
      if(length(MAIN) > 1){Vn2 <- lapply(VN2, function(z) z[[M]])}
      mainNames <- MAIN[M]
      hasNames <- sapply(Vn2, length)
      multiLev <- MULTILEV[M]
      nameSets <- hasNames/multiLev
      if(any(hasNames > 0)){
        Vn3 <- list(); Vnames <- list(); Vnew <- list()
        for(i in which(hasNames != 0)){
          if(any(nameSets > 1)){
            Vn3[[i]] <- Vnames[[i]] <- vector("list", length = nameSets[i])
            k <- rep(1:nameSets[i], each = multiLev)
            for(j in 1:nameSets[i]){
              Vn3[[i]][[j]] <- Vnames[[i]][[j]] <- Vn2[[i]][which(k == j)]
              Vn3[[i]][[j]] <- matrix(abs(V[[i]][,which(colnames(V[[i]]) %in% Vn3[[i]][[j]])]), ncol = multiLev)
              colnames(Vn3[[i]][[j]]) <- Vnames[[i]][[j]]
            }
            Vn3[[i]] <- lapply(Vn3[[i]], function(z) as.matrix(rowMeans(z), ncol = 1))
            for(j in 1:nameSets[i]){colnames(Vn3[[i]][[j]]) <- Vnames[[i]][[j]][1]}
            Vn3[[i]] <- do.call(cbind, Vn3[[i]])
            Vnew[[i]] <- matrix(V[[i]][,-grep(mainNames, colnames(V[[i]]))], nrow = nrow(V[[i]]))
            colnames(Vnew[[i]]) <- colnames(V[[i]])[-grep(mainNames, colnames(V[[i]]))]
            Vnew[[i]] <- cbind(Vnew[[i]], Vn3[[i]])
          }
          if(any(hasNames == 0)){for(m in which(hasNames == 0)){Vnew[[m]] <- V[[m]]}}
        }
        Vnew <- lapply(Vnew, abs)
        Vnew <- lapply(Vnew, colMeans)
        Vnew <- lapply(Vnew, function(z) t(as.matrix(z, nrow = 1)))
      }
      Vnames3 <- lapply(Vnew, colnames)
      Vnames3.2 <- lapply(Vnames3, strsplit, split = ":")
      for(first in 1:length(Vnames3.2)){
        for(second in 1:length(Vnames3.2[[first]])){
          for(third in 1:2){
            if(!grepl("[.]$", Vnames3.2[[first]][[second]][third])){
              Vnames3.2[[first]][[second]][third] <- gsub(".$", "", Vnames3.2[[first]][[second]][third])
            }
          }
        }
        Vnames3.2[[first]] <- unlist(lapply(Vnames3.2[[first]], function(z) paste0(z, collapse = ":")))
        colnames(Vnew[[first]]) <- Vnames3.2[[first]]
      }
      VNEW[[M]] <- Vnew
    }
  } else {
    Vnames <- ifelse(is.null(unlist(lapply(V, colnames))), list(lapply(V, names)), list(lapply(V, colnames)))[[1]]
    Vnames2 <- lapply(Vnames, strsplit, split = ":")
    for(i in 1:length(Vnames2)){
      for(j in 1:length(Vnames2[[i]])){
        for(k in 1:2){
          if(!grepl("[.]$", Vnames2[[i]][[j]][k])){
            Vnames2[[i]][[j]][k] <- gsub(".$", "", Vnames2[[i]][[j]][k])
          }
        }
      }
      Vnames2[[i]] <- unlist(lapply(Vnames2[[i]], function(z) paste0(z, collapse = ":")))
      if(is.null(colnames(V[[i]]))){names(V[[i]]) <- Vnames2[[i]]} else {colnames(V[[i]]) <- Vnames2[[i]]}
    }
    V <- lapply(1:length(V), function(z){
      if(!is.null(colnames(V[[z]]))){colMeans(abs(V[[z]]))} else {V[[z]]}
    })
    V <- lapply(V, function(z) t(as.matrix(z, nrow = 1)))
    VNEW <- list(V)
  }
  if(length(VNEW) > 1){
    vs <- list(); vs_num <- list(); rmVs <- list()
    for(j in 1:length(VNEW[[1]])){
      vs[[j]] <- unique(colnames(VNEW[[1]][[j]])[duplicated(colnames(VNEW[[1]][[j]]))])
    }
    for(i in which(sapply(vs, length) != 0)){
      vs_num[[i]] <- rep(NA, length(vs[[i]]))
      for(j in 1:length(vs[[i]])){
        vs_num[[i]][j] <- mean(t(matrix(VNEW[[1]][[i]][,which(colnames(VNEW[[1]][[i]]) == vs[[i]][j])])))
        names(vs_num[[i]])[j] <- vs[[i]][j]
      }
      vs_num[[i]] <- t(as.matrix(vs_num[[i]], nrow = 1))
      rmVs[[i]] <- paste0(unlist(vs), collapse = "|")
      VNEW[[1]][[i]] <- t(as.matrix(VNEW[[1]][[i]][,-grep(rmVs[[i]], colnames(VNEW[[1]][[i]]))], nrow = 1))
      VNEW[[1]][[i]] <- cbind(VNEW[[1]][[i]], vs_num[[i]])
    }
  }
  VNEW <- VNEW[[1]]
  VNEW
}

##### combine3way: Aggregate MGM moderator coefs across models
combine3way <- function(V, rule = "AND"){
  Vn <- list()
  for(i in 1:length(V)){
    colnames(V[[i]]) <- paste0("V", i, ".:", colnames(V[[i]]))
    Vn[[i]] <- strsplit(colnames(V[[i]]), split = ":")
  }
  names(V) <- paste0("V", 1:length(V), ".")
  co <- combn(names(V), 3)
  matches3 <- array(NA, dim = c(length(V), max(sapply(V, ncol)), ncol(co)), 
                    dimnames = list(c(1:length(V)), c(1:max(sapply(V, ncol))), c(apply(co, 2, function(z) paste0(z, collapse = ":")))))
  for(i in 1:ncol(co)){
    for(j in 1:length(V)){
      matches3[j,1:ncol(V[[j]]),i] <- unlist(lapply(1:length(Vn[[j]]), function(z) setequal(co[,i], Vn[[j]][[z]])))
    }
  }
  matched <- lapply(1:ncol(co), function(z) any(matches3[,,z] == TRUE))
  matches3 <- matches3[,,which(matched == TRUE)]
  mf1 <- list(); mf2 <- list()
  keep1 <- list(); keep2 <- list()
  for(i in 1:dim(matches3)[3]){
    mf1[[i]] <- matrix(NA, ncol = 2, nrow = 3)
    keep1[[i]] <- which(sapply(1:(dim(matches3)[1]), function(ABC) any(matches3[ABC,,i] == TRUE)))
    keep2[[i]] <- rep(NA, 3)
    for(k in 1:3){keep2[[i]][k] <- which(matches3[keep1[[i]][k],,i] == TRUE)}
    mf1[[i]][,1] <- keep1[[i]]
    mf1[[i]][,2] <- keep2[[i]]
    mf2[[i]] <- c(V[[mf1[[i]][1,1]]][mf1[[i]][1,2]], V[[mf1[[i]][2,1]]][mf1[[i]][2,2]], V[[mf1[[i]][3,1]]][mf1[[i]][3,2]])
  }
  names(mf2) <- attributes(matches3)$dimnames[[3]]
  if(rule == "AND"){
    zeros <- sapply(mf2, function(z) all(z != 0))
    mf2 <- mf2[zeros]
  }
  mf2
}

##### lagIt: Create lagged dataset for MVAR models
lagIt <- function(data, y, lags = 1){
  if(max(lags) > 1){
    lagged <- list()
    for(i in 1:length(lags)){
      lagged[[i]] <- data[-c((nrow(data) - lags[i] + 1):nrow(data)),]
      names(lagged[[i]]) <- paste0(names(lagged[[i]]), "lag", lags[i], ".")
    }
    check <- sapply(lagged, nrow)
    if(any(check != check[length(check)])){
      check2 <- which(check != check[length(check)])
      for(j in 1:length(check2)){
        lagged[[j]] <- lagged[[j]][-c(1:(nrow(lagged[[j]]) - nrow(lagged[[length(lagged)]]))),]
        rownames(lagged[[j]]) <- c(1:nrow(lagged[[j]]))
      }
    }
    lagged <- do.call(cbind, lagged)
  } else {
    lagged <- data[-nrow(data),]
    names(lagged) <- paste0(names(lagged), "lag1.")
  }
  response <- data[-c(1:max(lags)),]
  new <- data.frame(y = response[,y], lagged)
  new
}

##### lagMods: Creates separate models for each lag if lags > 1 is specified
lagMods <- function(mods, type, lag){
  newMods <- list()
  for(i in 1:length(mods)){
    if(type[i] != "g"){
      newMods[[i]] <- vector("list", length = level[i])
      for(j in 1:length(mods[[i]]$model)){
        newMods[[i]][[j]] <- as.matrix(mods[[i]]$model[[j]][grep(paste0("Intercept|", paste0("\\b", "lag", lag, "\\b")), 
                                                                 rownames(mods[[i]]$model[[j]])),], ncol = 1)
      }
    } else {
      newMods[[i]] <- as.matrix(mods[[i]]$model[grep(paste0("Intercept|", paste0("\\b", "lag", lag, "\\b")), 
                                                     rownames(mods[[i]]$model)),], ncol = 1)
    }
    newMods[[i]] <- list(lambda = mods[[i]]$lambda, tau = mods[[i]]$tau, model = newMods[[i]])
  }
  newMods
}

### ======================================================================== ###
### ================= Simulate mlGVAR data (mlVAR Version) ================= ###
### ======================================================================== ###
##### mlGVARsim: main workhorse for simulating VAR and mlGVAR data
mlGVARsim <- function(nTime = 50, nPerson = 10, nNode = 3, m = NULL, lag = 1, 
                      thetaVar = NULL, mu_SD = NULL, init_beta_SD = NULL, 
                      fixedMuSD = 1, manual = FALSE, shrink_fixed = 0.9, 
                      shrink_deviation = 0.9, contemporaneous = "wishart",
                      GGMsparsity = .5, propPos = .5, m0 = 1, m1 = .7, m2 = .25, 
                      m1SD = .1, m1_range = NULL, m2SD = .1, m2_range = NULL, 
                      mcenter = TRUE, getM = FALSE, skew = FALSE, skewErr = FALSE,
                      ordinal = FALSE, nLevels = 5, ordWithin = TRUE, minOrd = 3, 
                      thresholds = NULL, mseed = NULL, keepNets = FALSE){
  if(identical(m, FALSE)){m <- NULL}
  if(minOrd < 2){minOrd <- 2}
  if(is.numeric(ordinal)){nLevels <- ordinal; ordinal <- TRUE}
  if(is.null(thetaVar)){thetaVar <- rep(1, nNode)} ##### START
  if(is.null(mu_SD)){mu_SD <- c(1, 1)}
  if(is.null(init_beta_SD)){init_beta_SD <- c(.1, 1)}
  if(is.null(m1_range)){m1_range <- c(.1, .4)}
  if(is.null(m2_range)){m2_range <- c(.1, .3)}
  if(length(nTime) == 1){nTime <- rep(nTime, nPerson)}
  DF_theta <- nNode * 2
  nTemporal <- nNode^2 * lag
  contemporaneous <- match.arg(contemporaneous, c("wishart", "randomGGM", "fixed"))
  Omega_mu <- simCor(nNode, GGMsparsity)
  Omega_Beta <- simCor(nTemporal, GGMsparsity)
  Omega <- rbind(cbind(Omega_mu, matrix(0, nNode, nTemporal)), 
                 cbind(matrix(0, nTemporal, nNode), Omega_Beta))
  SD <- runif(nNode + nTemporal, c(rep(mu_SD[1], nNode), rep(init_beta_SD[1], nNode)), 
              c(rep(mu_SD[2], nNode), rep(init_beta_SD[2], nNode)))
  Omega <- diag(SD) %*% Omega %*% diag(SD)
  if(contemporaneous == "wishart"){
    Theta_fixed <- simCor(nNode, GGMsparsity)
    Theta_fixed <- diag(sqrt(thetaVar)) %*% Theta_fixed %*% diag(sqrt(thetaVar))
    Theta <- rWishart(nPerson, DF_theta, Theta_fixed/DF_theta)
  } else {
    if(contemporaneous == "randomGGM"){
      Theta <- lapply(1:nPerson, function(x) simCor(nNode, GGMsparsity))
      Theta <- do.call(abind::abind, c(Theta, along = 3))
      Theta_fixed <- apply(Theta, 1:2, mean)
    } else {
      Theta_fixed <- simCor(nNode, GGMsparsity)
      Theta <- lapply(1:nPerson, function(x) Theta_fixed)
      Theta <- do.call(abind::abind, c(Theta, along = 3))
    }
  }
  if(!is.null(m)){
    if(isTRUE(m)){m <- "random"}
    m <- match.arg(tolower(m), c(
      "fixed", "random", "mixed1", "mixed2", "ar", "binary", 
      "skewed", "mlm", "random0", "ordinal"), several.ok = TRUE)
    if(all(c("fixed", "random") %in% m)){m <- m[-which(m == "fixed")]}
    if("ar" %in% m){
      if(m0 >= 1){m0 <- .3}
      mm <- lapply(seq_len(nPerson), function(z){
        as.numeric(arima.sim(n = nTime[z] + 100, model = list(ar = m0)))})
    } else if("binary" %in% m){
      mm <- lapply(seq_len(nPerson), function(z) rbinom(nTime[z] + 100, 1, .5))
    } else if("skewed" %in% m){
      mm <- lapply(seq_len(nPerson), function(z) sn::rsn(nTime[z] + 100, 0, 1, m0))
    } else {
      mm <- lapply(seq_len(nPerson), function(z) rnorm(nTime[z] + 100, 0, m0))
    }
    if(!is.null(mseed)){set.seed(mseed)}
    if("ordinal" %in% m & FALSE){ # Shut down, moved near end of function
      mm <- lapply(mm, function(z){
        ord <- c()
        while(length(unique(ord)) < minOrd){
          ord <- as.numeric(cut(z, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
        }
        return(ord)
      })
    }
    if("mlm" %in% m){ # NOT WORKING YET
      m10 <- rnorm(nNode, 0, mean(m1_range))
      while(any(abs(m10) < min(m1_range)) | any(abs(m10) > max(m1_range))){
        m10 <- rnorm(nNode, 0, mean(m1_range))
      }
      m20 <- rnorm(nNode^2, 0, mean(m2_range))
      while(any(abs(m20) < min(m2_range)) | any(abs(m20) > max(m2_range))){
        m20 <- rnorm(nNode^2, 0, mean(m2_range))
      }
      m10[order(abs(m10))][1:round(nNode * (1 - m1))] <- 0
      m20[order(abs(m20))][1:round((nNode^2) * (1 - m2))] <- 0
      Omega_m1 <- simCor(nNode, (1 - m1))
      Omega_m2 <- simCor(nNode^2, (1 - m2))
      mmb1 <- mvtnorm::rmvnorm(nPerson, m10, Omega_m1)
      mmb2 <- mvtnorm::rmvnorm(nPerson, m20, Omega_m2) # NOPE......
      mmb1 <- lapply(apply(mmb1, 1, list), unlist)
      mmb2 <- lapply(lapply(apply(mmb2, 1, list), unlist), matrix, nNode, nNode)
    } else {
      if(isTRUE(m1) | is.null(m1)){m1 <- .7}
      if(m1 >= 0 & m1 < 1){
        m1 <- sample(c(0, 1), nNode, TRUE, prob = c(1 - m1, m1))
      } else if(m1 >= 1){
        m10 <- numeric(nNode)
        m10[sample(1:nNode, ifelse(m1 > nNode, nNode, round(m1)))] <- 1
        m1 <- m10
      }
      if(m2 >= 0 & m2 < 1){
        m2 <- matrix(sample(c(0, 1), nNode^2, TRUE, prob = c(1 - m2, m2)), nNode, nNode)
      } else if(m2 >= 1){
        m20 <- numeric(nNode^2)
        m20[sample(1:(nNode^2), ifelse(m2 > nNode^2, nNode^2, round(m2)))] <- 1
        m2 <- m20
      }
      m1 <- m1 * sample(c(-1, 1), nNode, TRUE, prob = c(1 - propPos, propPos))
      m2 <- m2 * sample(c(-1, 1), nNode^2, TRUE, prob = c(1 - propPos, propPos))
      if(any(c("fixed", "mixed1") %in% m)){
        m1 <- m1 * runif(nNode, min(m1_range), max(m1_range))
        mmb1 <- rep(list(m1), nPerson)
      } else if("random0" %in% m){
        mmb1 <- lapply(1:nPerson, function(z){
          m10 <- rnorm(nNode, 0, m1SD)
          while(any(abs(m10) < min(m1_range)) | any(abs(m10) > max(m1_range))){
            m10 <- rnorm(nNode, 0, m1SD)
          }
          m10 <- m1 * m10
          return(m10)
        })
      } else {
        mmb1 <- lapply(1:nPerson, function(z) m1 * runif(nNode, min(m1_range), max(m1_range)))
      }
      if(any(c("fixed", "mixed2") %in% m)){
        m2 <- matrix(m2 * runif(nNode^2, min(m2_range), max(m2_range)), ncol = nNode)
        mmb2 <- rep(list(m2), nPerson)
      } else if("random0" %in% m){
        mmb2 <- lapply(1:nPerson, function(z){
          m20 <- rnorm(nNode^2, 0, m2SD)
          while(any(abs(m20) < min(m2_range)) | any(abs(m20) > max(m2_range))){
            m20 <- rnorm(nNode^2, 0, m2SD)
          }
          m20 <- c(m2) * m20
          return(matrix(m20, ncol = nNode))
        })
      } else {
        mmb2 <- lapply(1:nPerson, function(z){
          matrix(c(m2) * runif(nNode^2, min(m2_range), max(m2_range)), ncol = nNode)})
      }
    }
    if(getM){return(list(m = mm[[1]], mb1 = mmb1[[1]], mb2 = mmb2[[1]]))}
  } else {
    mm <- mmb1 <- mmb2 <- NULL
  }
  mu_fixed <- rnorm(nNode, 0, fixedMuSD)
  beta_fixed <- rnorm(nTemporal, 0)
  beta_fixed[order(abs(beta_fixed))[1:round(nTemporal/2)]] <- 0
  mat <- matrix(0, nNode, nNode * lag)
  diag(mat) <- 1
  t1 <- 0
  beta_fixed[c(mat) == 1] <- runif(sum(c(mat) == 1), 0, 1) ##### STOP
  if(lag > 0){
    repeat { ##### START 2
      t1 <- t1 + 1
      if(all(skew == FALSE)){
        Pars <- mvtnorm::rmvnorm(nPerson, c(mu_fixed, beta_fixed), sigma = Omega)
      } else {
        if(is.numeric(skew)){skew <- rep(skew, nNode + nNode^2)}
        if(is.character(skew)){if(all(skew == 'random')){skew <- runif(nNode + nNode^2, -10000, 10000)}}
        if(isTRUE(skew) | length(skew) != (nNode + nNode^2)){skew <- rep(3, nNode + nNode^2)}
        Omega <- as.matrix(forceSymmetric(Omega))
        Pars <- sn::rmsn(nPerson, c(mu_fixed, beta_fixed), Omega, skew)
      }
      Mus <- Pars[, 1:nNode]
      Betas <- array(c(t(Pars[, -(1:nNode)])), c(nNode, nNode * lag, nPerson))
      if(lag > 1){
        under <- cbind(diag(nNode * (lag - 1)), matrix(0, nNode * (lag - 1), nNode))
        ev <- sapply(seq_len(nPerson), function(i){
          mat <- rbind(Betas[, , i], under)
          eigen(mat)$values
        })
      } else {
        ev <- sapply(seq_len(nPerson), function(i){eigen(Betas[, , i])$values})
      }
      allEV <- c(ev)
      if(all(Re(ev)^2 + Im(ev)^2 < 1)){
        if(manual){break}
        if(nPerson == 1){
          Mus <- matrix(mu_fixed, ncol = nNode, nrow = 1)
          Betas <- array(matrix(beta_fixed, nNode, nNode), c(nNode, nNode, 1))
          Theta <- array(Theta_fixed, c(nNode, nNode, 1))
        }
        DataList <- lapply(1:nPerson, function(p){
          parms <- lapply(seq_len(lag), function(l) array(c(Betas[, , p]), c(nNode, nNode, lag))[, , l])
          if(lag > 0){
            if(manual){
              pars("trevSimulateVAR", parms = parms, means = Mus[p, ], 
                   lags = 1, Nt = nTime[p], init = Mus[p, ], burnin = 100, 
                   residuals = Theta[, , p], m = mm[[p]], 
                   mb1 = mmb1[[p]], mb2 = mmb2[[p]])
            }
            if(keepNets){
              res <- list(parms = parms, means = Mus[p, ], lags = seq_len(lag),
                          Nt = nTime[p], init = Mus[p, ], burnin = 100, 
                          residuals = Theta[, , p], m = mm[[p]], mb1 = mmb1[[p]],
                          mb2 = mmb2[[p]], mcenter = mcenter, skewErr = skewErr)
              return(res)
            }
            res <- trevSimulateVAR(parms, means = Mus[p, ], lags = seq_len(lag), 
                                   Nt = nTime[p], init = Mus[p, ], burnin = 100, 
                                   residuals = Theta[, , p], m = mm[[p]], 
                                   mb1 = mmb1[[p]], mb2 = mmb2[[p]], 
                                   mcenter = mcenter, skewErr = skewErr)
          } else {
            res <- mvtnorm::rmvnorm(nTime[p], Mus[p, ], Theta[, , p])
          }
          colnames(res) <- paste0("V", 1:nNode)
          if(ordinal & ordWithin){
            for(vv in 1:ncol(res)){
              tick <- ord <- 0
              while(length(unique(ord)) < minOrd){
                thresh <- switch(2 - is.null(thresholds), rnorm(nLevels - 1), thresholds[[vv]])
                ord <- as.numeric(cut(res[, vv], sort(c(-Inf, thresh, Inf))))
                tick <- tick + 1
                if(tick == 10 & !is.null(thresholds)){thresholds <- NULL}
                if(tick == 20){break}
              }
              res[, vv] <- ord
            }
          }
          res$ID <- p
          res
        })
        if(keepNets){return(DataList[[1]])}
        Data <- do.call(rbind, DataList)
        if(!any(abs(Data[, 1:nNode]) > 100)){break}
      }
      beta_fixed <- beta_fixed * shrink_fixed
      D <- diag(sqrt(diag(Omega)))
      D[-(1:nNode), -(1:nNode)] <- shrink_deviation * D[-(1:nNode), -(1:nNode)]
      Omega <- D %*% cov2cor(Omega) %*% D
    } ##### STOP 2
  } else {
    Pars <- mvtnorm::rmvnorm(nPerson, mu_fixed, sigma = Omega)
    Mus <- Pars[, 1:nNode]
    Betas <- array(dim = c(0, 0, nPerson))
    DataList <- lapply(1:nPerson, function(p){
      res <- as.data.frame(mvtnorm::rmvnorm(nTime[p], Mus[p, ], Theta[, , p]))
      colnames(res) <- paste0("V", 1:nNode)
      res$ID <- p
      res
    })
    Data <- do.call(rbind, DataList)
  }
  model <- list(mu = trevModelArray(mean = mu_fixed, SD = mu_SD, subject = lapply(1:nrow(Mus), function(i) Mus[i, ])), 
                Beta = trevModelArray(
                  mean = array(beta_fixed, c(nNode, nNode, lag)), 
                  SD = array(sqrt(diag(Omega[-(1:nNode), -(1:nNode)])), c(nNode, nNode, lag)), 
                  subject = lapply(1:nPerson, function(p) array(Betas[, , p], c(nNode, nNode, lag)))), 
                Omega_mu = trevModelCov(cov = trevModelArray(mean = Omega[1:nNode, 1:nNode])), 
                Theta = trevModelCov(cov = trevModelArray(mean = Theta_fixed, subject = lapply(1:nPerson, function(p) Theta[, , p]))), 
                Omega = trevModelCov(cov = trevModelArray(mean = Omega)))
  kappa <- corpcor::pseudoinverse(Theta_fixed)
  beta <- matrix(beta_fixed, nNode, nNode)
  mod2 <- list(fixedKappa = kappa, fixedPCC = pcor2(Theta_fixed), fixedBeta = beta, 
               fixedPDC = trevPDC(beta, kappa), between = pcor2(Omega[1:nNode, 1:nNode]))
  if(ordinal & !ordWithin){
    for(vv in 1:(ncol(Data) - 1)){
      tick <- ord <- 0
      while(length(unique(ord)) < minOrd){
        thresh <- switch(2 - is.null(thresholds), rnorm(nLevels - 1), thresholds[[vv]])
        ord <- as.numeric(cut(Data[, vv], sort(c(-Inf, thresh, Inf))))
        tick <- tick + 1
        if(tick == 10 & !is.null(thresholds)){thresholds <- NULL}
        if(tick == 20){break}
      }
      Data[, vv] <- ord
    }
  }
  Results <- list(data = Data, vars = paste0("V", 1:nNode), idvar = "ID", lag = lag, model = model)
  Results <- append(list(data = Data), append(mod2, Results[-1]))
  if(!is.null(m)){
    if('ordinal' %in% m){
      mm <- lapply(mm, function(z){
        ord <- c()
        while(length(unique(ord)) < minOrd){
          ord <- as.numeric(cut(z, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
        }
        return(ord)
      })
    }
    mm <- lapply(mm, function(z) z[-(1:100)])
    dat2 <- data.frame(Data[, 1:nNode], M = unlist(mm), ID = Data[, "ID"])
    mb1 <- if(any(c("fixed", "mixed1") %in% m)){mmb1[[1]]} else {Reduce("+", mmb1)/nPerson}
    mb2 <- if(any(c("fixed", "mixed2") %in% m)){mmb2[[1]]} else {Reduce("+", mmb2)/nPerson}
    Results$mm <- list(data = dat2, m = m, subjects = list(mb1 = mmb1, mb2 = mmb2), mb1 = mb1, mb2 = mb2)
    if(nPerson == 1){Results$mm$data <- Results$mm$data[, -ncol(Results$mm$data)]}
  }
  if(nPerson == 1){
    Results$data <- Results$data[, -ncol(Results$data)]
    Results$between <- NULL
    attr(Results, "simMLgvar") <- TRUE
    class(Results) <- c('list', 'simMLgvar')
  } else {
    class(Results) <- attr(Results, "mlVARsim") <- "mlVARsim"
  }
  return(Results)
}

##### trevSimulateVAR: core single-subject VAR sampler
trevSimulateVAR <- function(parms, means = 0, lags = 1, Nt = 100, init, 
                            residuals = 0.1, burnin, m = NULL, mb1 = NULL, 
                            mb2 = NULL, mcenter = TRUE, skewErr = FALSE){
  if(is.matrix(parms)){parms <- list(parms)} ##### START 3
  if(any(sapply(parms, function(x) length(unique(dim(x))) > 1))){
    stop("non-square graph detected.")
  }
  if(missing(burnin)){burnin <- min(round(Nt/2), 100)}
  Ni <- ncol(parms[[1]])
  if(length(means) == 1){means <- rep(means, Ni)}
  maxLag <- max(lags)
  if(length(residuals) == 1){
    residuals <- diag(residuals, Ni)
  } else if(length(residuals) == Ni){
    residuals <- diag(residuals)
  }
  if(!is.matrix(residuals) && ncol(residuals) != Ni && nrow(residuals) != Ni){
    stop("'residuals' is not a square matrix")
  }
  totTime <- Nt + burnin
  if(missing(init)){init <- matrix(0, maxLag, Ni)}
  Res <- matrix(NA, totTime, Ni)
  if(!is.null(m)){
    if(!is(m, "matrix")){
      if(mcenter){m <- m - mean(m)}
      m <- matrix(rep(m, ncol(Res)), ncol = ncol(Res))
    }
  }
  skewErr <- ifelse(is.numeric(skewErr), skewErr, ifelse(isTRUE(skewErr), 3, FALSE))
  Res[1:maxLag, ] <- init ##### STOP 3
  for(t in (maxLag + 1):(totTime)){
    if(!is.null(m)){
      x1 <- rowSums(do.call(cbind, lapply(seq_along(lags), function(i){
        parms[[i]] %*% (Res[t - lags[i], ] - means)})))
      x2 <- mb1 * m[t - 1, ]
      x3 <- rowSums(mb2 %*% ((Res[t - 1, ] - means) * m[t - 1, ]))
      Res[t, ] <- means + x1 + x2 + x3
    } else {
      Res[t, ] <- means + rowSums(do.call(cbind, lapply(seq_along(lags),function(i){
        parms[[i]] %*% (Res[t - lags[i], ] - means)})))
    }
    e <- switch(
      2 - is.numeric(skewErr), 
      sn::rmsn(1, rep(0, Ni), residuals, rep(skewErr, Ni)), 
      mvtnorm::rmvnorm(1, rep(0, Ni), residuals)
    )
    Res[t, ] <- Res[t, ] + e
  }
  out <- as.data.frame(Res[-(1:burnin), ])
  return(out)
}

##### simPcor: creates network matrices
simPcor <- function(Nvar, sparsity = 0.5, parRange = c(0.5, 1), constant = 1.1, 
                    propPos = 0.5, precision = FALSE, finalRange = NULL){
  trueKappa <- matrix(0, Nvar, Nvar)
  kupper <- upper.tri(trueKappa)
  klower <- lower.tri(trueKappa)
  totEdges <- sum(kupper)
  trueKappa[kupper][sample(seq_len(totEdges), round((1 - sparsity) * totEdges))] <- 1
  vals <- sample(c(-1, 1), totEdges, TRUE, prob = c(propPos, 1 - propPos)) * runif(totEdges, min(parRange), max(parRange))
  trueKappa[kupper] <- trueKappa[kupper] * vals
  trueKappa[klower] <- t(trueKappa)[klower]
  diag(trueKappa) <- constant * rowSums(abs(trueKappa))
  diag(trueKappa) <- ifelse(diag(trueKappa) == 0, 1, diag(trueKappa))
  trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
  trueKappa <- (trueKappa + t(trueKappa))/2
  if(!precision){trueKappa <- as.matrix(qgraph::wi2net(trueKappa))}
  if(!is.null(finalRange)){
    if(max(abs(trueKappa)) > max(finalRange)){
      trueKappa <- trueKappa/(max(abs(abs(trueKappa)))/max(finalRange))
    }
    if(min(abs(trueKappa[trueKappa != 0])) < min(finalRange)){
      wmin <- abs(trueKappa[trueKappa != 0]) < min(finalRange)
      while(min(abs(trueKappa[trueKappa != 0])) < min(finalRange)){
        trueKappa[trueKappa != 0][wmin] <- trueKappa[trueKappa != 0][wmin] * constant
      }
    }
  }
  return(trueKappa)
}

##### simCor: generates correlation matrices based on sparse precision matrices
simCor <- function(x, sparsity = .5, maxiter = 100, ...){
  ow <- getOption("warn"); options(warn = 2)
  try <- 0
  while(try < maxiter){
    out <- tryCatch({
      cov2cor(solve(diag(x) - simPcor(x, sparsity, ...)))}, 
      error = function(e){TRUE})
    if(!isTRUE(out)){break}
    try <- try + 1
  }
  options(warn = ow)
  if(try == maxiter){stop("Estimate failed to converge")}
  return(out)
}

##### trevModelArray: straight from mlVAR
trevModelArray <- function(mean, subject, SE, P, lower, upper, SD){
  Results <- list()
  if(missing(subject)){subject <- NULL}
  if(is.null(subject) && missing(mean)){stop("Either 'subject' or 'mean' must be used")}
  if(missing(mean) && !is.null(subject)){
    N <- length(subject)
    mean <- Reduce("+", subject)/N
  }
  dim <- switch(2 - is.null(dim(mean)), length(mean), dim(mean))
  if(is.null(dim)){dim <- length(mean)}
  if(missing(SD)){
    if(!is.null(subject)){
      N <- length(subject)
      SD <- Reduce("+", lapply(subject, function(x) (x - mean)^2))/(N - 1)
    } else {
      SD <- array(NA, dim)
    }
  }
  if(missing(P)){
    if(!missing(SE)){
      P <- 2 * (1 - pnorm(abs(mean/SE)))
    } else {
      P <- array(NA, dim)
    }
  }
  if(missing(SE)){SE <- array(NA, dim)}
  if(missing(lower)){
    if(!missing(SE)){
      lower <- mean - 1.959964 * SE
    } else {
      lower <- array(NA, dim)
    }
  }
  if(missing(upper)){
    if(!missing(SE)){
      upper <- mean + 1.959964 * SE
    } else {
      upper <- array(NA, dim)
    }
  }
  if(is.null(subject)){subject <- NULL}
  Results[["mean"]] <- mean
  Results[["SD"]] <- SD
  Results[["lower"]] <- lower
  Results[["upper"]] <- upper
  Results[["SE"]] <- SE
  Results[["P"]] <- P
  Results[["subject"]] <- subject
  class(Results) <- c("mlVARarray", "list")
  return(Results)
}

##### trevModelCov: straight from mlVAR
trevModelCov <- function(cov, cor, prec, pcor){
  if(missing(cov)){stop("'cov' can not be missing")}
  if(!missing(cor) && !is(cor, "mlVARarray")){stop("'cor' must be missing or an object of class 'mlVARarray'")}
  if(!missing(prec) && !is(prec, "mlVARarray")){stop("'prec' must be missing or an object of class 'mlVARarray'")}
  if(!missing(pcor) && !is(pcor, "mlVARarray")){stop("'pcor' must be missing or an object of class 'mlVARarray'")}
  cov2corNA <- function(x){
    if(any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(cov2cor(x))
    }
  }
  solveNA <- function(x){
    if(any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(corpcor::pseudoinverse(x))
    }
  }
  cor2pcorNA <- function(x){
    if(any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(corpcor::cor2pcor(x))
    }
  }
  if(missing(cor)){
    cor <- trevModelArray(mean = cov2corNA(cov[["mean"]]), 
                          subject = lapply(cov[["subject"]], cov2corNA))
  }
  if(missing(prec)){
    prec <- trevModelArray(mean = solveNA(cov[["mean"]]), 
                           subject = lapply(cov[["subject"]], solveNA))
  }
  if(missing(pcor)){
    pcor <- trevModelArray(mean = cor2pcorNA(cov[["mean"]]), 
                           subject = lapply(cov[["subject"]], cor2pcorNA))
  }
  Results <- list(cov = cov, cor = cor, prec = prec, pcor = pcor)
  class(Results) <- c("mlVarCov", "list")
  return(Results)
}

##### runSim: not really necessary
runSim <- function(x = 1){
  assign("manual", TRUE, envir = .GlobalEnv)
  ff <- readLines("trevMLVARsim.R")
  f1.1 <- which(ff == "  if(is.null(thetaVar)){thetaVar <- rep(1, nNode)} ##### START")
  f1.2 <- which(ff == "  beta_fixed[c(mat) == 1] <- runif(sum(c(mat) == 1), 0, 1) ##### STOP")
  f2.1 <- which(ff == "    repeat { ##### START 2")
  f2.2 <- which(ff == "    } ##### STOP 2")
  f3.1 <- which(ff == "  if(is.matrix(parms)){parms <- list(parms)} ##### START 3")
  f3.2 <- which(ff == "  Res[1:maxLag, ] <- init ##### STOP 3")
  if(x == 1){
    slines("trevMLVARsim.R", c(f1.1:f1.2, f2.1:f2.2))
  } else if(x == 2){
    slines("trevMLVARsim.R", f3.1:f3.2)
  }
}

##### pcor2: creates partial correlation matrix
pcor2 <- function(x){
  x <- -corpcor::pseudoinverse(x)
  diag(x) <- -diag(x)
  x <- cov2cor(x) - diag(ncol(x))
  return((x + t(x))/2)
}


### ======================================================================== ###
### ======================= graphicalVAR Version =========================== ###
### ======================================================================== ###
##### simMLGVAR: simulate mlGVAR data
simMLGVAR <- function(nTime = 50, nPerson = 20, nVar = 3, propPositive = 0.5, 
                      kappaRange = NULL, betaRange = NULL, betweenRange = NULL, 
                      rewireWithin = 0, betweenVar = 1, withinVar = 0.25, 
                      temporalOffset = 2, m = NULL, m0 = 1, m1 = .3, 
                      m2 = .3, m2pos = .7, m2max = .3, mcenter = TRUE){
  kbb <- list(kappaRange = kappaRange, betaRange = betaRange, betweenRange = betweenRange)
  if(all(sapply(kbb, is.null))){
    kappaRange <- betaRange <- betweenRange <- c(.25, .5)
  } else if(any(sapply(kbb, is.null))){
    kbb[which(sapply(kbb, is.null))] <- c(.25, .5)
    invisible(lapply(1:3, function(i) assign(names(kbb)[i], kbb[[i]], pos = 1)))
  }
  repeat{
    trueKappa <- as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1, nVar, 1, 0)))
    trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] * sample(c(-1, 1), sum(upper.tri(trueKappa)), TRUE, prob = c(propPositive, 1 - propPositive))
    trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]
    trueBeta <- diag(1, nVar)
    for(i in 1:nVar){trueBeta[(i + (temporalOffset - 1))%%nVar + 1, i] <- sample(c(-1, 1), 1, propPositive)}
    trueBetween <- as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1, nVar, 1, 1)))
    trueBetween[upper.tri(trueBetween)] <- trueBetween[upper.tri(trueBetween)] * sample(c(-1, 1), sum(upper.tri(trueBetween)), TRUE, prob = c(propPositive, 1 - propPositive))
    trueBetween[upper.tri(trueBetween)] <- runif(sum(upper.tri(trueBetween)), betweenRange[1], betweenRange[2]) * trueBetween[upper.tri(trueBetween)]
    trueBetween[lower.tri(trueBetween)] <- t(trueBetween)[lower.tri(trueBetween)]
    diag(trueBetween) <- 1
    evK <- round(eigen(trueKappa)$values, 10)
    evB <- round(eigen(trueBeta)$values, 10)
    evBet <- round(eigen(trueBetween)$values, 10)
    if(all(evBet > 0)){break}
  }
  Sigma <- cov2cor(solve(trueBetween))
  D <- diag(sqrt(betweenVar), nVar)
  Means <- mvtnorm::rmvnorm(nPerson, sigma = D %*% Sigma %*% D)
  if(!is.null(m)){
    if(is(m, "list")){
      mm <- m
    } else {
      mm <- list()
      for(j in 1:nPerson){
        if(tolower(m) == "ar"){
          if(m0 >= 1){m0 <- .3}
          mm[[j]] <- as.numeric(arima.sim(n = nTime + 100, list(ar = m0)))
        } else {
          mm[[j]] <- rnorm(nTime + 100, mean = 0, sd = m0)
        }
      }
    }
    if(length(m1) != nVar){m1 <- runif(nVar, -m1, m1)}
    mmb1 <- rep(list(m1), nPerson)
    if(!is(m2, "matrix")){
      m2 <- matrix(sample(c(0, 1), nVar^2, TRUE, prob = c(1 - m2, m2)), nVar, nVar)
      m2 <- m2 * sample(c(-1, 1), nVar^2, TRUE, prob = c(1 - m2pos, m2pos)) * runif(nVar^2, 0, m2max)
    }
    mmb2 <- rep(list(m2), nPerson)
  } else {
    mm <- mmb1 <- mmb2 <- NULL
  }
  SubjectData <- lapply(1:nPerson, function(i){
    try <- 1
    maxtry <- 10
    repeat{
      kappa <- trueKappa
      kappa[upper.tri(kappa)] <- runif(sum(upper.tri(kappa)), kappaRange[1], kappaRange[2]) * kappa[upper.tri(kappa)]
      kappa[lower.tri(kappa)] <- t(kappa)[lower.tri(kappa)]
      diag(kappa) <- 1
      beta <- trueBeta * runif(nVar^2, betaRange[1], betaRange[2])
      kappa <- trevRewire(kappa, rewireWithin)
      beta <- trevRewire(beta, rewireWithin)
      evK <- eigen(kappa)$values
      evB <- eigen(beta)$values
      while(any(Re(evB)^2 + Im(evB)^2 > 1)){
        warning("Shrinking parameters")
        beta <- 0.95 * beta
        evB <- eigen(beta)$values
      }
      if(all(evK > 0) & all(Re(evB)^2 + Im(evB)^2 < 1)){break}
      try <- try + 1
      if(try > maxtry){stop("Maximum number of tries reached.")}
    }
    D <- diag(sqrt(withinVar), nVar)
    Delta <- diag(1/sqrt(diag(solve(kappa))))
    kappa <- solve(D) %*% solve(Delta) %*% kappa %*% solve(Delta) %*% solve(D)
    Data <- as.data.frame(
      trevSimGVAR(nTime, beta, kappa, mean = Means[i, ], 
                  m = mm[[i]], mb1 = mmb1[[i]], mb2 = mmb2[[i]], 
                  mcenter = mcenter))
    Data$ID <- i
    return(list(kappa = kappa, beta = beta, PCC = trevPCC(kappa), 
                PDC = trevPDC(beta, kappa), data = Data))
  })
  fixedKappa <- Reduce("+", lapply(SubjectData, "[[", "kappa"))/nPerson
  fixedBeta <- Reduce("+", lapply(SubjectData, "[[", "beta"))/nPerson
  allData <- do.call(rbind, lapply(SubjectData, "[[", "data"))
  Results <- list(data = allData, fixedKappa = fixedKappa, 
                  fixedPCC = trevPCC(fixedKappa), fixedBeta = fixedBeta, 
                  fixedPDC = trevPDC(fixedBeta, fixedKappa), between = Sigma,
                  Omega = trueBetween, means = Means, personData = SubjectData, 
                  idvar = "ID", vars = names(allData)[names(allData) != "ID"])
  if(!is.null(m)){
    mm <- lapply(mm, function(z) z[-(1:100)])
    dat2 <- data.frame(allData[, 1:nVar], M = unlist(mm), ID = allData[, "ID"])
    Results$mm <- list(data = dat2, mb1 = mmb1[[1]], mb2 = mmb2[[1]])
  }
  if(nPerson == 1){
    Results[c("personData", "idvar")] <- NULL
    Results$data <- Results$data[, -which(colnames(Results$data) == "ID")]
    if(!is.null(m)){Results$mm$data <- Results$mm$data[, -ncol(Results$mm$data)]}
  }
  attr(Results, "simMLgvar") <- TRUE
  return(Results)
}

##### trevSimGVAR
trevSimGVAR <- function(nTime = 100, beta = NULL, kappa = NULL, mean = 0, 
                        init = NULL, warmup = 100, lbound = -Inf, ubound = Inf, 
                        m = NULL, mb1 = NULL, mb2 = NULL, mcenter = TRUE, ...){
  ubound <- Inf; lbound <- -Inf; warmup <- 100; init <- NULL; mm <- NULL
  if(length(mean) == 1){mean <- rep(mean, ncol(kappa))}
  if(length(lbound) == 1){lbound <- rep(lbound, ncol(kappa))}
  if(length(ubound) == 1){ubound <- rep(ubound, ncol(kappa))}
  if(is.null(init)){init <- mean}
  stopifnot(!is.null(beta))
  stopifnot(!is.null(kappa))
  Nvar <- ncol(kappa)
  init <- rep(init, length = Nvar)
  totTime <- nTime + warmup
  Data <- t(matrix(init, Nvar, totTime))
  if(!is.null(m)){
    if(!is.numeric(m)){
      mm <- mlGVARsim(nTime = nTime, nPerson = 1, nNode = Nvar, 
                      m = m, mcenter = mcenter, getM = TRUE, ...)
      m <- mm$m
      mb1 <- mm$mb1
      mb2 <- mm$mb2
    }
    if(!is(m, "matrix")){
      if(mcenter){m <- m - mean(m)}
      m <- matrix(rep(m, ncol(Data)), ncol = ncol(Data))
    }
  }
  Sigma <- solve(kappa)
  for(t in 2:totTime){
    if(!is.null(m)){
      x1 <- t(beta %*% (Data[t - 1, ] - mean))
      x2 <- mb1 * m[t - 1, ]
      x3 <- t(mb2 %*% ((Data[t - 1, ] - mean) * m[t - 1, ]))
      x4 <- mvtnorm::rmvnorm(1, rep(0, Nvar), Sigma)
      Data[t, ] <- mean + x1 + x2 + x3 + x4
    } else {
      Data[t, ] <- mean + t(beta %*% (Data[t - 1, ] - mean)) + mvtnorm::rmvnorm(1, rep(0, Nvar), Sigma)
    }
    Data[t, ] <- ifelse(Data[t, ] < lbound, lbound, Data[t, ])
    Data[t, ] <- ifelse(Data[t, ] > ubound, ubound, Data[t, ])
  }
  if(!is.null(mm)){Data <- cbind(Data, mm$m)}
  output <- Data[-seq_len(warmup), , drop = FALSE]
  if(!is.null(mm)){output <- list(data = output, mb1 = mb1, mb2 = mb2)}
  return(output)
}

##### trevRewire
trevRewire <- function(x, p, directed){
  if(missing(directed)){directed <- !all(x == t(x))}
  ind <- if(directed){diag(1, ncol(x)) != 1} else {upper.tri(x)}
  curEdges <- which(x != 0 & ind, arr.ind = TRUE)
  toRewire <- which(runif(nrow(curEdges)) < p)
  if(any(x == 0)){
    for(i in seq_along(toRewire)){
      curZeros <- which(x == 0 & ind, arr.ind = TRUE)
      dest <- sample(seq_len(nrow(curZeros)), 1)
      x[curZeros[dest, 1], curZeros[dest, 2]] <- x[curEdges[toRewire[i], 1], curEdges[toRewire[i], 2]]
      x[curEdges[toRewire[i], 1], curEdges[toRewire[i], 2]] <- 0
      if(!directed){
        x[curZeros[dest, 2], curZeros[dest, 1]] <- x[curEdges[toRewire[i], 2], curEdges[toRewire[i], 1]]
        x[curEdges[toRewire[i], 2], curEdges[toRewire[i], 1]] <- 0
      }
    }
  }
  return(x)
}

##### trevPDC
trevPDC <- function(beta, kappa){
  if(ncol(beta) == nrow(beta) + 1){beta <- beta[, -1, drop = FALSE]}
  sigma <- solve(kappa)
  t(beta/sqrt(diag(sigma) %o% diag(kappa) + beta^2))
}

##### trevPCC
trevPCC <- function(x){
  x <- -cov2cor(x)
  diag(x) <- 0
  return(as.matrix(Matrix::forceSymmetric(x)))
}

##### randomGVAR: simulate the matrices only
randomGVAR <- function(Nvar, probKappaEdge = 0.1, probBetaEdge = 0.1, 
                       probKappaPositive = 0.5, probBetaPositive = 0.5, 
                       maxtry = 10, kappaConstant = 1.1){
  try <- 0
  repeat {
    kappaRange = c(0.5, 1)
    trueKappa <- matrix(0, Nvar, Nvar)
    trueKappa[upper.tri(trueKappa)] <- sample(c(0, 1), sum(upper.tri(trueKappa)), TRUE, prob = c(1 - probKappaEdge, probKappaEdge))
    trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] * sample(
      c(1, -1), sum(upper.tri(trueKappa)), TRUE, prob = c(1 - probKappaPositive, probKappaPositive)) * runif(
        sum(upper.tri(trueKappa)), min(kappaRange), max(kappaRange))
    trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]
    diag(trueKappa) <- kappaConstant * rowSums(abs(trueKappa))
    diag(trueKappa) <- ifelse(diag(trueKappa) == 0, 1, diag(trueKappa))
    trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
    trueKappa <- (trueKappa + t(trueKappa))/2
    Vmin <- min(abs(trueKappa[trueKappa != 0]))
    trueBeta <- matrix(sample(c(0, 1), Nvar^2, TRUE, prob = c(1 - probBetaEdge, probBetaEdge)), Nvar, Nvar)
    trueBeta <- trueBeta * sample(c(-1, 1), Nvar^2, TRUE, prob = c(1 - probBetaPositive, probBetaPositive)) * runif(Nvar^2, Vmin, 1)
    diag(trueBeta) <- Vmin
    evK <- eigen(trueKappa)$values
    evB <- eigen(trueBeta)$values
    while(any(Re(evB)^2 + Im(evB)^2 > 1)){
      trueBeta <- 0.95 * trueBeta
      evB <- eigen(trueBeta)$values
    }
    if(all(evK > 0) & all(Re(evB)^2 + Im(evB)^2 < 1)){break}
    try <- try + 1
    if(try > maxtry){stop("Maximum number of tries reached.")}
  }
  Res <- list(kappa = trueKappa, beta = trueBeta, PCC = trevPCC(trueKappa), 
              PDC = trevPDC(trueBeta, trueKappa))
  class(Res) <- "gVARmodel"
  return(Res)
}

##### cor0: deal with sd == 0
cor0 <- function(x, y){
  if(all(is.na(x)) | all(is.na(y))){return(NA)}
  if(sd(x, na.rm = TRUE) == 0 | sd(y, na.rm = TRUE) == 0){return(0)}
  return(cor(x, y, use = 'pairwise.complete.obs'))
}


### ======================================================================== ###
### ============================ CROSS-SECTIONAL =========================== ###
### ======================================================================== ###
##### mgmSim: simulates data with interactions from mgm
mgmSim <- function(N = 500, p = 5, prop = .2, bmin = .1, bmax = .2, exo = NULL,
                   propPos = .5, omega_mu = 1, muSD = 1, maxiter = 10, mgm2 = FALSE){
  if(mgm2){source('editing/extras/mgmsampler2.R')}
  ff <- t(combn(p, 2))
  n <- nrow(ff)
  k <- floor(n * prop)
  edges <- ff[sample(1:n, k, replace = FALSE), ]
  if(!is.null(exo) & mgm2){
    while(exo %in% edges){
      edges <- ff[sample(1:n, k, replace = FALSE), ]
    }
  }
  factors <- list(edges, matrix(c(2, 3, 4), ncol = 3, byrow = T))
  ints0 <- runif(k, bmin, bmax) * sample(
    c(-1, 1), k, TRUE, prob = c(propPos, 1 - propPos))
  interactions <- list(lapply(as.list(ints0), array, dim = c(1, 1)), 
                       list(array(.15, dim = c(1, 1, 1))))
  intercepts <- as.list(rnorm(p, 0, muSD))
  sds <- if(length(muSD) == 1){rep(muSD, p)} else {muSD[1:k]}
  if(any(is.na(sds))){sds[is.na(sds)] <- 1}
  oldw <- getOption("warn"); options(warn = 2)
  out <- iter <- TRUE
  while(isTRUE(out)){
    out <- tryCatch({
      if(mgm2){
        mgmsampler2(factors = factors, interactions = interactions, 
                    thresholds = intercepts, sds = sds, type = rep("g", p), 
                    level = rep(1, p), N = N, nIter = 100, pbar = TRUE,
                    exogenous = exo)
      } else {
        mgm::mgmsampler(factors = factors, interactions = interactions, 
                        thresholds = intercepts, sds = sds, type = rep("g", p), 
                        level = rep(1, p), N = N, nIter = 100, pbar = TRUE)
      }
    }, error = function(e){TRUE})
    if(isTRUE(out)){
      interactions[[2]] <- lapply(interactions[[2]], function(z){
        z[, , 1] <- z[, , 1] * .9
        return(z)
      })
    }
    iter <- as.numeric(iter) + 1
    if(iter > maxiter){break}
  }
  out$chains <- NULL
  out$trueNet <- diag(0, p)
  for(i in seq_len(k)){out$trueNet[edges[i, 1], edges[i, 2]] <- ints0[i]}
  out$trueNet <- out$trueNet + t(out$trueNet)
  class(out) <- c('list', 'mgmSim')
  options(warn = oldw)
  return(out)
}

##### mgm2: just quick 'mgm' wrapper
mgm2 <- function(data, m = NULL, sel = 'EBIC', gamma = .25, ...){
  if(identical(m, FALSE)){m <- NULL}
  type <- rep('g', ncol(data))
  level <- rep(1, ncol(data))
  mgm::mgm(data = data, type = type, level = level, lambdaSel = sel, 
           lambdaGam = gamma, moderators = m, ...)
}


### ======================================================================== ###
### ======================================================================== ###
setupModnets <- function(folders = '.', rm = FALSE){
  if(any(folders %in% c('.', ''))){folders <- '.'}
  fold <- "~/C:\\ Desktop/COMPS/METHODS/CODE/modnets/"
  mods <- paste0(c("functions", "ggm", "testing", "centrality", "sim", 
                   "mlGVAR", "simGVAR", "penalized", "power", "plots"), ".R")
  for(j in seq_along(folders)){
    if(!folders[j] %in% dir() & folders[j] != '.'){
      stopifnot(length(folders) == 1)
      folders[j] <- '.'
    }
    if(isTRUE(rm)){system(paste0('rm -rf ', folders[j], '/modnets'))}
    system(paste0('mkdir ', folders[j], '/modnets'))
    for(k in seq_along(mods)){system(paste0("cp ", fold, mods[k], " ./", folders[j], "/modnets/"))}
    f0 <- readLines(paste0(folders[j], "/modnets/functions.R"))
    f1 <- which(f0 == "setwd('~/C: Desktop/COMPS/METHODS/CODE/modnets')")
    f2 <- which(f0 == "files <- paste0('./', c('ggm', 'testing', 'centrality', 'sim', 'mlGVAR', 'simGVAR', 'penalized', 'power', 'plots'),'.R')")
    f0[f1] <- paste0("#", f0[f1])
    f0[f2] <- gsub("'./'", "'modnets/'", f0[f2])
    f0 <- paste0(paste0(f0, collapse = "\n"), "\n")
    cat(f0, file = paste0(folders[j], "/modnets/functions.R"))
  }
}

### ======================================================================== ###
### ========================= Simulate MLVAR data ========================== ###
### ======================================================================== ###
##### VARsim: simulate VAR data for single subject
VARsim <- function(N, p, b1, b2 = NULL, c1 = NULL, c2 = NULL, covariate = NULL, 
                   b0 = 0, eps = 0.2, burnin = 100, getDat = TRUE, 
                   fullDat = FALSE, tvInit = FALSE){
  N <- N + burnin
  y <- matrix(0, N, p)
  if(!is(eps, "matrix")){
    sigma2 <- matrix(eps, p, p)
    diag(sigma2) <- 1
  } else {
    sigma2 <- eps
  }
  if(all(b0 == 0)){b0 <- sapply(rep(0, p), rep, N)}
  else if(length(b0) == 1){b0 <- sapply(rep(b0, p), rep, N)}
  else if(length(b0) == p){b0 <- sapply(b0, rep, N)}
  else if(length(b0) == N){b0 <- matrix(b0, N, p)}
  if(is.null(b2)){b2 <- rep(0, p)} else {if(length(b2) != p * (p - 1)){stop("Length of b2 needs to equal p * (p - 1)")}}
  if(is.null(c1)){c1 <- rep(0, p)} else {if(length(c1) != p){stop("Length of c1 needs to be equal to p")}}
  if(is.null(c2)){c2 <- rep(0, p)} else {if(length(c2) != p * p){stop("Length of c2 needs to be equal to p^2")}}
  if(is.null(covariate)){covariate <- rep(0, N)} else if(length(covariate) == 1){
    if(covariate == "c"){covariate <- sample(c(0, 1), N, replace = T)} 
    else if(covariate == "g" | covariate == TRUE){covariate <- rnorm(N)}
    else if(covariate == "rmean"){covariate <- rnorm(n = N, mean = runif(1))}
    else if(covariate == "ar"){
      covariate <- rnorm(1)
      if(is(eps, "matrix")){eps <- 0.2}
      for(i in 2:N){covariate[i] <- .1 * covariate[i - 1] + 1.5 * rnorm(n = 1, sd = sqrt(eps))}
    }
  }
  covariate2 <- matrix(NA, N, p)
  covariate3 <- covariate
  for(i in 1:p){covariate2[, i] <- covariate}
  if(class(b1) != "matrix"){
    if(length(b1) != p^2){stop("Length of b1 needs to equal p^2")}
    b1 <- sapply(b1, rep, N)
  }
  b2 <- sapply(b2, rep, N)
  c2 <- sapply(c2, rep, N)
  if(tvInit){
    pphi <- kronecker(matrix(b1[1, ], p, p, byrow = T), matrix(b1[1, ], p, p, byrow = T))
    mu1 <- solve(diag(p) - matrix(b1[1, ], p, p, byrow = T)) %*% matrix(b0[1, ], p, 1)
    s1 <- matrix(solve(diag(p * p) - pphi) %*% matrix(sigma2, ncol = 1), p, p)
  } else {
    mu1 <- rep(0, p)
    s1 <- sigma2
  }
  y[1, ] <- rockchalk::mvrnorm(1, mu = mu1, Sigma = s1)
  for(t in 2:N){
    a0 <- matrix(b0[t, ], p, 1)
    a1 <- matrix(b1[t, ], p, p, T) %*% y[t - 1, ]
    a2 <- matrix(b2[t, ], p, (p - 1), T) %*% (y[t - 1, 1] * y[t - 1, -1])
    a3 <- matrix(c1 * covariate2[t - 1, ], p, 1)
    a4 <- matrix(c2[t, ], p, p, T) %*% (y[t - 1, ] * covariate2[t - 1, ])
    e <- matrix(rockchalk::mvrnorm(1, rep(0, p), sigma2), p, 1)
    y[t, ] <- a0 + a1 + a2 + a3 + a4 + e
    if(any(is.na(y) | is.infinite(y))){stop(paste0("Failed at t = ", t))}
  }
  colnames(y) <- paste0("y", 1:p)
  if(all(covariate == 0)){covariate <- NULL} else {covariate <- covariate[-N]}
  betas <- list(b0 = matrix(b0[1, ], p, 1), b1 = matrix(b1[1, ], p, p, T),
                b2 = matrix(b2[1, ], p, p - 1, T), c1 = matrix(c1, p, 1),
                c2 = matrix(c2[1, ], p, p, T))
  if(all(b0 == 0)){betas$b0 <- NULL}
  if(all(b2 == 0)){betas$b2 <- NULL}
  if(all(c1 == 0)){betas$c1 <- NULL}
  if(all(c2 == 0)){betas$c2 <- NULL}
  out <- list(Y = y[-(1:(burnin + 1)), ], X = data.frame(y[-c(1:burnin, N), ]), betas = betas)
  colnames(out$X) <- paste0(colnames(out$X), "L")
  if(p == 1){colnames(out$X) <- "y1L"}
  if(!is.null(covariate)){out$X <- data.frame(out$X, covariate[-(1:burnin)])}
  if(getDat){
    if(p == 1){
      out$dat <- data.frame(y1 = out$Y, out$X)
    } else {
      out$dat <- list()
      for(i in 1:p){
        out$dat[[i]] <- data.frame(y = out$Y[, i], out$X)
        colnames(out$dat[[i]])[1] <- paste0("y", i)
      }
      names(out$dat) <- paste0("dat", 1:p)
    }
  }
  if(fullDat){
    if(!is.null(covariate)){bottom <- c(out$Y[(N - burnin) - 1, ], covariate3[N])}
    else{bottom <- c(out$Y[(N - burnin) - 1, ])}
    full <- data.frame(rbind(out$X, bottom))
    colnames(full) <- paste0("V", 1:ncol(full))
    out$full <- full
  }
  out
}

##### MLVARsim: simulate VAR data for multiple subjects
MLVARsim <- function(n, t, covariate = "rmean", intercept = "random", p = 2, 
                     b1 = NULL, c1 = NULL, c2 = NULL, getBetas = FALSE, eps = 0.2){
  covariate <- match.arg(covariate, c("true", "rmean", "c", "g", "false", "ar"))
  if(covariate == "true"){covariate <- as.logical(covariate)}
  if(covariate == "false"){covariate <- as.logical(covariate)}
  if(class(intercept) != "character"){b0 <- intercept}
  if(is.null(b1)){b1 <- c(.3, .2, .1, .4)} else if(all(b1 == 0)){b1 <- rep(0, p * p)}
  if(is.null(c1)){c1 <- c(.2, .15)} else if(all(c1 == 0)){c1 <- rep(0, p)}
  if(is.null(c2)){c2 <- c(0, .2, 0, .2)} else if(all(c2 == 0)){c2 <- rep(0, p * p)}
  if(getBetas){
    B <- rbind(cbind(matrix(b1, p, p, T), matrix(c1, p, 1), matrix(c2, p, p, T)), 0)
    colnames(B) <- c("V1.1", "V2.1", "V3.1", "V1.1:V3.1", "V2.1:V3.1")
    rownames(B) <- c("V1", "V2", "V3")
    return(B)
  }
  stopifnot(length(b1) == p * p)
  stopifnot(length(c1) == p)
  stopifnot(length(c2) == p * p)
  subjectData <- list()
  for(i in 1:n){
    if(class(intercept) == "character"){b0 <- rnorm(2)}
    dat <- VARsim(N = t, p = p, b1 = b1, covariate = covariate, b0 = b0, 
                  c1 = c1, c2 = c2, eps = eps, fullDat = TRUE)
    subjectData[[i]] <- dat$full
  }
  if(n == 1){return(dat)}
  ids <- rep(1:n, each = t)
  finalDat <- cbind.data.frame(do.call(rbind, subjectData), ID = ids)
  vars <- colnames(finalDat)[!colnames(finalDat) %in% "ID"]
  out <- list(data = finalDat, vars = vars, idvar = "ID")
  out
}


### ======================================================================== ###
### =========================== Fit GVAR models ============================ ###
### ======================================================================== ###
##### setupVAR: prepare data for analysis
setupVAR <- function(data, idvar = NULL, method = c("gvar", "lmer", "all"),
                     center = TRUE, scale = TRUE, vars = NULL, 
                     centerWithin = TRUE, scaleWithin = FALSE){
  method <- match.arg(method)
  if(class(data) == "list"){
    if(!"data" %in% names(data)){stop("Must supply data frame")}
    data <- as.data.frame(data[["data"]])
  }
  if(is.null(idvar)){
    if(!any(grepl("ID", colnames(data)))){stop("Must supply 'idvar'")}
    idvar <- colnames(data)[grep("ID", colnames(data))]
  }
  data <- data.frame(data[, -which(colnames(data) == idvar)], ID = data[, idvar])
  if(is.null(vars)){vars <- colnames(data)[!colnames(data) %in% "ID"]}
  ids <- unique(data[, "ID"])
  binary <- apply(data[, -which(colnames(data) == "ID")], 2, function(z) length(table(z)) <= 2)
  binary <- ifelse(any(binary), list(names(which(binary))), list(NULL))[[1]]
  vars0 <- setdiff(vars, binary)
  if(center){data[, vars0] <- apply(data[, vars0], 2, scale, center, scale)}
  dataByID <- lapply(ids, function(z) data[data[, "ID"] == z, seq_along(vars)])
  N <- rep(seq_along(ids), (sapply(dataByID, nrow) - ifelse(method == "all", 0, 1)))
  dataMeans <- do.call(rbind, lapply(dataByID, colMeans))[N, ]
  colnames(dataMeans) <- paste0(colnames(dataMeans), ".m")
  if(method != "lmer"){
    data0 <- lapply(dataByID, function(z){if(centerWithin){
      z[, vars0] <- apply(z[, vars0], 2, scale, centerWithin, scaleWithin)}
      return(z)
    })
    if(method == "all"){return(data.frame(do.call(rbind, data0), dataMeans, ID = ids[N]))}
    Y <- do.call(rbind, lapply(data0, function(z) z[-1, ]))
    X <- do.call(rbind, lapply(data0, function(z) z[-nrow(z), ]))
  } else {
    Y <- do.call(rbind, lapply(dataByID, function(z) z[-1, ]))
    X <- do.call(rbind, lapply(dataByID, function(z){
      z <- z[-nrow(z), ]
      if(centerWithin){z[, vars0] <- apply(
        z[, vars0], 2, scale, centerWithin, scaleWithin)}
      return(z)
    }))
  }
  dat <- data.frame(Y, X, dataMeans, ID = ids[N])
  dat
}

##### mlGVAR: fit GVAR models with multilevel data
mlGVAR <- function(data, m = NULL, selectFUN = NULL, subjectNets = FALSE, idvar = 'ID',
                   exogenous = TRUE, center = TRUE, scale = TRUE, fixedType = 'g', 
                   betweenType = 'g', centerWithin = TRUE, scaleWithin = FALSE,
                   rule = 'OR', threshold = 'none', verbose = TRUE, pcor = FALSE, 
                   fixedArgs = NULL, betweenArgs = NULL, bm = FALSE, ...){
  t1 <- Sys.time()
  mnames <- mi <- m
  args0 <- list(...)
  if(!bm){
    bm <- list(moderators = NULL)
    betweenArgs <- switch(
      2 - is.null(betweenArgs), bm, 
      if(!'moderators' %in% names(betweenArgs)){append(betweenArgs, bm)} else {betweenArgs})
  }
  if(is.character(threshold)){
    threshold <- sapply(match.arg(tolower(
      threshold), c('none', 'pcc', 'between', 'fixed', 'all'), several.ok = TRUE), 
      function(z) switch(z, none = FALSE, all = TRUE, z), USE.NAMES = FALSE)
  }
  if(!idvar %in% colnames(data)){stop('Must supply idvar')}
  if(!is.null(m)){
    mnames <- switch(2 - is.character(m), m, colnames(data)[m])
  } else if(!is.null(selectFUN)){
    a0 <- names(args0)
    if(!'method' %in% a0){
      args0$method <- 'glmnet'
    } else if(args0$method == 'glinternet'){
      mnames <- m <- setdiff(colnames(data), idvar)
    }
    if(!'criterion' %in% a0){args0$criterion <- 'AIC'}
    if(any(c('resample', 'split') %in% selectFUN) & !'sampMethod' %in% a0){
      args0$sampMethod <- 'split'
    }
  }
  data <- data.frame(data[, -which(colnames(data) == idvar)], ID = data[, idvar])
  vars <- setdiff(colnames(data), 'ID')
  dat <- setupVAR(data = data, idvar = 'ID', method = 'all', center = center,
                  scale = scale, vars = vars, centerWithin = centerWithin,
                  scaleWithin = scaleWithin)
  fixedDat <- dat[, vars]
  samp_ind <- as.numeric(cumsum(table(dat[, 'ID'])))
  attr(fixedDat, 'samp_ind') <- (1:nrow(dat))[-samp_ind]
  if(!is.null(m)){
    m <- which(vars %in% mnames)
    if(length(m) >= ncol(fixedDat) - 1){exogenous <- FALSE}
    if(!is.null(selectFUN)){args0$method <- 'glinternet'}
  }
  if(verbose){message('Estimating fixed networks')}
  fixedThresh <- ifelse(!is.character(threshold), threshold, ifelse(
    'fixed' %in% threshold, TRUE, ifelse('pcc' %in% threshold, 'PCC', FALSE)))
  fitNetArgs <- setdiff(formalArgs('fitNetwork'), '...')
  args1 <- list(data = fixedDat, moderators = m, type = fixedType, lags = 1,
                exogenous = exogenous, center = FALSE, scale = FALSE, pcor = pcor,
                rule = rule, threshold = fixedThresh, verbose = verbose)
  if(length(args0) > 0){args0 <- args0[setdiff(names(args0), names(args1))]}
  if(!is.null(fixedArgs)){
    fix1 <- intersect(names(fixedArgs), names(args1))
    if(length(fix1) > 0){args1 <- replace(args1, fix1, fixedArgs[fix1])}
    fixedArgs <- fixedArgs[setdiff(names(fixedArgs), names(args1))]
    args1 <- append(args1, fixedArgs[intersect(fitNetArgs, names(fixedArgs))])
  }
  if(is.null(selectFUN)){
    args1 <- append(args1, args0[intersect(fitNetArgs, names(args0))])
    fixedNets <- do.call(fitNetwork, args1)
  } else {
    if(isTRUE(selectFUN)){selectFUN <- 'varSelect'}
    if(length(selectFUN) == 2){betweenType <- isTRUE(selectFUN[2] == FALSE)}
    selectFUN <- match.arg(selectFUN[1], c('varSelect', 'resample', 'stability', 'bootstrap', 'split'))
    if(!selectFUN %in% c('varSelect', 'resample')){
      args0$sampMethod <- selectFUN
      selectFUN <- 'resample'
    }
    FUNargs <- setdiff(formalArgs(selectFUN), '...')
    FUN <- match.fun(selectFUN)
    args1 <- setNames(args1, gsub('moderators', 'm', names(args1)))
    args1.1 <- append(args1, args0)
    args1.1 <- args1.1[intersect(FUNargs, names(args1.1))]
    if('criterion' %in% names(args1.1)){
      if(length(args1.1$criterion) > 1){
        args1.1$criterion <- args1.1$criterion[1]}}
    if('gamma' %in% names(args1.1)){
      if(length(args1.1$gamma) > 1){args1.1$gamma <- args1.1$gamma[1]}
    }
    fixedType <- tryCatch({do.call(FUN, args1.1)}, error = function(e){TRUE})
    if(selectFUN == 'varSelect' | isTRUE(fixedType)){
      args1.2 <- append(replace(args1, c('type', 'verbose'), list(
        type = 'g', verbose = FALSE)), args0)
      names(args1.2)[names(args1.2) == 'm'] <- 'moderators'
      args1.2 <- args1.2[intersect(fitNetArgs, names(args1.2))]
      fixedNets <- tryCatch({do.call(fitNetwork, replace(args1.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, replace(args1.2, c('type', 'saveMods'), c('g', FALSE)))})
    } else {
      if('fit0' %in% names(fixedType)){fixedType$fit0 <- NULL}
      if(is(fixedType, 'list')){attr(fixedType[[2]], 'threshold') <- fixedThresh}
      args1.2 <- append(list(obj = fixedType, data = args1.1$data, fit = TRUE), 
                        args0[intersect(c('select', 'thresh'), names(args0))])
      fixedNets <- tryCatch({do.call(modSelect, replace(args1.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, setNames(args1, gsub('^m$', 'moderators', names(args1))))
      })
      if(verbose){cat('\n')}
    }
  }
  ids <- unique(dat[, 'ID'])
  meansDat <- data.frame(do.call(rbind, lapply(ids, function(z){
    z <- dat[dat[, 'ID'] == z, paste0(vars, '.m')]
    return(z[1, ])
  })))
  colnames(meansDat) <- vars
  rownames(meansDat) <- 1:nrow(meansDat)
  if(verbose){message('Estimating between-subjects network')}
  if(nrow(meansDat) <= ncol(meansDat)){
    warning('Not enough subjects to fit unregularized between-subjects network')
  }
  betThresh <- ifelse(!is.character(threshold), threshold, 'between' %in% threshold)
  args2 <- list(data = meansDat, moderators = m, type = betweenType, center = FALSE,
                scale = FALSE, exogenous = exogenous, lags = NULL, rule = rule, 
                threshold = betThresh, pcor = pcor, verbose = verbose)
  if(!is.null(betweenArgs)){
    if('covariates' %in% names(betweenArgs) & !'moderators' %in% names(betweenArgs)){
      betweenArgs <- append(betweenArgs, list(moderators = NULL))
    }
    fix2 <- intersect(names(betweenArgs), names(args2))
    if(length(fix2) > 0){
      if('moderators' %in% fix2){
        if(is.null(betweenArgs$moderators) & !is.null(args2$moderators)){
          if(!'covariates' %in% names(betweenArgs) & exogenous){
            args2$data <- args2$data[, -args2$moderators]
          }
        }
      }
      args2 <- replace(args2, fix2, betweenArgs[fix2])
    }
    betweenArgs <- betweenArgs[setdiff(names(betweenArgs), names(args2))]
    args2 <- append(args2, betweenArgs[intersect(fitNetArgs, names(betweenArgs))])
  }
  if(is.null(selectFUN) | isTRUE(betweenType)){
    if(isTRUE(betweenType)){args2$type <- betweenType <- 'g'}
    args2 <- append(args2, args0[intersect(fitNetArgs, names(args0))])
    betNet <- do.call(fitNetwork, args2)
  } else {
    if(is.logical(betweenType)){
      betweenType <- 'g'
      selectFUN <- 'varSelect'
      FUNargs <- formalArgs('varSelect')
      FUN <- match.fun('varSelect')
    }
    args2 <- setNames(args2, gsub('moderators', 'm', names(args2)))
    args2.1 <- append(args2, args0)
    if(is.null(args2.1$m) & args2.1$method == 'glinternet'){args2.1$method <- 'glmnet'}
    args2.1 <- args2.1[intersect(FUNargs, names(args2.1))]
    if('criterion' %in% names(args2.1)){
      if(length(args2.1$criterion) > 1){
        args2.1$criterion <- args2.1$criterion[2]}}
    if('gamma' %in% names(args2.1)){
      if(length(args2.1$gamma) > 1){args2.1$gamma <- args2.1$gamma[2]}
    }
    betweenType <- tryCatch({do.call(FUN, args2.1)}, error = function(e){TRUE})
    if(selectFUN == 'varSelect' | isTRUE(betweenType)){
      args2.2 <- append(replace(args2, c('type', 'verbose'), list(
        type = 'g', verbose = FALSE)), args0)
      names(args2.2)[names(args2.2) == 'm'] <- 'moderators'
      args2.2 <- args2.2[intersect(fitNetArgs, names(args2.2))]
      betNet <- tryCatch({do.call(fitNetwork, replace(args2.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, replace(args2.2, c('type', 'saveMods'), c('g', FALSE)))})
    } else {
      if('fit0' %in% names(betweenType)){betweenType$fit0 <- NULL}
      if(is(betweenType, 'list')){attr(betweenType[[2]], 'threshold') <- betThresh}
      args2.2 <- append(list(obj = betweenType, data = args2.1$data, fit = TRUE), 
                        args0[intersect(c('select', 'thresh'), names(args0))])
      betNet <- tryCatch({do.call(modSelect, replace(args2.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, setNames(args2, gsub('^m$', 'moderators', names(args2))))
      })
    }
  }
  if(all(subjectNets != FALSE)){
    subjectSelect <- ifelse(all(subjectNets == 'select'), TRUE, FALSE)
    if(isTRUE(subjectNets) | isTRUE(subjectSelect)){subjectNets <- ids}
    if(any(!subjectNets %in% ids)){subjectNets <- intersect(subjectNets, ids)}
    if(verbose){
      message('Estimating subject-specific networks')
      pb <- txtProgressBar(max = length(subjectNets), style = 3)
    }
    indNets <- lapply(seq_along(subjectNets), function(i){
      dati <- dat[dat[, 'ID'] == subjectNets[i], vars]
      if(!is.null(m)){
        mi <- sapply(seq_along(m), function(j) length(unique(dati[, m[j]])) != 1)
        mi <- ifelse(all(!mi), list(NULL), list(m[mi]))[[1]]
        if(!identical(m, mi)){dati <- dati[, -setdiff(m, mi)]}
      }
      p0 <- ncol(dati)
      fiti <- nrow(dati) > (2 * p0) + 1
      if(!is.null(mi)){
        p1 <- sum(sapply(seq_along(mi), function(p) p0 - p))
        fiti <- ifelse(exogenous & length(mi) < ncol(dati) - 1, 
                       nrow(dati) > p1 + (2 * p0) - (length(mi) - 1), 
                       nrow(dati) > p1 + (2 * p0) + 1)
      }
      if(fiti){
        indi <- fitNetwork(data = dati, moderators = mi, threshold = fixedThresh,
                           exogenous = exogenous, center = center, type = 'g', 
                           lags = 1, scale = ifelse(center, scale, FALSE), 
                           saveMods = FALSE, ...)
      }
      if(verbose){setTxtProgressBar(pb, i)}
      if(fiti){return(indi)} else {return(list())}
    })
    names(indNets) <- paste0('subject', subjectNets)
    if(any(sapply(indNets, length) == 0)){
      inderrs <- subjectNets[sapply(indNets, length) == 0]
      indNets <- indNets[which(sapply(indNets, length) > 0)]
      if(length(indNets) > 0){
        message(paste0('Too few time points to estimate networks for subject',
                       ifelse(length(inderrs) == 1, ': ', 's: '),
                       paste(inderrs, collapse = ', ')))
      } else {
        subjectNets <- FALSE
        message('Too few time points to estimate subject-specific networks')
      }
    }
  }
  outcall <- list(fixedType = fixedType, betweenType = betweenType, m = mnames, 
                  center = center, scale = scale, exogenous = exogenous, 
                  threshold = threshold, centerWithin = centerWithin, 
                  scaleWithin = scaleWithin, rule = rule)
  if(length(args0) > 0){outcall <- append(outcall, args0)}
  if(!is.null(fixedArgs) & length(fixedArgs) > 0){outcall <- append(outcall, list(fixedArgs = fixedArgs))}
  if(!is.null(betweenArgs) & length(betweenArgs) > 0){outcall <- append(outcall, list(betweenArgs = betweenArgs))}
  out <- list(call = outcall, fixedNets = fixedNets, betweenNet = betNet, ids = ids)
  if(!is.null(selectFUN)){
    out$ids <- NULL
    out$varMods <- list(fixedMods = fixedType, betweenMods = betweenType)
    out$call$fixedType <- fixedNets$call$type
    out$call$betweenType <- betNet$call$type
    out$call$selectFUN <- selectFUN
    out$ids <- ids
  }
  if(all(subjectNets != FALSE)){out$subjectNets <- indNets}
  out$netData <- list(samp_ind = samp_ind, fixedDat = fixedDat, meansDat = meansDat)
  out$dat <- dat
  attr(out, 'mlGVAR') <- TRUE
  class(out) <- c('list', 'mlGVAR')
  attr(out, 'time') <- t2 <- Sys.time() - t1
  if(verbose){cat('\nCOMPLETE:', round(t2, 2), attr(t2, 'units')); cat('\n')}
  return(out)
}


### ======================================================================== ###
### ========================= Fit lmerVAR models =========================== ###
### ======================================================================== ###
##### lmerVAR: fit lmerVAR models
lmerVAR <- function(data, m = NULL, temporal = "default", contemp = "default",
                    idvar = "ID", intvars = NULL, center = TRUE, scale = TRUE, 
                    centerWithin = TRUE, scaleWithin = FALSE, exogenous = TRUE,
                    covariates = NULL, fix = NULL, warnings = FALSE, 
                    verbose = TRUE){
  t1 <- Sys.time()
  suppressMessages(invisible(c(require(lme4), require(lmerTest))))
  if(!warnings){oldw <- getOption("warn"); options(warn = -1)}
  mnames <- m
  data <- data.frame(data)
  vars <- colnames(data)
  if(!idvar %in% vars){stop("Must supply 'idvar'")}
  if(!is.null(m)){mnames <- switch(2 - is.character(m), m, vars[m])}
  data <- data.frame(data[, -which(vars == idvar)], ID = data[, idvar])
  vars <- setdiff(vars, idvar)
  dat <- setupVAR(data = data, idvar = "ID", method = "all", center = center,
                  scale = scale, centerWithin = FALSE, scaleWithin = FALSE)
  dat0 <- dat[, vars]
  samp_ind <- as.numeric(cumsum(table(dat[, "ID"])))
  nid <- attr(dat0, "samp_ind") <- setdiff(1:nrow(dat), samp_ind)
  ids0 <- unique(dat[, "ID"])
  ids <- dat[nid, "ID"]
  if(!is.null(m)){
    m <- match(mnames, vars)
    exogenous <- ifelse(length(m) >= length(vars) - 1, FALSE, exogenous)
  }
  dat0 <- lagMat(data = dat0, m = m, center = FALSE, scale = FALSE, 
                 exogenous = exogenous, lags = 1, checkType = TRUE, 
                 covariates = covariates)
  if(centerWithin){
    dat0$X[, vars] <- do.call(rbind, lapply(seq_along(ids0), function(i){
      apply(subset(dat0$X[, vars], ids == ids0[i]), 2, scale, TRUE, scaleWithin)
    }))
  }
  yvars <- colnames(dat0$Y)
  mvars <- paste0(vars, ".m")
  dat1 <- data.frame(dat0$Y, dat0$X, dat[nid, mvars], ID = ids, check.names = FALSE)
  if(is.null(intvars) & !is.null(m)){intvars <- colnames(dat1)[grep(":", colnames(dat1))]}
  temporal <- match.arg(temporal, c("default", "correlated", "orthogonal", "fixed", "intfixed"))
  if(temporal == "default"){temporal <- ifelse(length(yvars) > 6, "orthogonal", "correlated")}
  if(verbose){
    message(paste0("Estimating temporal and between-subject networks (", temporal, ")"))
    pb <- txtProgressBar(min = 0, max = length(yvars), style = 3)
  }
  tempForms <- lapply(yvars, function(y){
    x <- paste0(vars, collapse = " + ")
    means <- paste0(paste0(setdiff(
      gsub("[.]m$", "", mvars), 
      gsub("[.]y$", "", y)), ".m"), 
      collapse = " + ")
    ints <- paste0(intvars, collapse = " + ")
    ints <- switch(2 - isTRUE(ints == ""), NULL, ints)
    fixed <- paste0(c(x, means, ints), collapse = " + ")
    if(temporal == "correlated"){
      rands <- paste0("(", paste0(c(x, ints), collapse = " + "), " | ID)")
    } else {
      rands <- "(1 | ID)"
      randvars <- switch(2 - isTRUE(temporal == "intfixed"), vars, c(vars, intvars))
      if(!is.null(fix)){randvars <- setdiff(randvars, fix)}
      if(temporal != "fixed"){
        rands <- paste0("(", rands, " + ", paste0(
          paste0("(0 + ", randvars, " | ID)"), 
          collapse = " + "), ")")
      }
    }
    return(as.formula(paste0(y, " ~ ", fixed, " + ", rands)))
  })
  tempMods <- setNames(lapply(seq_along(yvars), function(z){
    tm <- suppressMessages(lmer(tempForms[[z]], data = dat1, REML = FALSE))
    if(verbose){setTxtProgressBar(pb, z)}
    return(tm)
  }), yvars)
  resDat <- structure(do.call(data.frame, lapply(tempMods, resid)), names = yvars)
  contemp <- match.arg(contemp, c("default", "correlated", "orthogonal"))
  if(contemp == "default"){contemp <- ifelse(length(yvars) > 6, "orthogonal", "correlated")}
  if(verbose){
    message(paste0("\nEstimating contemporaneous network (", contemp, ")"))
    pb <- txtProgressBar(min = 0, max = length(yvars), style = 3)
  }
  contempForms <- lapply(yvars, function(y){
    x <- setdiff(yvars, y)
    fixed <- paste0(c(0, x), collapse = " + ")
    if(contemp == "correlated"){
      rands <- paste0("(", fixed, " | ID)")
    } else {
      rands <- paste0("((1 | ID) + ", paste0(
        paste0("(0 + ", x, " | ID)"), collapse = " + "), ")")
    }
    return(as.formula(paste0(y, " ~ ", fixed, " + ", rands)))
  })
  resDat$ID <- ids
  contempMods <- setNames(lapply(seq_along(yvars), function(z){
    cm <- suppressMessages(lmer(contempForms[[z]], data = resDat, REML = FALSE))
    if(verbose){setTxtProgressBar(pb, z)}
    return(cm)
  }), yvars)
  fit <- do.call(rbind, lapply(tempMods, function(z) c(AIC(z), BIC(z))))
  fit <- structure(cbind.data.frame(yvars, fit), names = c("var", "aic", "bic"))
  inds <- list(yvars = yvars, vars = vars, mvars = mvars, intvars = intvars)
  outcall <- list(m = mnames, temporal = temporal, contemp = contemp, 
                  exogenous = exogenous, center = center, scale = scale, 
                  centerWithin = centerWithin, scaleWithin = scaleWithin)
  if(!is.null(covariates)){outcall <- append(outcall, list(covariates = covariates))}
  model <- list(tempMods = tempMods, contempMods = contempMods, fit = fit)
  out <- tryCatch({lmerNets(model = model, inds = inds, m = mnames)}, 
                  error = function(e){list(inds = inds)})
  for(i in seq_along(yvars)){
    attributes(model$tempMods[[i]])$formula <- tempForms[[i]]
    attributes(model$contempMods[[i]])$formula <- contempForms[[i]]
  }
  out <- append(list(call = outcall), append(out, list(mods = model, data = dat1)))
  attr(out, "temporal") <- paste0(temporal, ifelse(
    !is.null(intvars), " (interaction)", ifelse(
      !is.null(covariates), " (covariate)", "")))
  attr(out, "contemporaneous") <- contemp
  attr(out, "lmerVAR") <- TRUE
  class(out) <- c('list', 'lmerVAR')
  attr(out, "time") <- t2 <- Sys.time() - t1
  if(verbose){cat("\n"); print(Sys.time() - t1)}
  if(!warnings){options(warn = oldw)}
  return(out)
}

##### lmerNets: create networks from lmerVAR models
lmerNets <- function(model, inds, m = NULL, threshold = FALSE, 
                     rule = "OR", ggm = "pcor"){
  rule <- match.arg(tolower(rule), c("or", "and"))
  ggm <- match.arg(tolower(ggm), c("pcor", "cor", "cov", "prec"))
  forcePositive <- function(x){
    x <- (x + t(x))/2
    if(any(eigen(x)$values < 0)){
      x <- x - diag(nrow(x)) * min(eigen(x)$values) - 0.001
    }
    return(x)
  }
  y <- inds$yvars
  x <- inds$vars
  mvars <- inds$mvars
  intvars <- inds$intvars
  k <- length(y)
  y1 <- list(y, x)
  y2 <- rep(list(y), 2)
  ### BETA
  beta <- beta0 <- do.call(rbind, lapply(
    model$tempMods, function(z) fixef(z)[x]))
  betaSE <- betaSE0 <- do.call(rbind, lapply(
    model$tempMods, function(z) arm::se.fixef(z)[x]))
  beta_pvals <- betaPs0 <- (1 - pnorm(abs(beta/betaSE))) * 2
  dimnames(beta) <- dimnames(betaSE) <- dimnames(beta_pvals) <- y1
  dimnames(beta0) <- dimnames(betaSE0) <- dimnames(betaPs0) <- y1
  if(ncol(beta) != nrow(beta)){
    beta <- beta[y, match(y, paste0(x, ".y"))]
    betaSE <- betaSE[y, match(y, paste0(x, ".y"))]
    beta_pvals <- beta_pvals[y, match(y, paste0(x, ".y"))]
  }
  ### GAMMA THETA
  gammaTheta <- lapply(model$contempMods, fixef)
  gammaThetaSE <- lapply(model$contempMods, arm::se.fixef)
  gt1 <- gt2 <- matrix(NA, k, k)
  for(i in 1:k){
    gt1[i, match(names(gammaTheta[[i]]), y)] <- gammaTheta[[i]]
    gt2[i, match(names(gammaThetaSE[[i]]), y)] <- gammaThetaSE[[i]]
  }
  diag(gt1) <- diag(gt2) <- 0
  gammaTheta <- gt1
  gammaThetaSE <- gt2
  gammaTheta_pvals <- (1 - pnorm(abs(gammaTheta/gammaThetaSE))) * 2
  diag(gammaTheta_pvals) <- 1
  dimnames(gammaTheta) <- dimnames(gammaThetaSE) <- dimnames(gammaTheta_pvals) <- y2
  ### THETA AND PDC
  D <- diag(1/sapply(model$contempMods, sigma)^2)
  inv <- tryCatch({forcePositive(D %*% (diag(k) - gammaTheta))}, error = function(e){diag(k)})
  if(identical(inv, diag(k))){message("\nNon-convergent estimate of Theta")}
  Theta_prec <- inv
  Theta_cov <- corpcor::pseudoinverse(inv)
  Theta_cor <- cov2cor(Theta_cov)
  d <- 1/sqrt(diag(inv))
  Theta <- -t(d * inv) * d
  diag(Theta) <- 0
  dimnames(Theta) <- dimnames(Theta_cov) <- dimnames(Theta_cor) <- dimnames(Theta_prec) <- y2
  PDC <- beta/(sqrt(diag(Theta_cov) %o% diag(Theta_prec) + beta^2))
  colnames(beta) <- colnames(PDC) <- colnames(betaSE) <- colnames(beta_pvals) <- paste0(colnames(beta), ".lag1.")
  ### GAMMA OMEGA
  between <- TRUE
  if(length(y) != length(x)){
    if(!is.null(m)){
      if(length(setdiff(x, m)) < 3){between <- FALSE}
      if(between){mvars <- setdiff(mvars, paste0(m, ".m"))}
    } else {
      mvars <- mvars[gsub("[.]m$", ".y", mvars) %in% y]
    }
  }
  gammaOmega <- lapply(model$tempMods, function(z){
    z2 <- fixef(z)
    z2 <- z2[intersect(names(z2), mvars)]
    return(z2)
  })
  gammaOmegaSE <- lapply(model$tempMods, function(z){
    z2 <- arm::se.fixef(z)
    z2 <- z2[intersect(names(z2), mvars)]
    return(z2)
  })
  if(between){
    go1 <- go2 <- matrix(NA, k, k)
    for(i in 1:k){
      go1[i, match(names(gammaOmega[[i]]), mvars)] <- gammaOmega[[i]]
      go2[i, match(names(gammaOmegaSE[[i]]), mvars)] <- gammaOmegaSE[[i]]
    }
    diag(go1) <- diag(go2) <- 0
    gammaOmega <- go1
    gammaOmegaSE <- go2
    gammaOmega_pvals <- (1 - pnorm(abs(gammaOmega/gammaOmegaSE))) * 2
    diag(gammaOmega_pvals) <- 1
    dimnames(gammaOmega) <- dimnames(gammaOmegaSE) <- dimnames(gammaOmega_pvals) <- y2
    ### OMEGA
    mu_SD <- sapply(model$tempMods, function(z) attr(lme4::VarCorr(z)[[1]], "stddev")[1])
    D <- diag(1/mu_SD^2)
    inv <- tryCatch({forcePositive(D %*% (diag(k) - gammaOmega))}, error = function(e){diag(k)})
    if(identical(inv, diag(k))){message("\nNon-convergent estimate of Omega")}
    Omega_prec <- inv
    Omega_cov <- corpcor::pseudoinverse(inv)
    Omega_cor <- cov2cor(Omega_cov)
    d <- 1/sqrt(diag(inv))
    Omega <- -t(d * inv) * d
    diag(Omega) <- 0
  } else {
    gammaOmega <- gammaOmegaSE <- gammaOmega_pvals <- Omega <- Omega_cov <- Omega_cor <- Omega_prec <- diag(0, k)
    dimnames(gammaOmega) <- dimnames(gammaOmegaSE) <- dimnames(gammaOmega_pvals) <- y2
  }
  dimnames(Omega) <- dimnames(Omega_cov) <- dimnames(Omega_cor) <- dimnames(Omega_prec) <- y2
  ### COLLECT RESULTS
  out <- list(Beta = list(mu = beta, SE = betaSE, Pvals = beta_pvals), 
              Theta = list(cor = Theta_cor, cov = Theta_cov, pcor = Theta, prec = Theta_prec), 
              Omega = list(cor = Omega_cor, cov = Omega_cov, pcor = Omega, prec = Omega_prec), 
              gammaTheta = list(mu = gammaTheta, SE = gammaThetaSE, Pvals = gammaTheta_pvals),
              gammaOmega = list(mu = gammaOmega, SE = gammaOmegaSE, Pvals = gammaOmega_pvals))
  if(!is.null(m)){
    ints <- do.call(rbind, lapply(model$tempMods, function(z){
      fixef(z)[grepl(":", names(fixef(z)))]}))
    intsSE <- do.call(rbind, lapply(model$tempMods, function(z){
      arm::se.fixef(z)[grepl(":", names(arm::se.fixef(z)))]}))
    intsPvals <- (1 - pnorm(abs(ints/intsSE))) * 2
    out$ints <- list(coefs = ints, SE = intsSE, Pvals = intsPvals) 
  }
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    dimnames(colMat) <- dimnames(adjMat)
    colMat
  }
  ### THRESHOLDING
  if(threshold != FALSE){
    if(threshold == TRUE){threshold <- 0.05}
    thresh <- function(obj, alpha = 0.05, rule = "or", ggm = "pcor"){
      beta <- obj$Beta$mu
      theta <- obj$Theta[[ggm]]
      omega <- obj$Omega[[ggm]]
      bp <- obj$Beta$Pvals
      tp <- obj$gammaTheta$Pvals
      op <- obj$gammaOmega$Pvals
      diag(tp) <- diag(op) <- 1
      beta2 <- beta * ifelse(bp <= alpha, 1, 0)
      if("PDC" %in% names(obj)){pdc2 <- obj$PDC * ifelse(bp <= alpha, 1, 0)}
      if(rule == "or"){
        theta2 <- theta * ifelse(tp <= alpha | t(tp) <= alpha, 1, 0)
        omega2 <- omega * ifelse(op <= alpha | t(op) <= alpha, 1, 0)
      } else if(rule == "and"){
        theta2 <- theta * ifelse(tp <= alpha & t(tp) <= alpha, 1, 0)
        omega2 <- omega * ifelse(op <= alpha & t(op) <= alpha, 1, 0)
      }
      res <- list(beta = beta2, theta = theta2, omega = omega2)
      if("PDC" %in% names(obj)){res$PDC <- pdc2}
      return(res)
    }
    out3 <- thresh(obj = append(out, list(PDC = PDC)), 
                   alpha = threshold, rule = rule, ggm = ggm)
    beta <- out3$beta
    Theta <- out3$theta
    Omega <- out3$omega
    PDC <- out3$PDC
  }
  ### OUTPUT
  out2 <- list(temporal = list(adjMat = beta, edgeColors = getEdgeColors(beta)),
               contemporaneous = list(adjMat = Theta, edgeColors = getEdgeColors(Theta)),
               between = list(adjMat = Omega, edgeColors = getEdgeColors(Omega)))
  out2$temporal <- append(out2$temporal, list(PDC = list(adjMat = PDC, edgeColors = getEdgeColors(PDC))))
  out2$contemporaneous$kappa <- Theta_prec
  output <- append(out2, list(coefs = out))
  #output <- append(out2, list(coefs = out, mods = model[grepl("Mods|fit", names(model))]))
  attr(output, "lmerVAR") <- TRUE
  return(output)
}

##### compareVAR: compare lmerVAR models
compareVAR <- function(m1, m2, m3 = NULL, anova = NULL, type = "tempMods"){
  if(is.null(m3)){
    if(all(unlist(lapply(list(m1, m2), function(z) "temporal" %in% names(attributes(z)))))){
      cat(paste0("Model 1: ", attributes(m1)$temporal), "\n")
      cat(paste0("Model 2: ", attributes(m2)$temporal), "\n\n")
    }
    m1 <- m1$mods; m2 <- m2$mods
    aic <- apply(cbind(m1$fit$aic, m2$fit$aic), 1, which.min)
    bic <- apply(cbind(m1$fit$bic, m2$fit$bic), 1, which.min)
    if(all.equal(m1$fit$aic, m2$fit$aic) == TRUE){aic <- rep(0, nrow(m1$fit))}
    if(all.equal(m1$fit$bic, m2$fit$bic) == TRUE){bic <- rep(0, nrow(m1$fit))}
    return(data.frame(aic, bic))
  } else if(is.null(anova)){
    M <- list(m1, m2, m3)
    if(all(unlist(lapply(M, function(z) "temporal" %in% names(attributes(z)))))){
      names(M) <- c("Model 1", "Model 2", "Model 3")
      for(i in 1:3){cat(paste0(names(M)[i], ": ", lapply(M, attr, "temporal")[[i]]), "\n")}
      cat("\n")
    }
    m1 <- m1$mods; m2 <- m2$mods; m3 <- m3$mods
    vars <- m1$fit$var
    out <- list()
    for(i in 1:length(vars)){
      out[[i]] <- rbind(m1$fit[m1$fit$var == vars[i], ], m2$fit[m2$fit$var == vars[i], ],
                        m3$fit[m3$fit$var == vars[i], ])
    }
    return(cbind.data.frame(var = vars, do.call(rbind, lapply(out, function(z) apply(z[,-1], 2, which.min)))))
  } else if(is.numeric(anova) & anova <= length(m1$mods$tempMods)){
    m <- anova
    type <- match.arg(type, c("tempMods", "contempMods"))
    M1 <- m1$mods[[type]][[m]]
    M2 <- m2$mods[[type]][[m]]
    M3 <- m3$mods[[type]][[m]]
    return(anova(M1, M2, M3))
  }
}

compareVAR2 <- function(mods, p){
  vs <- paste0(names(mods), "$mods$tempMods[[", p, "]]")
  obj <- paste0("anova(", paste0(vs, collapse = ", "), ")")
  for(i in 1:length(vs)){assign(names(mods)[i], mods[[i]], pos = 1)}
  out <- eval(parse(text = obj))
  out
}


### ======================================================================== ###
### ============================= PLOTS ==================================== ###
### ======================================================================== ###
##### makePlot:
makePlot <- function(model, data, zxy = NULL, res = 50, radius = .03, 
                     points = TRUE, ppoints = TRUE, save = FALSE){
  if(is.null(zxy)){zxy <- as.list(colnames(data))}
  z <- zxy[[1]]
  x <- zxy[[2]]
  y <- zxy[[3]]
  if(is.null(names(zxy))){
    zlab <- z
    xlab <- x
    ylab <- y
  } else {
    zlab <- names(zxy)[1]
    xlab <- names(zxy)[2]
    ylab <- names(zxy)[3]
  }
  require(rgl)
  predictgrid <- function(model, xvar, yvar, zvar, res, type = NULL){
    xr <- range(model$model[[xvar]])
    yr <- range(model$model[[yvar]])
    newdata <- expand.grid(x = seq(xr[1], xr[2], length = res), y = seq(yr[1], yr[2], length = res))
    names(newdata) <- c(xvar, yvar)
    newdata[[zvar]] <- predict(model, newdata = newdata, type = type)
    newdata
  }
  df2mat <- function(p, xvar = NULL, yvar = NULL, zvar = NULL){
    if(is.null(xvar)){xvar <- names(p)[1]}
    if(is.null(yvar)){yvar <- names(p)[2]}
    if(is.null(zvar)){zvar <- names(p)[3]}
    x <- unique(p[[xvar]])
    y <- unique(p[[yvar]])
    z <- unique(p[[zvar]], nrow = length(y), ncol = length(x))
    m <- list(x, y, z)
    names(m) <- c(xvar, yvar, zvar)
    m
  }
  interleave <- function(v1, v2){as.vector(rbind(v1, v2))}
  data$pred <- predict(model)
  grid_df <- predictgrid(model, x, y, z, res = res)
  grid_list <- df2mat(grid_df)
  plot3d(data[,x], data[,y], data[,z], xlab = "", ylab = "", zlab = "", 
         axes = F, type = ifelse(points, "s", "n"), size = .4, lit = F)
  if(ppoints & points){
    spheres3d(data[,x], data[,y], data$pred, radius = radius, alpha = .4, type = "s", size = .5, lit = F)
    segments3d(interleave(data[,x], data[,x]), interleave(data[,y], data[,y]), 
               interleave(data[,z], data$pred), alpha = .4, col = "red")
  }
  surface3d(grid_list[[x]], grid_list[[y]], grid_list[[z]], 
            alpha = .4, front = "lines", back = "lines")
  rgl.bbox(color = "grey50", emission = "grey50", xlen = 0, ylen = 0, zlen = 0)
  rgl.material(color = "black")
  axes3d(edges = c("x--", "y+-", "z--"), ntick = 6, cex = .75)
  mtext3d(xlab, edge = "x--", line = 2)
  mtext3d(ylab, edge = "y+-", line = 3)
  mtext3d(zlab, edge = "z--", line = 3)
  if(save == TRUE){rgl.snapshot("~/Desktop/fuckery.png", fmt = "png")}
}

##### makePlot2:
makePlot2 <- function(data, zxy = NULL, center = TRUE, res = 50, effect = TRUE){
  if(is.null(zxy)){zxy <- as.list(colnames(data))}
  require(ggplot2)
  require(margins)
  if(is.null(names(zxy))){
    zlab <- zxy[[1]]
    xlab <- zxy[[2]]
    ylab <- zxy[[3]]
  } else {
    zlab <- names(zxy)[1]
    xlab <- names(zxy)[2]
    ylab <- names(zxy)[3]
  }
  z <- data[,zxy[[1]]]
  x <- data[,zxy[[2]]]
  y <- data[,zxy[[3]]]
  if(center){
    x <- x - mean(x)
    y <- y - mean(y)
  }
  dat <- data.frame(z, x, y)
  mod <- lm(z ~ x * y, data = dat)
  xvals <- seq(min(x), max(x), length = res)
  lowy <- mean(y) - sd(y)
  highy <- mean(y) + sd(y)
  lowdf <- data.frame(x = xvals, y = lowy)
  highdf <- data.frame(x = xvals, y = highy)
  pred_low <- predict(mod, data.frame(x = lowdf$x, y = lowdf$y), se.fit = T)
  yhat_low <- pred_low$fit
  se_low <- pred_low$se.fit * qnorm(.975)
  pred_low <- data.frame(xvals = xvals, zvals = yhat_low, upper = yhat_low + se_low, 
                         lower = yhat_low - se_low, y = -1)
  pred_high <- predict(mod, data.frame(x = highdf$x, y = highdf$y), se.fit = T)
  yhat_high <- pred_high$fit
  se_high <- pred_high$se.fit * qnorm(.975)
  pred_high <- data.frame(xvals = xvals, zvals = yhat_high, upper = yhat_high + se_high,
                          lower = yhat_high - se_high, y = 1)
  pred_mod <- rbind(pred_low, pred_high)
  pred_mod$y <- factor(pred_mod$y)
  colnames(pred_mod)[2] <- "zvals"
  levels(pred_mod$y) <- c(paste0("Low ", ylab), paste0("High ", ylab))
  pred_plot <- ggplot(data = pred_mod, aes(x = xvals, group = y)) +
    geom_line(aes(y = zvals, color = y)) + 
    geom_ribbon(alpha = .2, aes(ymin = lower, ymax = upper)) + 
    labs(title = paste0("Effect of ", xlab, " on ", zlab, " Moderated by ", ylab),
         subtitle = paste0("Simple Slopes at +/- 1 SD of Mean-Centered ", ylab),
         y = paste0("Predicted Values of ", zlab),
         x = paste0("Mean Centered Value of ", xlab),
         colour = "") + 
    theme_bw() +
    theme(legend.position = c(.1, .85))
  marg_mod <- margins::cplot(mod, dx = "x", x = "y", what = "effect", data = dat, draw = FALSE)
  marg_plot <- ggplot(data = marg_mod, aes(x = xvals, y = yvals)) +
    geom_line(color = "red") +
    geom_ribbon(alpha = .2, aes(ymin = lower, ymax = upper)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste0("Marginal Effect of ", xlab, " on ", zlab, " as a Function of ", ylab),
         subtitle = paste0("Marginal Effect of ", xlab, " Across Range of Mean-Centered ", ylab),
         y = paste0("Estimated Effect of ", xlab, " on ", zlab),
         x = paste0("Mean Centered Values of ", ylab),
         color = "") +
    theme_bw()
  if(effect){
    plot(marg_plot)
  } else {
    cowplot::plot_grid(pred_plot, marg_plot, ncol = 2, align = "h")
  }
}

##### moderationPlot:
moderationPlot <- function(data, vars = list("Y", "X", "M"), plot = TRUE, 
                           n = 20, mtype = "g", mean = TRUE){
  if(!is(vars, 'list')){vars <- as.list(vars)}
  if(length(vars) != 3){stop("Must provide 3 variable names")}
  if(is.null(names(vars))){names(vars) <- c("Y", "X", "M"); message("Assumed vars order: Y, X, M")}
  varCols <- sapply(vars, function(x) which(colnames(data) == x))
  colnames(data)[varCols] <- names(vars)
  X <- seq(min(data[,"X"]), max(data[,"X"]), length = n)
  form <- "Y ~ X * M"
  if(all(data$M %in% c(0, 1))){mtype <- "c"}
  if(mtype == "c"){
    dat <- list(datZero = data.frame(X = X, M = 0),
                datOne = data.frame(X = X, M = 1))
  } else {
    Mm <- mean(data[,"M"])
    Mlo <- Mm - sd(data[,"M"])
    Mhi <- Mm + sd(data[,"M"])
    dat <- list(datLo = data.frame(X = X, M = Mlo), 
                datM = data.frame(X = X, M = Mm), 
                datHi = data.frame(X = X, M = Mhi))
  }
  if(ncol(data) > 3){
    colnames(data)[-varCols] <- paste0("Z", 1:(ncol(data) - 3))
    zs <- grep("Z", colnames(data))
    for(i in 1:length(zs)){
      zs[i] <- mean(data[, zs[i]])
      dat <- lapply(dat, function(z){
        zz <- cbind(z, zs[i])
        colnames(zz)[i + 2] <- paste0("Z", i)
        return(zz)
      })
      form <- paste0(form, " + Z", i)
    }
  }
  mod <- lm(form, data = data)
  preds <- lapply(dat, function(z){
    zz <- predict(mod, z, se.fit = TRUE)
    zz$se.fit <- zz$se.fit * qnorm(.975)
    zz <- data.frame(X = X, Y = zz$fit, upper = zz$fit + zz$se.fit, lower = zz$fit - zz$se.fit)
    return(zz)
  })
  if(mtype == "c"){
    predMod <- cbind(do.call(rbind, preds), M = c(rep(0, n), rep(1, n)))
  } else {
    predMod <- cbind(do.call(rbind, preds), M = c(rep("Low", n), rep("Mean", n), rep("High", n)))
    if(mean == FALSE){predMod <- predMod[predMod$M != "Mean",]}
  }
  predMod$M <- factor(predMod$M)
  rownames(predMod) <- 1:nrow(predMod)
  if(plot == TRUE){
    require(ggplot2)
    predPlot <- ggplot(data = predMod, aes(x = X, group = M)) +
      geom_line(aes(y = Y, color = M)) +
      geom_ribbon(alpha = .2, aes(ymin = lower, ymax = upper)) +
      labs(title = "Effect of X on Y Moderated by M",
           subtitle = ifelse(mtype == "c", "Simple slopes at levels of M", "Simple slopes at +/- 1 SD of M"),
           y = "Predicted values of Y",
           x = "X",
           color = ifelse(mtype == "c", "M", "")) +
      theme_bw() +
      theme(legend.position = c(.1, .85))
    predPlot
  } else {
    predMod
  }
}

##### mnetPowerSim: power simulator for cross-sectional and idiographic networks
mnetPowerSim <- function(niter = 10, N = 100, p = 5, m = FALSE, m1 = 0, m2 = .1, sparsity = .5,
                         lags = NULL, trueNet = NULL, threshold = NULL, rule = 'OR', avg = TRUE, 
                         maxiter = 100, saveFits = FALSE, saveData = FALSE, intercepts = NULL, 
                         mbinary = FALSE, select = NULL, vargs = list(), type = 'g', 
                         gibbs = TRUE, ordinal = FALSE, mord = FALSE, nLevels = 5, 
                         minOrd = 3, div = 1000, modType = 'none', m1_range = NULL, 
                         m2_range = c(.1, .3), time = TRUE, skewErr = FALSE, 
                         nCores = 1, cluster = 'SOCK', ...){
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
      args1 <- list(Nvar = p, sparsity = sparsity, finalRange = m2_range)
      args1 <- args1[intersect(names(args1), formalArgs('simPcor'))]
      args$b2 <- trueNet2 <- do.call(simPcor, replace(args1, 'sparsity', 1 - m2))
      args1$finalRange <- m1_range
      if(!identical(m, FALSE)){
        if(is.numeric(modType)){modType <- switch(modType + 1, 'none', 'full', 'partial')}
        modType <- match.arg(tolower(modType), c('none', 'full', 'partial'))
      } else {
        modType <- 'none'
      }
      trueNet1 <- do.call(simPcor, args1)
      en1 <- eigen(diag(p) - trueNet1)$values
      while(any(round(en1, 14) == 1)){
        trueNet1 <- do.call(simPcor, args1)
        en1 <- eigen(diag(p) - trueNet1)$values
      }
      if(modType != 'none'){
        while(ifelse(modType == 'full', !all(trueNet1[trueNet2 != 0] == 0), any(trueNet1[trueNet2 != 0] == 0))){
          trueNet1 <- do.call(simPcor, args1)
          en1 <- eigen(diag(p) - trueNet1)$values
          while(any(round(en1, 14) == 1)){
            trueNet1 <- do.call(simPcor, args1)
            en1 <- eigen(diag(p) - trueNet1)$values
          }
        }
      }
      args$b1 <- trueNet1
    } else {
      avg <- FALSE
      m0 <- ifelse(identical(m, FALSE), FALSE, ifelse(mbinary, 'binary', ifelse(mord, 'ordinal', m)))
      out1 <- mlGVARsim(nTime = N[1], nPerson = 1, nNode = p, m = m0, m2 = m2, 
                        keepNets = TRUE, ordinal = ordinal, nLevels = nLevels, 
                        GGMsparsity = sparsity, m1_range = m1_range, 
                        m2_range = m2_range, minOrd = minOrd, 
                        skewErr = skewErr, m1 = m1)
      if(any(eigen(out1$residuals)$values == 1)){
        while(any(eigen(out1$residuals)$values == 1)){
          out1 <- mlGVARsim(nTime = N[1], nPerson = 1, nNode = p, m = m0, m2 = m2, 
                            keepNets = TRUE, ordinal = ordinal, nLevels = nLevels, 
                            GGMsparsity = sparsity, m1_range = m1_range, 
                            m2_range = m2_range, minOrd = minOrd, 
                            skewErr = skewErr, m1 = m1)
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
  args2 <- args[intersect(names(args), setdiff(formalArgs(FUN), '...'))]
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
        parallel::clusterExport(cl, c(FUN, 'fitNetwork', objects), envir = environment())
      }
    } else {
      cl <- nCores
    }
  }
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
          Data <- tryCatch({do.call(match.fun(FUN), args2)}, error = function(e){TRUE})
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
          lags = lags)
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
          Data[[i]] <- tryCatch({do.call(match.fun(FUN), args2)}, error = function(e){TRUE})
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
          lags = lags)
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
          Data <- tryCatch({
            do.call(match.fun(FUN), replace(args2, colnames(vars), as.list(vars0[i])))},
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
          lags = lags)
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
            out0 <- tryCatch({
              do.call(match.fun(FUN), replace(args2, colnames(vars), as.list(vars[j, ])))}, 
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
            lags = lags)
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
  output$args <- append(list(call = as.list(match.call())[-1]), args)
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
plotPower <- function(x, by = NULL, yvar = 'default', yadd = NULL, hline = .8,
                      xlab = 'Number of cases', title = NULL, ...){
  args <- list(...)
  if(is(x, 'list')){x <- x$Results}
  if(identical(yvar, 'default')){yvar <- c('sensitivity', 'specificity', 'correlation')}
  if(any(colnames(x) == 'cor')){colnames(x)[colnames(x) == 'cor'] <- 'correlation'}
  if(!is.null(yadd)){yvar <- c(yvar, yadd)}
  cents <- c('Strength', 'EI')
  xvar0 <- c('N', 'p', 'm2', 'index', 'iter', 'type', 'network')
  yvar0 <- setdiff(c(colnames(x), cents), xvar0)
  yvar <- match.arg(tolower(yvar), unique(tolower(yvar0)), several.ok = TRUE)
  if('fdr' %in% tolower(yvar)){
    x$precision <- 1 - x$precision
    colnames(x)[colnames(x) == 'precision'] <- 'FDR'
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

### ======================================================================== ###
### ======================================================================== ###
##### centTable: mimics centralityTable
centTable <- function(Wmats, scale = TRUE, which.net = "temporal", labels = NULL, 
                      relative = FALSE, weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(Wmats, "mlGVAR"))){
    Wmats <- switch(which.net, between = Wmats$betweenNet, Wmats$fixedNets)}
  if("SURnet" %in% c(names(Wmats), names(attributes(Wmats)))){
    if("SURnet" %in% names(Wmats)){Wmats <- Wmats$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p")[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc"))
    Wmats <- Wmats[[ifelse(which.net == "contemporaneous", "contemporaneous", "temporal")]]
    if(which.net == "pdc"){Wmats <- Wmats$PDC}
  }
  if("adjMat" %in% names(Wmats)){Wmats <- t(Wmats$adjMat)}
  if(any(grepl("lag", dimnames(Wmats)))){dimnames(Wmats) <- lapply(dimnames(Wmats), function(z) gsub("[.]lag1.*|[.]y$", "", z))}
  if(!is.list(Wmats)){Wmats <- list(Wmats)}
  if(any(sapply(Wmats, ncol) == 1)){stop("Not supported for single-node graphs")}
  names(Wmats) <- qgraph:::fixnames(Wmats, "graph ")
  centOut <- lapply(Wmats, centAuto, which.net = which.net, weighted = weighted, signed = signed)
  for(g in seq_along(centOut)){
    if(!is(centOut[[g]], "centrality_auto")){
      names(centOut[[g]]) <- qgraph:::fixnames(centOut[[g]], "type ")
      for(t in seq_along(centOut[[g]])){
        if(!is.null(labels)){
          centOut[[g]][[t]][["node.centrality"]][["node"]] <- labels
        } else if(!is.null(rownames(centOut[[g]][[t]][["node.centrality"]]))){
          centOut[[g]][[t]][["node.centrality"]][["node"]] <- rownames(centOut[[g]][[t]][["node.centrality"]])
        } else {
          centOut[[g]][[t]][["node.centrality"]][["node"]] <- paste("Node", seq_len(nrow(centOut[[g]][[t]][["node.centrality"]])))
        }
        centOut[[g]][[t]]$node.centrality$graph <- names(centOut)[g]
        centOut[[g]][[t]]$node.centrality$type <- names(centOut[[g]])[t]
      }
    } else {
      centOut[[g]]$node.centrality$graph <- names(centOut)[g]
      if(!is.null(labels)){
        centOut[[g]][["node.centrality"]][["node"]] <- labels
      } else if(!is.null(rownames(centOut[[g]][["node.centrality"]]))){
        centOut[[g]][["node.centrality"]][["node"]] <- rownames(centOut[[g]][["node.centrality"]])
      } else {
        centOut[[g]][["node.centrality"]][["node"]] <- paste("Node", seq_len(nrow(centOut[[g]][["node.centrality"]])))
      }
    }
  }
  isList <- sapply(centOut, function(x) !"centrality_auto" %in% class(x))
  if(any(isList)){
    for(l in which(isList)){centOut <- c(centOut, centOut[[l]])}
    centOut <- centOut[-which(isList)]
  }
  for(i in seq_along(centOut)){
    if(relative | scale){
      if(relative & scale){warning("Using 'relative' and 'scale' together is not recommended")}
      for(j in which(sapply(centOut[[i]][["node.centrality"]], mode) == "numeric")){
        if(scale){centOut[[i]][["node.centrality"]][, j] <- qgraph:::scale2(centOut[[i]][["node.centrality"]][, j])}
        if(relative){
          mx <- max(abs(centOut[[i]][["node.centrality"]][, j]), na.rm = TRUE)
          if(mx != 0){centOut[[i]][["node.centrality"]][, j] <- centOut[[i]][["node.centrality"]][, j]/mx}
        }
        attributes(centOut[[i]][["node.centrality"]][, j]) <- NULL
      }
    }
  }
  wideCent <- plyr::rbind.fill(lapply(centOut, "[[", "node.centrality"))
  if(is.null(wideCent$type)){wideCent$type <- NA}
  longCent <- reshape2::melt(wideCent, variable.name = "measure", id.var = c("graph", "type", "node"))
  if(any(is.nan(longCent$value))){warning("NaN detected in centrality measures. Try relative = FALSE")}
  return(longCent)
}

##### centAuto: mimics centrality_auto
centAuto <- function(x, which.net = "temporal", weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(x, "mlGVAR"))){
    x <- switch(which.net, between = x$betweenNet, x$fixedNets)}
  if("SURnet" %in% c(names(x), names(attributes(x)))){
    if("SURnet" %in% names(x)){x <- x$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p")[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc"))
    x <- x[[ifelse(which.net == "contemporaneous", "contemporaneous", "temporal")]]
    if(which.net == "pdc"){x <- x$PDC}
  }
  if("adjMat" %in% names(x)){x <- t(x$adjMat)}
  if(any(grepl("lag", dimnames(x)))){dimnames(x) <- lapply(dimnames(x), function(z) gsub("[.]lag1.*|[.]y$", "", z))}
  if(is.list(x)){return(lapply(x, centAuto, which.net = which.net, weighted = weighted, signed = signed))}
  if(!weighted){x <- sign(x)}
  if(!signed){x <- abs(x)}
  if(!is.matrix(x)){stop("The input network must be an adjacency or weights matrix")}
  diag(x) <- 0
  directed.gr <- ifelse(isSymmetric.matrix(object = x, tol = 0.000000000001), FALSE, TRUE)
  weighted.gr <- ifelse(all(qgraph::mat2vec(x) %in% c(0, 1)), FALSE, TRUE)
  net_qg <- qgraph::qgraph(x, diag = FALSE, labels = colnames(x), DoNotPlot = TRUE, minimum = 0)
  centr <- qgraph::centrality(net_qg)
  if(directed.gr & !weighted.gr){
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness, Closeness = centr$Closeness, 
      InDegree = centr$InDegree, OutDegree = centr$OutDegree, 
      OutExpectedInfluence = centr$OutExpectedInfluence, 
      InExpectedInfluence = centr$InExpectedInfluence))
  }
  if(directed.gr & weighted.gr){
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness, Closeness = centr$Closeness, 
      InStrength = centr$InDegree, OutStrength = centr$OutDegree, 
      OutExpectedInfluence = centr$OutExpectedInfluence, 
      InExpectedInfluence = centr$InExpectedInfluence))
  }
  if(!directed.gr & !weighted.gr){
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness/2, Closeness = centr$Closeness, 
      Degree = centr$OutDegree, ExpectedInfluence = centr$OutExpectedInfluence))
  }
  if(!directed.gr & weighted.gr){ 
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness/2, Closeness = centr$Closeness, 
      Strength = centr$OutDegree, ExpectedInfluence = centr$OutExpectedInfluence))
  }
  row.names(centr1) <- colnames(x)
  log <- capture.output({
    graph <- igraph::graph.adjacency(
      adjmatrix = 1 * (x != 0), 
      mode = ifelse(directed.gr, "directed", "undirected"))
    comps <- igraph::components(graph)
    largcomp <- comps$membership == which.max(comps$csize)
  })
  if(sum(largcomp) < ncol(x) & sum(largcomp) > 1){
    x2 <- x[largcomp, largcomp]
    clos <- qgraph::centrality(qgraph::qgraph(
      x2, diag = FALSE, labels = colnames(x)[largcomp], 
      DoNotPlot = TRUE, minimum = 0))$Closeness
    centr1$Closeness[largcomp] <- clos
    centr1$Closeness[!largcomp] <- NA
  }
  net_ig_abs <- igraph::graph.adjacency(
    adjmatrix = abs(1/x), mode = ifelse(directed.gr, "directed", "undirected"), 
    weighted = ifelse(weighted.gr, list(TRUE), list(NULL))[[1]], diag = FALSE)
  edgebet <- igraph::edge.betweenness(graph = net_ig_abs, directed = directed.gr)
  el <- data.frame(igraph::get.edgelist(graph = net_ig_abs), stringsAsFactors = FALSE)
  edgebet <- merge(el, edgebet, by = 0)
  edgebet$Row.names <- NULL
  names(edgebet) <- c("from", "to", "edgebetweenness")
  edgebet <- edgebet[order(edgebet$edgebetweenness, decreasing = TRUE), ]
  ShortestPathLengths <- centr$ShortestPathLengths
  rownames(ShortestPathLengths) <- colnames(ShortestPathLengths) <- colnames(x)
  Res <- list(node.centrality = centr1, edge.betweenness.centrality = edgebet, 
              ShortestPathLengths = ShortestPathLengths)
  class(Res) <- c("list", "centrality_auto")
  return(Res)
}

##### centPlot: mimics centralityPlot
centPlot <- function(Wmats, scale = c("z-scores", "raw", "raw0", "relative"), 
                     which.net = "temporal", include = "all", labels = NULL, 
                     orderBy = NULL, decreasing = FALSE, plot = TRUE, 
                     verbose = TRUE, weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(Wmats, "mlGVAR"))){
    Wmats <- switch(which.net, between = Wmats$betweenNet, Wmats$fixedNets)}
  if(is.logical(scale)){scale <- ifelse(scale, "z-scores", "raw")}
  invisible(suppressMessages(require(ggplot2)))
  measure <- value <- node <- type <- NULL
  scale <- match.arg(scale)
  include0 <- c("Degree", "Strength", "OutDegree", "InDegree", "OutStrength", 
                "InStrength", "Closeness", "Betweenness", "ExpectedInfluence", 
                "OutExpectedInfluence", "InExpectedInfluence")
  if(all(tolower(include) == "all")){
    include <- include0
  } else if(isTRUE(attr(Wmats, "SURnet")) & which.net != "contemporaneous"){
    include0 <- include0[!grepl("Degree|^S|^E", include0)]
    include <- include0[grep(paste(tolower(include), collapse = "|"), tolower(include0))]
  } else if(which.net %in% c("between", "contemporaneous")){
    include0 <- include0[!grepl("Degree|^Out|^In", include0)]
    include <- include0[grep(paste(tolower(include), collapse = "|"), tolower(include0))]
  }
  include <- match.arg(include, c("Degree", "Strength", "OutDegree", "InDegree", "OutStrength", 
                                  "InStrength", "Closeness", "Betweenness", "ExpectedInfluence", 
                                  "OutExpectedInfluence", "InExpectedInfluence"), several.ok = TRUE)
  if(scale == "z-scores" & verbose & plot){message("Note: z-scores are shown on x-axis.")}
  if(scale == "relative" & verbose & plot){message("Note: relative centrality indices are shown on x-axis.")}
  Long <- centTable(Wmats = Wmats, scale = (scale == "z-scores"), labels = labels, 
                    which.net = which.net, relative = (scale == "relative"), 
                    weighted = weighted, signed = signed)
  Long <- subset(Long, measure %in% include)
  Long$measure <- factor(Long$measure)
  if(ifelse(is.null(orderBy), FALSE, ifelse(orderBy == "default", TRUE, FALSE))){
    nodeLevels <- unique(gtools::mixedsort(
      as.character(Long$node), decreasing = decreasing))
  } else if(!is.null(orderBy)){
    nodeLevels <- names(sort(tapply(
      Long$value[Long$measure == orderBy], 
      Long$node[Long$measure == orderBy], mean), 
      decreasing = decreasing))
  } else {
    nodeLevels <- rev(unique(as.character(Long$node)))
  }
  Long$node <- factor(as.character(Long$node), levels = nodeLevels)
  Long <- Long[gtools::mixedorder(Long$node), ]
  if(length(unique(Long$type)) > 1){
    g <- ggplot(Long, aes(x = value, y = node, group = type, colour = type))
  } else {
    g <- ggplot(Long, aes(x = value, y = node, group = type))
  }
  g <- g + geom_path() + xlab("") + ylab("") + geom_point() + theme_bw()
  if(length(unique(Long$graph)) > 1){
    g <- g + facet_grid(graph ~ measure, scales = "free")
  } else {
    g <- g + facet_grid(~measure, scales = "free")
  }
  if(scale == "raw0"){g <- g + xlim(0, NA)}
  if(plot){plot(g)} else {invisible(g)}
}


### ======================================================================== ###
### ======================================================================== ###
##### clustTable: mimics clusteringTable
clustTable <- function(Wmats, scale = TRUE, labels = NULL, 
                       relative = FALSE, signed = TRUE){
  if("SURnet" %in% c(names(Wmats), names(attributes(Wmats)))){
    if("SURnet" %in% names(Wmats)){Wmats <- Wmats$SURnet}
    Wmats <- Wmats$contemporaneous$adjMat
  } else if("adjMat" %in% names(Wmats)){
    Wmats <- Wmats$adjMat
  }
  if(any(grepl("lag", dimnames(Wmats)))){
    dimnames(Wmats) <- lapply(dimnames(Wmats), function(z) gsub("[.]lag1.*|[.]y$", "", z))
  }
  if(!is.list(Wmats)){Wmats <- list(Wmats)}
  if(any(sapply(Wmats, ncol) == 1)){stop("Not supported for single-node graphs")}
  syms <- sapply(Wmats, isSymmetric)
  if(any(!syms)){
    if(all(!syms)){stop("No symmetrical graphs detected")}
    warning(paste0(sum(!syms), " Nonsymmetrical graph", ifelse(sum(!syms) > 1, "s ", " "), "removed"))
    Wmats <- Wmats[-which(!syms)]
  }
  names(Wmats) <- qgraph:::fixnames(Wmats, "graph ")
  clustOut <- lapply(Wmats, clustAuto)
  for(g in seq_along(clustOut)){
    if(!is(clustOut[[g]], "clustcoef_auto")){
      names(clustOut[[g]]) <- qgraph:::fixnames(clustOut[[g]], "type ")
      for(t in seq_along(clustOut[[g]])){
        if(!is.null(labels)){
          clustOut[[g]][[t]][["node"]] <- labels
        } else if(!is.null(rownames(clustOut[[g]][[t]]))){
          clustOut[[g]][[t]][["node"]] <- rownames(clustOut[[g]][[t]])
        } else {
          clustOut[[g]][[t]][["node"]] <- paste("Node", seq_len(nrow(clustOut[[g]][[t]])))
        }
        clustOut[[g]][[t]]$graph <- names(clustOut)[g]
        clustOut[[g]][[t]]$type <- names(clustOut[[g]])[t]
      }
    } else {
      clustOut[[g]]$graph <- names(clustOut)[g]
      if(!is.null(labels)){
        clustOut[[g]][["node"]] <- labels
      } else if(!is.null(rownames(clustOut[[g]]))){
        clustOut[[g]][["node"]] <- rownames(clustOut[[g]])
      } else {
        clustOut[[g]][["node"]] <- paste("Node", seq_len(nrow(clustOut[[g]])))
      }
    }
  }
  isList <- sapply(clustOut, function(x) !"clustcoef_auto" %in% class(x))
  if(any(isList)){
    for(l in which(isList)){clustOut <- c(clustOut, clustOut[[l]])}
    clustOut <- clustOut[-which(isList)]
  }
  for(i in seq_along(clustOut)){
    if(any(grepl("signed_", names(clustOut[[i]])))){
      clustOut[[i]] <- clustOut[[i]][, sapply(clustOut[[i]], mode) != "numeric" | grepl("signed_", names(clustOut[[i]])) == signed]
      names(clustOut[[i]]) <- gsub("signed_", "", names(clustOut[[i]]))
    }
    names(clustOut[[i]]) <- gsub("clust", "", names(clustOut[[i]]))
    if(relative | scale){
      if(relative & scale){warning("Using 'relative' and 'scale' together is not recommended")}
      for(j in which(sapply(clustOut[[i]], mode) == "numeric")){
        if(scale){clustOut[[i]][, j] <- qgraph:::scale2(clustOut[[i]][, j])}
        if(relative){
          mx <- max(abs(clustOut[[i]][, j]), na.rm = TRUE)
          if(mx != 0){clustOut[[i]][, j] <- clustOut[[i]][, j]/mx}
        }
        attributes(clustOut[[i]][, j]) <- NULL
      }
    }
  }
  WideCent <- plyr::rbind.fill(clustOut)
  if(is.null(WideCent$type)){WideCent$type <- NA}
  LongCent <- reshape2::melt(WideCent, variable.name = "measure", id.var = c("graph", "type", "node"))
  return(LongCent)
}

##### clustAuto: mimics clustcoef_auto
clustAuto <- function(x, thresholdWS = 0, thresholdON = 0){
  if("SURnet" %in% c(names(x), names(attributes(x)))){
    if("SURnet" %in% names(x)){x <- x$SURnet}
    x <- x$contemporaneous$adjMat
  } else if("adjMat" %in% names(x)){
    x <- x$adjMat
  }
  if(any(grepl("lag", dimnames(x)))){
    dimnames(x) <- lapply(dimnames(x), function(z) gsub("[.]lag1.*|[.]y$", "", z))
  }
  if(is.list(x)){
    return(lapply(x, clustAuto, thresholdWS = thresholdWS, thresholdON = thresholdWS))
  }
  dim = dim(x)
  if(is.null(dim) || length(dim) != 2){stop("adjacency is not two-dimensional")}
  if(!is.numeric(x)){stop("adjacency is not numeric")}
  if(dim[1] != dim[2]){stop("adjacency is not square")}
  if(max(abs(x - t(x)), na.rm = TRUE) > 0.000000000001){stop("adjacency is not symmetric")}
  if(min(x, na.rm = TRUE) < -1 || max(x, na.rm = TRUE) > 1){x <- x/max(abs(x))}
  weighted.gr <- ifelse(all(abs(x) %in% c(0, 1)), FALSE, TRUE)
  signed.gr <- ifelse(all(x >= 0), FALSE, TRUE)
  net_ig <- igraph::graph.adjacency(
    adjmatrix = abs(x), mode = "undirected", 
    weighted = ifelse(weighted.gr, list(TRUE), list(NULL))[[1]], diag = FALSE)
  cb <- igraph::transitivity(net_ig, type = "barrat", isolates = "zero")
  cw <- qgraph:::clustWS(x, thresholdWS)
  cz <- qgraph:::clustZhang(x)
  co <- qgraph:::clustOnnela(x, thresholdON)
  if(!signed.gr & !weighted.gr){output <- cbind(clustWS = cw[, 1])}
  if(!signed.gr & weighted.gr){
    output <- cbind(clustWS = cw[, 1], clustZhang = cz[, 1], 
                    clustOnnela = co[, 1], clustBarrat = cb)
  }
  if(signed.gr & !weighted.gr){
    output <- cbind(clustWS = cw[, 1], signed_clustWS = cw[, 2])
  }
  if(signed.gr & weighted.gr){
    output <- cbind(clustWS = cw[, 1], signed_clustWS = cw[, 2], 
                    clustZhang = cz[, 1], signed_clustZhang = cz[, 2], 
                    clustOnnela = co[, 1], signed_clustOnnela = co[, 2], 
                    clustBarrat = cb)
  }
  output[is.na(output)] <- 0
  Res <- data.frame(output)
  class(Res) <- c("data.frame", "clustcoef_auto")
  rownames(Res) <- colnames(x)
  return(Res)
}

##### clustPlot: mimics clusteringPlot
clustPlot <- function(Wmats, scale = c("z-scores", "raw", "raw0", "relative"), 
                      include = "all", labels = NULL, orderBy = NULL, 
                      decreasing = FALSE, plot = TRUE, signed = TRUE, 
                      verbose = TRUE){
  if(is.logical(scale)){scale <- ifelse(scale, "z-scores", "raw")}
  invisible(suppressMessages(require(ggplot2)))
  measure <- value <- node <- type <- NULL
  scale <- match.arg(scale)
  if(scale == "z-scores" & verbose & plot){message("Note: z-scores are shown on x-axis.")}
  if(scale == "relative" & verbose & plot){message("Note: relative centrality indices are shown on x-axis.")}
  Long <- clustTable(Wmats = Wmats, scale = (scale == "z-scores"), labels = labels, 
                     relative = (scale == "relative"), signed = signed)
  Long$value[!is.finite(Long$value)] <- 0
  if(all(include == "all")){include <- c("WS", "Zhang", "Onnela", "Barrat")}
  include <- match.arg(include, c("WS", "Zhang", "Onnela", "Barrat"), several.ok = TRUE)
  Long <- subset(Long, measure %in% include)
  Long$measure <- factor(Long$measure)
  if(ifelse(is.null(orderBy), FALSE, ifelse(orderBy == "default", TRUE, FALSE))){
    nodeLevels <- unique(gtools::mixedsort(
      as.character(Long$node), decreasing = decreasing))
  } else if(!is.null(orderBy)){
    nodeLevels <- names(sort(tapply(
      Long$value[Long$measure == orderBy], 
      Long$node[Long$measure == orderBy], mean), 
      decreasing = decreasing))
  } else {
    nodeLevels <- rev(unique(as.character(Long$node)))
  }
  Long$node <- factor(as.character(Long$node), levels = nodeLevels)
  Long <- Long[gtools::mixedorder(Long$node), ]
  if(length(unique(Long$type)) > 1){
    g <- ggplot(Long, aes(x = value, y = node, group = type, colour = type))
  } else {
    g <- ggplot(Long, aes(x = value, y = node, group = type))
  }
  g <- g + geom_path() + xlab("") + ylab("") + geom_point() + theme_bw()
  if(length(unique(Long$graph)) > 1){
    g <- g + facet_grid(graph ~ measure, scales = "free")
  } else {
    g <- g + facet_grid(~measure, scales = "free")
  }
  if(scale == "raw0"){g <- g + xlim(0, NA)}
  if(plot){plot(g)} else {invisible(g)}
}


### ======================================================================== ###
### ======================================================================== ###
##### plotCentrality: plot centrality (and clustering) for multiple networks
plotCentrality <- function(Wmats, which.net = "temporal", scale = TRUE, 
                           labels = NULL, plot = TRUE, centrality = "all", 
                           clustering = "Zhang"){
  if(any(c("ggm", "SURnet", "mlGVAR") %in% names(attributes(Wmats)))){Wmats <- list(net1 = Wmats)}
  if(all(sapply(Wmats, function(z) isTRUE(attr(z, "mlGVAR"))))){
    Wmats <- lapply(Wmats, function(z) switch(
      which.net, between = z$betweenNet, z$fixedNets))
  }
  if(any(grepl("ggm", lapply(Wmats, function(z) names(attributes(z)))))){which.net <- "contemporaneous"}
  if(length(unique(lapply(Wmats, checkInclude))) != 1){stop("All networks must be of the same type")}
  if(is.null(names(Wmats))){names(Wmats) <- paste0("net", seq_along(Wmats))}
  which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc"))
  c0 <- c01 <- do.call(rbind, lapply(seq_along(Wmats), function(z){
    cbind.data.frame(
      centTable(Wmats[[z]], scale = scale, which.net = which.net, 
                labels = labels), group = names(Wmats)[z])
  }))
  if(all(centrality != "all")){
    include0 <- checkInclude(Wmats[[1]], which.net = which.net)
    include0 <- include0[!grepl(ifelse(
      which.net != "contemporaneous", "Degree|^S|^E", 
      "Degree|^Out|^In"), include0)]
    centrality <- include0[grep(paste(tolower(
      centrality), collapse = "|"), tolower(include0))]
    c0 <- c01 <- subset(c0, measure %in% centrality)
  }
  if(which.net == "contemporaneous" & clustering != FALSE){
    c1 <- do.call(rbind, lapply(seq_along(Wmats), function(z){
      z1 <- clustTable(Wmats[[z]], scale = scale, labels = labels)
      z1 <- z1[z1$measure == ifelse(is.logical(clustering), "Zhang", clustering), ]
      z1$measure <- "Clust. coef."
      z1$node <- as.character(z1$node)
      rownames(z1) <- 1:nrow(z1)
      return(cbind.data.frame(z1, group = names(Wmats)[z]))
    }))
    c01 <- rbind(c0, c1)
  }
  c01 <- c01[order(c01$node), ]
  c01 <- c01[order(c01$group), ]
  rownames(c01) <- 1:nrow(c01)
  c01$node <- stringr::str_sub(c01$node, 1, 6)
  if(!plot){
    list2env(list(c0 = c0, c1 = c1, c01 = c01), .GlobalEnv)
  } else {
    invisible(suppressMessages(require(ggplot2)))
    g1 <- ggplot(c01, aes(x = value, y = node, group = group, color = group, shape = group)) +
      geom_path(alpha = 1, size = 1) + geom_point(size = 2) + xlab("") + ylab("") + theme_bw() + 
      facet_grid(. ~ measure, scales = "free") + scale_x_continuous(breaks = c(-1, 0, 1)) + 
      theme(axis.line.x = element_line(colour = "black"),
            axis.ticks.x = element_line(colour = "black"),
            axis.ticks.y = element_line(colour = "white", size = 0),
            axis.text.y = element_text(colour = "black"),
            axis.text.x = element_text(angle = 45, colour = "black"))
    g1
  }
}

##### checkInclude
checkInclude <- function(x, which.net = "temporal"){
  if("SURnet" %in% c(names(x), names(attributes(x)))){
    if("SURnet" %in% names(x)){x <- x$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p")[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc"))
    x <- x[[ifelse(which.net == "contemporaneous", "contemporaneous", "temporal")]]
    if(which.net == "pdc"){x <- x$PDC}
  }
  if("adjMat" %in% names(x)){x <- t(x$adjMat)}
  directed <- !isTRUE(all.equal(x, t(x), check.attributes = FALSE))
  weighted <- any(!x %in% c(0, 1))
  include <- c("Degree", "Strength", "OutDegree", "InDegree", "OutStrength", 
               "InStrength", "Closeness", "Betweenness", "ExpectedInfluence", 
               "OutExpectedInfluence", "InExpectedInfluence")
  if(directed){
    include <- c("OutStrength", "InStrength", "Closeness", "Betweenness",
                 "OutExpectedInfluence", "InExpectedInfluence")
  } else {
    include <- c("Strength", "Closeness", "Betweenness", "ExpectedInfluence")
  }
  if(!weighted){include <- gsub("Strength", "Degree", include)}
  return(include)
}


