if(!exists("pkgs")){pkgs <- c("glmnet", "mgm", "qgraph", "grpregOverlap", "systemfit", "leaps", "hierNet", "glinternet")}
if(!"none" %in% pkgs){invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))}
if(exists("clear")){if(clear == TRUE){rm(list = ls()); clear <- TRUE}} else {rm(pkgs)}
if(".modnets" %in% search()){detach(".modnets")}
source("~/C: Desktop/COMPS/METHODS/CODE/modnets/ggm.R")
source("~/C: Desktop/COMPS/METHODS/CODE/modnets/stabilityPaths.R")

################################################################################
# ~~~~~~~ MODNETS: Estimating and Analyzing Moderated Temporal Networks ~~~~~~ #
################################################################################
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
  if(d == TRUE){
    whatData <- c("sur1", "sur2", "sur3", "mvar", "mgm", "gauss1", "gauss2", 
                  "autism", "large", "symptom", "short", "resting", "big5", 
                  "m3", "m5", "mod", "PTSD", "fried", "fried2", "wichers", 
                  "new", "dep1", "dep2", "fried3", "bfi1", "bfi2", "bfi3", "obama")
    whatData <- sort(whatData)
    return(data.frame(whatData))
  }
  if(!is.null(dat)){
    if("resample" %in% c(dat, moderators)){
      poop = list(niter = 3, m = NULL, criterion = "AIC", sampMethod = "split", 
                  gamma = .25, nfolds = 10, nlam = 50, which.lam = "min", rule = "OR",
                  varMethod = "glinternet", threshold = FALSE, bonf = TRUE, alpha = .05,
                  split = .5, center = TRUE, varSeed = NULL, seeds = NULL, verbose = TRUE,
                  exogenous = TRUE, binary = NULL, type = "gaussian", saveMods = FALSE)
      list2env(poop, envir = .GlobalEnv)
      if(dat == "resample"){return(message("settings loaded for resample2"))}
    }
    if(grepl("sur", dat)){
      if(dat == "sur1"){
        S <- matrix(c(1, .4, .3, .4, 1, .6, .3, .6, 1), 3, 3)
        beta <- matrix(c(.01, .014, 0, .005, .07, .3, .01, .09, .1, .23, .3, .001), 
                       ncol = 4, nrow = 3, byrow = TRUE)
        beta2 <- matrix(c(.01, .05, .1, .1, .01, .005), ncol = 2, byrow = T)
        poop <- list(B = cbind(beta, beta2), S = S)
        list2env(poop, envir = .GlobalEnv)
        return(message("B and S in .GlobalEnv"))
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
        list2env(poop, envir = .GlobalEnv)
        return(message("B and S in .GlobalEnv"))
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
        list2env(poop, envir = .GlobalEnv)
        return(message("B and S in .GlobalEnv"))
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
      if(dat == "gauss1"){data <- readRDS("datG.RDS")$data}
      if(dat == "gauss2"){data <- readRDS("dataG.RDS")$data}
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
      data <- Fried2015$data
      type <- Fried2015$type
      level <- Fried2015$level
    }
    if(dat == "fried2"){
      data <- readRDS("Fried2015_nD.RDS")$data
      colnames(data) <- readRDS("Fried2015_nD.RDS")$names
      type <- c(rep("g", 11), "c")
      level <- c(rep(1, 11), 2)
    }
    if(dat == "fried3"){
      data <- readRDS("Fried2015_nD.RDS")$data
      colnames(data) <- readRDS("Fried2015_nD.RDS")$names
      mod <- list(loss = data[, 12])
      dat <- list(Y = data[, 1:11], X = data)
      data <- data[, 1:11]
      poop <- list(data = data.frame(data, mod))
      list2env(poop, envir = .GlobalEnv)
      return(message("data in .GlobalEnv"))
    }
    if(dat == "obama"){
      path <- "~/C: Documents/TREVOR/Computing/attitudes/"
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
      data <- readRDS("Wicherts2016_Mood.RDS")$data_mood
      colnames(data) <- readRDS("Wicherts2016_Mood.RDS")$labels
      type <- rep("g", 9)
      level <- rep(1, 9)
    }
    if(dat == "new"){
      data <- readRDS("newData.RDS")
      type <- data$call$type
      level <- data$call$level
      data <- data$data
    }
    if(grepl("bfi", dat)){
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
    poop <- list(data = data, type = type, level = level, lambda = NULL, 
                 lags = lags, rule = "AND", moderators = moderators, std = TRUE,
                 which.lam = "lambda.min", gamma = 0.25, seed = 666, alpha = 1,
                 threeWay = TRUE, folds = 10, family = NULL, threshold = TRUE,
                 measure = "deviance", saveData = TRUE, group = NULL, scale = FALSE,
                 adaMethod = "ridge", adaGam = NULL, grPenalty = "grLasso", time = FALSE, 
                 covariates = NULL, center = TRUE, verbose = TRUE, exogenous = NULL)
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

################################################################################
################################################################################
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
      betas <- lapply(1:length(allCoefs), function(z) as.matrix(c("(Intercept)" = fit$betahat[[z + 1]][1], allCoefs[[z]]), ncol = 1))
      preds <- lapply(1:length(betas), function(z) cbind(1, X) %*% betas[[z]])
      s2 <- sapply(1:length(betas), function(z) sum((y - preds[[z]])^2)/length(y))
      LL_models <- sapply(1:length(preds), function(z) sum(dnorm(y, mean = preds[[z]], sd = sqrt(s2[z]), log = TRUE)))
      deviance <- sapply(1:length(preds), function(z) sum((y - preds[[z]][, 1])^2))
      ic_lambda <- -2 * LL_models + n_neighbors * ifelse(
        "AIC" %in% lambda, 2, log(n)) + ifelse(
          "EBIC" %in% lambda, 2 * gamma * n_neighbors * log(p), 0)
      n_neighbors <- n_neighbors[which.min(ic_lambda)]
      betas <- betas[[which.min(ic_lambda)]]
      LL_model <- LL_models[which.min(ic_lambda)]
      deviance <- deviance[which.min(ic_lambda)]
      lambda_min <- fit$lambda[which.min(ic_lambda) + 1]
      modFitIndex <- ic_lambda[which.min(ic_lambda)]
      fit <- glinternet(X = X2, Y = y, numLevels = rep(1, ncol(X2)), 
                        interactionCandidates = m, lambda = lambda_min,
                        family = fam)
      fit <- list(fit = fit)
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
                        which.lam = lam, measure = measure, lambda = lambda, std = std)
      }
      betas <- coef(fit, s = lam)
      deviance <- (1 - fit$glmnet.fit$dev.ratio) * fit$glmnet.fit$nulldev
      deviance <- deviance[which(fit$lambda == fit[[lam]])]
    } else {
      if(!is.null(adaGam)){warning("Adaptive lasso not implemented for group lasso")}
      if(fam == "multinomial"){fam <- "binomial"}
      fit <- cv.grpregOverlap(X = X, y = y, family = fam, group = group, dfmax = length(group)*3,
                              nfolds = folds, alpha = alpha, penalty = grPenalty)
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
      criterion == "EBIC", 2 * gamma * n_neighbors * log(p), 0)
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

##### makeInts: create all interaction terms from a design matrix
makeInts <- function(x, mat = TRUE){
  k <- ncol(x)
  stopifnot(k > 1)
  k2 <- (k * (k - 1))/2
  x2 <- data.frame(matrix(NA, ncol = k2, nrow = nrow(x)))
  ints <- t(combn(x = c(1:k), m = 2))
  for(i in 1:k2){x2[, i] <- x[, ints[i, 1]] * x[, ints[i, 2]]}
  colnames(x2) <- paste0(colnames(x)[ints[,1]], ":", colnames(x)[ints[,2]])
  if(mat){x2 <- as.matrix(x2)}
  x2
}

################################################################################
################################################################################
##### combineMods: Feeds model results to 'combine' function; creates n lists for n lags
combineMods <- function(data, mods, type, moderators = NULL, threshold = TRUE,
                        rule = "AND", lags = NULL, scale = TRUE){
  if(length(type) == 1){
    typetype <- type
    type <- rep("c", ncol(data))
    attr(type, "family") <- typetype
  }
  if(is.null(lags)){
    output <- combine(data = data, mods = mods, type = type, threshold = threshold,
                      moderators = moderators, rule = rule, scale = scale)
  } else {
    if(length(lags) == 1){
      output <- combine(data = data, mods = mods, type = type, threshold = threshold,
                        moderators = moderators, lags = lags, scale = scale)
    } else {
      n_lags <- length(lags)
      lag_mods <- list(); output <- list()
      for(i in 1:n_lags){
        lag_mods[[i]] <- lagMods(mods = mods, type = type, lag = lags[i])
        output[[i]] <- combine(data = data, mods = lag_mods[[i]], type = type, threshold = threshold,
                               moderators = moderators, lags = lags[i], scale = scale)
      }
    }
  }
  output
}

##### combine: Main function used to extract and aggregate nodewise coefs for MGMs and mVARs
combine <- function(data, mods, type, moderators = NULL, threshold = TRUE,
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

################################################################################
################################################################################
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

################################################################################
################################################################################
##### fitNetwork: Top-level function that ties everything together
fitNetwork <- function(data, moderators = NULL, type = "gaussian", seed = NULL,
                       lags = NULL, lambda = NULL, folds = 10, gamma = 0.25,
                       which.lam = 'lambda.min', rule = "OR", threshold = FALSE,
                       scale = FALSE, std = TRUE, group = NULL, grPenalty = "gel",
                       adaGam = NULL, adaMethod = "ridge", measure = "deviance",
                       alpha = 1, saveData = TRUE, center = TRUE, covariates = NULL, 
                       verbose = TRUE, exogenous = TRUE, binary = NULL, ...){
  t1 <- Sys.time()
  data <- data.frame(data)
  if(!is.null(lags)){if(lags == TRUE){lags <- 1}}
  if(!is.null(moderators)){if(all(moderators == 0)){moderators <- NULL}}
  if(is.null(lambda) & ncol(data) >= nrow(data)){lambda <- "EBIC"}
  if(is.null(lambda)){
    if(grepl("min", which.lam)){which.lam <- "min"} else {which.lam <- "1se"}
    output <- list()
    output$call <- list(type = type, moderators = moderators, lags = lags,
                        which.lam = which.lam, rule = rule, threshold = threshold,
                        center = center, scale = scale)
    args <- list(...)
    if(length(args) != 0){output$call <- append(output$call, args)}
    if(is.null(type)){
      output$call[c("which.lam")] <- NULL
    } else if(class(type) == "list"){
      output$call$type <- lapply(type, function(z) unlist(z[[ifelse(which.lam == "min", 1, 2)]]))
      if(attr(type, "criterion") != "CV"){output$call$which.lam <- attr(type, "criterion")}
      if(attr(type, "criterion") == "EBIC"){output$call$gamma <- attr(type, "gamma")}
    } else {
      if(all(type == "gaussian")){type <- rep("g", ncol(data))}
      if(all(type == "binomial")){
        type <- rep("c", ncol(data))
        output$call$center <- output$call$scale <- center <- scale <- FALSE
      }
      type <- ifelse(length(type) == 1, ifelse(
        type %in% c("g", "gaussian"), list(rep("g", ncol(data))), 
        list(rep("c", ncol(data)))), list(type))[[1]]
      output$call$type <- type
    }
    if(is.null(lags)){
      output$call["lags"] <- NULL
      if(!is.null(covariates)){
        if(class(covariates) == "list"){covs <- names(covariates)}
        if(class(covariates) %in% c("numeric", "integer")){covs <- colnames(data)[covariates]}
        output$call$covariates <- covs
      }
      if(!is.null(moderators)){
        output$call$moderators <- ifelse(is.list(moderators), list(names(moderators)), list(colnames(data)[moderators]))[[1]]
        if("pcor" %in% names(output$call)){output$call$pcor <- NULL}
      } else {
        if("pcor" %in% names(output$call)){
          if(threshold == FALSE){output$call$pcor <- TRUE}
          output$call$rule <- NULL
        }
      }
      mods0 <- nodewise(data = data, mods = moderators, varMods = type, 
                        lambda = which.lam, center = center, scale = scale,
                        covariates = covariates, exogenous = exogenous)
      out <- modNet(models = mods0, threshold = threshold, rule = rule, ...)
      output <- append(output, out)
      output$mods0 <- mods0
      if("moderator" %in% names(attributes(out))){
        attr(output, "moderator") <- attr(out, "moderator")
        output$call$exogenous <- attr(mods0$models, "exogenous")
      }
      if("mval" %in% names(attributes(out))){attr(output, "mval") <- attr(out, "mval")}
      attributes(output)$ggm <- TRUE
      if(verbose){print(Sys.time() - t1)}
      return(output)
    } else {
      output$call[c("rule", "threshold")] <- NULL
      if(!is.null(moderators)){
        if(!is.null(moderators)){if(length(moderators) > 1){stop("Can only specify one moderator")}}
        if(moderators != 1){
          xn <- colnames(data)
          data <- cbind.data.frame(data[, moderators], data[, -moderators])
          if(!is.null(xn)){colnames(data) <- c(xn[moderators], xn[-moderators])}
          output$call$moderators <- moderators <- 1
        }
      }
      xn <- colnames(data)
      dat <- datMat(data = data, type = "g", m = moderators, lags = 1, scale = scale, center = center)
      if(!is.null(moderators)){
        colnames(dat$full)[(ncol(data) + 1):ncol(dat$full)] <- colnames(dat$X) <- c(xn, paste0(xn[1], ":", xn[-1]))
      }
      if(class(type) == "character"){type <- NULL}
      fit <- SURfit(dat = dat, varMods = type, m = moderators, mod = which.lam)
      net <- SURnet(fit = fit, dat = dat)
      output$SURnet <- net
      output$SURfit <- fit
      output$dat <- dat
      if(verbose){print(Sys.time() - t1)}
      return(output)
    }
  }
  if(!is.null(moderators)){
    if(class(moderators) == "list"){
      data <- data.frame(data, moderators)
      moderators <- ncol(data)
    }
  } else {
    if("glinternet" %in% lambda){moderators <- 1:ncol(data)}
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
  if(!all(lambda %in% c("CV", "glinternet")) & is.null(adaGam)){output$call[c("seed", "which.lam", "folds", "measure")] <- NULL}
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
  output
}

##### plotNet: Plot model results
plotNet <- function(object, which.net = "temporal", predict = NULL, names = TRUE, 
                    scale = FALSE, lag = NULL, nodewise = FALSE, con = "adjR2", 
                    cat = "nCC", layout = "spring", elabs = FALSE, elsize = 1, 
                    mnet = FALSE, plot = TRUE, ...){
  if(class(object) != "list"){
    stopifnot(ncol(object) == nrow(object))
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
    predict <- NULL
    nodewise <- FALSE
    names <- colnames(object)
    object <- list(adjMat = object, edgeColors = getEdgeColors(object))
  }
  if("SURnet" %in% names(object)){object <- object$SURnet}
  if(!is.null(names)){if(any(names == TRUE)){
    names <- tryCatch({colnames(object$data)}, error = function(e){c(1:length(object$mods))})
  } else if(any(names == FALSE)){
    names <- c(1:length(object$mods))}
  }
  if("temporal" %in% names(object)){
    if(!is.null(which.net)){
      if(class(which.net) == "numeric"){which.net <- c("t", "c", "b")[which.net]}
      net <- match.arg(which.net, c("temporal", "contemporaneous", "between"))
      object$adjMat <- object[[net]]$adjMat
      object$edgeColors <- object[[net]]$edgeColors
    } else {
      stop("Use the argument 'which.net' to indicate which network to plot")
    }
  }
  if(is.null(lag) & "adjMats" %in% names(object)){
    stop("More than one lag modeled; need to specify which to plot")
  } else if(!"adjMats" %in% names(object)){
    if(any(grepl("lag", colnames(object$adjMat)))){lag <- 1}
  }
  if(!is.null(predict) & !mnet){
    con <- match.arg(con, choices = c("R2", "adjR2", "MSE", "RMSE"))
    cat <- match.arg(cat, choices = c("nCC", "CC", "CCmarg"))
    type <- object$call$type
    if("ggm" %in% names(attributes(object))){
      if(is.list(type)){
        type <- unname(sapply(object$fitobj, attr, "family"))
        if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
        if("binomial" %in% type){type[type == "binomial"] <- "c"}
      } else {
        type <- ifelse(length(type) == 1, ifelse(
          type %in% c("g", "gaussian"), list(rep("g", ncol(data))), 
          list(rep("c", ncol(data)))), list(type))[[1]]
      }
    }
    if(class(predict) == "list"){
      predict <- list(predictNet(predict, all = FALSE, scale = scale), 
                      predictNet(object, all = FALSE, scale = scale))
      stopifnot(length(predict) == 2)
      pie <- list(); pieColor <- list()
      for(i in 1:length(type)){
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
      for(i in 1:length(type)){
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
    } else {
      qgraph(input = t(object$adjMat), layout = layout, edge.color = t(object$edgeColors),
             labels = names, DoNotPlot = !plot, pie = pie, pieColor = pieColor, 
             edge.labels = elabs, edge.label.cex = elsize, ...)
    }
  } else {
    if(nodewise == FALSE){
      if(mnet == FALSE){
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
plotNet2 <- function(object, whichNets = NULL, titles = c("PDC ", "PCC "), ...){
  if("SURnet" %in% names(object)){object <- object$SURnet}
  if(!is.null(whichNets)){
    tt <- ifelse(all(tolower(whichNets) == "all"), 3, 2)
    if(tt == 3){whichNets <- c("temporal", "contemporaneous", "between")}
    if(tt == 2){
      l <- averageLayout(plotNet(object, which.net = whichNets[1], plot = FALSE),
                         plotNet(object, which.net = whichNets[2], plot = FALSE))
    } else if(tt == 3){
      l <- averageLayout(plotNet(object, which.net = whichNets[1], plot = FALSE),
                         plotNet(object, which.net = whichNets[2], plot = FALSE),
                         plotNet(object, which.net = whichNets[3], plot = FALSE))
      if(length(titles) == 2){titles <- c("Temporal Effects", "Contemporaneous Effects", "Between-subject Effects")}
    }
  } else {
    tt <- 2
    l <- averageLayout(plotNet(object, which.net = "temporal", plot = FALSE), 
                       plotNet(object, which.net = "contemporaneous", plot = FALSE))
  }
  layout(t(1:tt))
  if(tt == 2){
    if(all(titles == c("PDC ", "PCC "))){
      if("lmer" %in% names(attributes(object)) & all(whichNets == c("temporal", "contemporaneous"))){
        title1 <- "Temporal Effects"
        title2 <- "Contemporaneous Effects"
        titles <- c(title1, title2)
      } else {
        title1 <- "Partial Directed Correlations"
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
    plotNet(object, which.net = "temporal", layout = l, title = title1, ...)
    plotNet(object, which.net = "contemporaneous", layout = l, title = title2, ...)
  }
}

##### predictNet: Calculate prediction error
predictNet <- function(object, data = NULL, all = FALSE, scale = FALSE){
  if("SURnet" %in% names(object)){object <- object$SURnet}
  #if("ggm" %in% names(attributes(object))){object$call <- NULL}
  if(is.null(data)){
    x <- try(data <- object$data)
    if(class(x) == "try-error"){stop("Must supply dataset")}
  }
  mods <- object$mods
  type <- object$call$type
  moderators <- object$call$moderators
  if("lags" %in% names(object$call)){lags <- object$call$lags} else {lags <- NULL}
  if(!is.null(moderators)){if(length(moderators) == 1){if(moderators == 0){moderators <- NULL}}}
  if(is.null(colnames(data))){colnames(data) <- paste0("V", 1:ncol(data))}
  predObj <- list(); Y <- list(); probObj <- vector("list", length = ncol(data))
  if("ggm" %in% names(attributes(object))){
    mods <- lapply(mods, function(z) list(model = z$model))
    if(is.list(type)){
      type <- unname(sapply(object$fitobj, attr, "family"))
      if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
      if("binomial" %in% type){type[type == "binomial"] <- "c"}
    } else {
      type <- ifelse(length(type) == 1, ifelse(
        type %in% c("g", "gaussian"), list(rep("g", ncol(data))), 
        list(rep("c", ncol(data)))), list(type))[[1]]
    }
    for(i in 1:length(mods)){
      Y[[i]] <- data[, i]
      predObj[[i]] <- predict(object = object$fitobj[[i]])
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
  } else {
    output <- errors
  }
  output
}

##### predictDelta: Calculate change in prediction for across nested models
predictDelta <- function(mod2, mod1, scale = FALSE){
  #if("ggm" %in% names(attributes(mod1))){mod1$call <- NULL}
  #if("ggm" %in% names(attributes(mod2))){mod2$call <- NULL}
  type <- mod2$call$type
  if("ggm" %in% names(attributes(mod2))){
    if(is.list(type)){
      type <- unname(sapply(object$fitobj, attr, "family"))
      if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
      if("binomial" %in% type){type[type == "binomial"] <- "c"}
    } else {
      type <- ifelse(length(type) == 1, ifelse(
        type %in% c("g", "gaussian"), rep("g", ncol(data)), rep("c", ncol(data))), type)
    }
  } else {stopifnot(mod2$call$type == mod1$call$type)}
  #if("ggm" %in% names(attributes(mod2))){type <- rep("g", ncol(mod2$data))}
  errors <- list(predictNet(mod1, all = FALSE, scale = scale), predictNet(mod2, all = FALSE, scale = scale))
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

################################################################################
################################################################################
##### SURsampler: Sample data for fitting system of lagged SURs
SURsampler <- function(B = NULL, S, n, seed = NULL, beta, beta2 = NULL,
                       full = TRUE, mu = 0, cholesky = FALSE, time = FALSE){
  p <- ncol(S)
  if(!is.null(B)){
    if(class(B) == "list"){B <- do.call(rbind, B)}
    if(missing(beta)){
      if(ncol(B) == (p + 1)){
        beta <- B
      } else {
        beta <- B[,1:(p + 1)]
        beta2 <- B[,-c(1:(p + 1)), drop = FALSE]
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
  if(cholesky == FALSE){
    R <- MASS::mvrnorm(n = n, mu = mu, Sigma = S)
  } else {
    R <- matrix(rnorm(p * n), ncol = p) %*% chol(S)
  }
  X <- matrix(NA, ncol = p, nrow = n + 1)
  if(!is.null(seed)){set.seed(seed)}
  X[1,] <- rnorm(p)
  it <- 0
  if(is.null(beta2)){
    repeat{
      it <- it + 1
      for(i in 1:p){
        X[it + 1, i] <- beta[i, 1] + sum(X[it,] * beta[i, -1]) + R[it, i]
      }
      if(it == n){break}
    }
    Y <- X[-1,]
  } else {
    Y <- X
    X2 <- X[,1] * X[,-1]
    X <- cbind(X, X2)
    beta <- cbind(beta, beta2)
    repeat{
      it <- it + 1
      for(i in 1:p){
        Y[it + 1, i] <- beta[i, 1] + sum(X[it,] * beta[i, -1]) + R[it, i]
        X <- cbind(Y, (Y[,1] * Y[,-1]))
      }
      if(it == n){break}
    }
    Y <- Y[-1,]
  }
  tx <- Sys.time() - tx
  if(time == TRUE){cat("Completed in", round(tx, 3), attr(tx, "units"), "\n")}
  colnames(Y) <- paste0("X", 1:ncol(Y), ".y")
  X <- X[-nrow(X),]
  colnames(X) <- paste0("X", 1:ncol(X))
  if(!is.null(beta2)){
    ints <- (ncol(beta) - ncol(beta2)):ncol(X)
    for(i in ints){colnames(X)[i] <- paste0("X1:X", which(ints == i) + 1)}
  }
  dat <- list(Y = Y, X = X)
  if(full == TRUE){
    dat <- list(Y = Y, X = X, full = data.frame(do.call(cbind, dat)))
    if(!is.null(beta2)){colnames(dat$full)[ncol(dat$Y) + ints] <- colnames(dat$X)[ints]}
  }
  if(time == TRUE){attributes(dat)$time <- tx}
  dat
}

##### SURfit: fit SUR model with or without constraints
SURfit <- function(dat, varMods = NULL, mod = "min", maxiter = 1000, m = NULL, ...){
  mod <- match.arg(mod, choices = c("min", "1se"))
  if(is.null(m)){m <- ifelse(any(grepl(":", names(dat$full))), 1, 0)}
  eqs <- makeEqs(dat = dat, varMods = varMods, mod = mod, m = m)
  fit <- systemfit::systemfit(formula = eqs, method = "SUR", data = dat$full, 
                              maxiter = maxiter, ...)
  fit
}

##### SURll: log-likelihood of SUR model
SURll <- function(fit, s = "res"){
  s <- match.arg(tolower(s), choices = c("res", "dfres", "sigma"))
  resid <- residuals(fit)
  residCov <- getCoefs(fit = fit, mat = s)
  residCovInv <- solve(residCov)
  resid <- as.matrix(resid)
  nEq <- ncol(resid)
  if((nrow(resid) * ncol(resid)) <= 40000){
    ll <- sum((-nEq/2) * log(2 * pi) - .5 * log(det(residCov)) - .5 * diag(resid %*% residCovInv %*% t(resid)))
  } else {
    ll <- 0
    for(i in 1:nrow(resid)){
      ll <- ll - (nEq/2) * log(2 * pi) - .5 * log(det(residCov)) - .5 * resid[i, , drop = FALSE] %*% residCovInv %*% t(resid[i, , drop = FALSE])
    }
  }
  df <- fit$rank + (nEq * (nEq + 1))/2
  c(logLik = as.numeric(ll), df = df)
}

##### SURtable: print fit indices for 3 SUR models
SURtable <- function(fits, ind = "RMSE", verbose = TRUE, d = 3,
                     labs = c("fitOLS", "fit0", "fit1se")){
  stopifnot(class(fits) == "list")
  if(ind %in% c("RSS", "rss")){ind <- "SSR"}
  ind <- match.arg(ind, choices = c("RMSE", "MSE", "SSR", "R2", "adjR2", "omnibus", "lrtest"))
  measure <- switch(ind, "RMSE" = "sigma", "MSE" = "sigma", "R2" = "r.squared", "lrtest" = "lrtest",
                    "adjR2" = "adj.r.squared", "SSR" = "ssr", "omnibus" = "omnibus")
  if(is.null(names(fits))){
    if(all(labs == c("fitOLS", "fit0", "fit1se"))){
      names(fits) <- labs[rev(order(sapply(fits, function(z) z$rank)))]
    } else {
      names(fits) <- labs
    }
  }
  p <- length(fits[[1]]$eq)
  if(measure == "omnibus"){
    dat <- matrix(NA, ncol = 8, nrow = 3)
    dat[,1] <- sapply(fits, logLik)
    dat[,2] <- sapply(fits, function(z) summary(z)$df[1] + (p * (p + 1))/2)
    dat[,3] <- sapply(fits, AIC)
    dat[,4] <- sapply(fits, BIC)
    dat[,5] <- sapply(fits, function(z) sum(colSums(residuals(z)^2)))
    dat[,6] <- sapply(fits, function(z) summary(z)$detResidCov)
    dat[,7] <- sapply(fits, function(z) summary(z)$ols.r.squared)
    dat[,8] <- sapply(fits, function(z) summary(z)$mcelroy.r.squared)
    colnames(dat) <- c("LL", "df", "AIC", "BIC", "SSR", "detSigma", "OLS.R2", "McElroy.R2")
    rownames(dat) <- names(fits)
  } else if(measure == "lrtest"){
    lrdat <- t(sapply(fits, SURll))
    lrdat <- lrdat[order(lrdat[,"df"], decreasing = TRUE),]
    x <- t(combn(1:3, 2))
    dat <- data.frame(LL1 = rownames(lrdat)[x[,1]], LL2 = rownames(lrdat)[x[,2]], 
                      Chisq = NA, Df = NA, p.value = NA, sig = NA)
    for(i in 1:3){
      dat[i,3] <- (2 * (lrdat[x[i,1], 1] - lrdat[x[i,2], 1]))
      dat[i,4] <- abs(lrdat[x[i,1], 2] - lrdat[x[i,2], 2])
      dat[i,5] <- pchisq(dat[i,3], dat[i,4], lower.tail = FALSE)
      dat[i,6] <- ifelse(dat[i,5] < .001, "***", 
                         ifelse(dat[i,5] < .01, "**", 
                                ifelse(dat[i,5] < .05, "*", 
                                       ifelse(dat[i,5] < .1, ".", ""))))
    }
    dat[,c(3,5)] <- round(dat[,c(3,5)], digits = d)
    RMSEA <- function(X2, df, N){
      rmsea <- sqrt(max(c(((X2/N)/df) - (1/N), 0)))
      lower.l <- function(lambda){(pchisq(X2, df = df, ncp = lambda) - .95)}
      lambda.l <- tryCatch({uniroot(f = lower.l, lower = 0, upper = X2)$root}, error = function(e){0})
      rmsea.lower <- sqrt(lambda.l/(N * df))
      upper.l <- function(lambda){(pchisq(X2, df = df, ncp = lambda) - .05)}
      lambda.u <- tryCatch({uniroot(f = upper.l, lower = 0, upper = max(N, X2 * 4))$root}, error = function(e){1})
      rmsea.upper <- sqrt(lambda.u/(N * df))
      rmsea.pvalue <- 1 - pchisq(X2, df = df, ncp = (N * df * (.05^2)))
      out <- c(lower = rmsea.lower, RMSEA = rmsea, upper = rmsea.upper, p.value = rmsea.pvalue)
      out
    }
    N <- fits[[1]]$df.residual + fits[[1]]$rank
    rmsea <- round(t(mapply(RMSEA, dat$Chisq, dat$Df, N)), 3)
    rownames(rmsea) <- rownames(dat)
    lrdat <- cbind(lrdat, AIC = (-2 * lrdat[,1]) + (2 * lrdat[,2]), 
                   BIC = (-2 * lrdat[,1]) + (lrdat[,2] * log(N)))
    dat <- list(LRtest = dat, RMSEA = rmsea, LLs = lrdat)
    if(verbose == TRUE){
      cat("\nIf difference is NOT significant, 'LL2' is preferred\nIf difference IS significant, 'LL1' is preferred\n\n")}
  } else {
    fs <- lapply(fits, function(z) lapply(z$eq, summary))
    dat <- sapply(fs, function(z) sapply(z, function(y) y[[measure]]))
    if(ind == "MSE"){dat <- dat^2}
    rownames(dat) <- sapply(1:p, function(z) as.character(fits[[1]]$eq[[z]]$terms[[2]]))
    colnames(dat) <- paste0(colnames(dat), ".", ind)
  }
  if(measure != "lrtest"){dat <- round(dat, digits = d)}
  dat
}

##### SURnet: create temporal and contemporaneous network of SUR results
SURnet <- function(fit, dat, s = "sigma", m = NULL, type = "g"){
  if(class(dat) == "list"){data <- dat$X} else {data <- dat}
  if(any(grepl(":", colnames(data)))){data <- data[,!grepl(":", colnames(data))]}
  if(length(type) == 1 & ncol(data) > 1){type <- rep(type, ncol(data))}
  if(is.null(m)){m <- ifelse(any(grepl(":", names(dat$full))), 1, 0)}
  s <- match.arg(tolower(s), choices = c("sigma", "res", "dfres"))
  call <- list(type = type, moderators = m, lags = 1, residMatType = s)
  beta <- sparsify(getCoefs(fit = fit, mat = "beta"))
  kappa <- sparsify(solve(getCoefs(fit = fit, mat = s)))
  mainEffects <- beta[,2:(ncol(kappa) + 1)]
  interactions <- if(m == 1){beta[,(ncol(kappa) + 2):ncol(beta)]} else {NA}
  PCC <- sparsify(getCoefs(fit = fit, mat = "pcor"))
  PDC <- sparsify(mainEffects/(sqrt(diag(solve(kappa)) %o% diag(kappa) + mainEffects^2)))
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    if(!any(grepl("lag", colnames(adjMat)))){
      colnames(colMat) <- paste0(colnames(adjMat), ".lag1.")
    } else {
      colnames(colMat) <- colnames(adjMat)
    }
    rownames(colMat) <- rownames(adjMat)
    colMat
  }
  temporal <- list(adjMat = PDC, edgeColors = getEdgeColors(PDC), beta = beta)
  contemporaneous <- list(adjMat = PCC, edgeColors = getEdgeColors(PCC), kappa = kappa)
  colnames(temporal$adjMat) <- colnames(temporal$edgeColors)
  colnames(contemporaneous$adjMat) <- colnames(contemporaneous$edgeColors)
  mods <- lapply(1:nrow(beta), function(z) list(model = as.matrix(beta[z,], ncol = 1)))
  surNet <- list(call = call, temporal = temporal, contemporaneous = contemporaneous,
                 interactions = interactions, mods = mods, data = data)
  if(m == 0){surNet$interactions <- NULL}
  surNet
}

################################################################################
################################################################################
##### results: print or save SUR model output (true+estimated matrices for beta/sigma/cor)
results <- function(true, est, x = "beta", d = NULL, labels = NULL,
                    intercepts = TRUE, keep = FALSE){
  x <- match.arg(x, choices = c("beta", "sigma", "cor", "res", "dfres", "pcor"))
  if(dim(true)[1] == dim(true)[2] & all(diag(true) == 1) & x == "beta"){x <- "sigma"}
  if(dim(true)[1] != dim(true)[2] & x != "beta"){x <- "beta"}
  x2 <- ifelse(x == "beta", "BETA", ifelse(x %in% c("sigma", "res", "dfres"), "RESID.COV", 
                                           ifelse(x == "cor", "RESID.COR", "RESID.PCOR")))
  if(x == "pcor"){
    true <- tryCatch({-solve(true)}, error = function(e){
      -corpcor::pseudoinverse(true)})
    diag(true) <- -diag(true)
    delta <- (1/sqrt(diag(true)))
    true <- t(delta * true) * delta
    if(is.null(d)){d <- 2}
  }
  if(is.null(d)){d <- ifelse(x == "beta", 2, 1)}
  int <- ifelse(x == "beta", intercepts, TRUE)
  if(is.null(labels)){
    ilabs <- rownames(getCoefs(fit = est, mat = x))
    jlabs <- colnames(getCoefs(fit = est, mat = x))
  } else if(class(labels) == "systemfit"){
    ilabs <- rownames(getCoefs(labels, mat = x))
    jlabs <- colnames(getCoefs(labels, mat = x))
  } else {ilabs <- labels$ilabs; jlabs <- labels$jlabs}
  if(keep == FALSE){
    if(int == TRUE){
      cat(paste0("\nTRUE ", x2, "\n"))
      print(sparsify(mat = true, d = d, ilabs = ilabs, jlabs = jlabs))
      cat(paste0("\nESTIMATED ", x2, "\n"))
      print(sparsify(mat = getCoefs(fit = est, mat = x), d = d, 
                     ilabs = ilabs, jlabs = jlabs))
    } else {
      cat(paste0("\nTRUE ", x2, "\n"))
      print(sparsify(mat = true, d = d, ilabs = ilabs, jlabs = jlabs)[,-1])
      cat(paste0("\nESTIMATED ", x2, "\n"))
      print(sparsify(mat = getCoefs(fit = est, mat = x), d = d,
                     ilabs = ilabs, jlabs = jlabs)[,-1])
    }
  } else {
    m1 <- sparsify(mat = true, d = d, ilabs = ilabs, jlabs = jlabs)
    m2 <- sparsify(mat = getCoefs(fit = est, mat = x), d = d, 
                   ilabs = ilabs, jlabs = jlabs)
    if(int == FALSE){m1 <- m1[,-1]; m2 <- m2[,-1]}
    m <- list(m1, m2)
    names(m) <- c(paste0("true", Hmisc::capitalize(x)), paste0("est", Hmisc::capitalize(x)))
    m
  }
}

##### getCoefs: extract beta matrix and/or sigma matrix from systemfit model
getCoefs <- function(fit, mat = "beta", d = NULL, sparse = FALSE, ints = TRUE, data = NULL){
  if("SURfit" %in% names(fit)){fit <- fit$SURfit}
  if(class(fit) == "systemfit"){
    mat <- match.arg(arg = mat, choices = c("beta", "sigma", "res", "dfres", "cor", "pcor"))
    if(mat == "beta"){
      ynames <- c(); n <- list()
      for(i in 1:length(fit$eq)){
        ynames[i] <- as.character(fit$eq[[i]]$terms[[2]])
        n[[i]] <- names(coef(fit$eq[[i]]))
      }
      N <- ifelse(any(grep(":", names(coef(fit)))), length(fit$eq) * 2, length(fit$eq) + 1)
      b <- t(matrix(0, ncol = length(fit$eq), nrow = N))
      rownames(b) <- ynames
      if(max(sapply(n, length)) != ncol(b)){
        data <- tryCatch({eval(fit$call$data)}, error = function(e){
          if(!is.null(data)){
            data
          } else {stop("Need to supply original data to create beta matrix")}
        })
        if(class(data) == "list"){data <- data$full}
        colnames(b) <- c("(Intercept)", colnames(data)[!colnames(data) %in% ynames])
      } else {
        colnames(b) <- n[[which.max(sapply(n, length))]]
      }
      for(i in 1:nrow(b)){
        for(j in 1:length(n[[i]])){
          b[i, names(coef(fit$eq[[i]]))[j]] <- coef(fit$eq[[i]])[j]
        }
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
  if(!is.null(d)){b <- round(b, digits = d)}
  if(sparse == TRUE){b <- sparsify(b)}
  if(ints == FALSE & mat == "beta"){b <- b[, -1, drop = FALSE]}
  b
}

##### sparsify: convert any (i x j) matrix into a sparse matrix
sparsify <- function(mat, d = NULL, ilabs = NULL, jlabs = NULL){
  stopifnot(length(dim(mat)) == 2)
  p <- ncol(mat)
  n <- nrow(mat)
  if(!is.null(rownames(mat)) & is.null(ilabs)){ilabs <- rownames(mat)}
  if(!is.null(colnames(mat)) & is.null(jlabs)){jlabs <- colnames(mat)}
  if(class(mat) == "data.frame"){mat <- as.matrix(mat, ncol = p, nrow = n)}
  if(!is.null(d)){mat <- round(mat, digits = d)}
  dat <- data.frame(i = rep(1:n, p), j = rep(1:p, each = n), x = as.vector(mat))
  dat <- dat[dat$x != 0 & !is.na(dat$x), ]
  dat2 <- Matrix::sparseMatrix(i = dat$i, j = dat$j, x = dat$x, dims = c(n, p))
  if(!is.null(ilabs)){rownames(dat2) <- ilabs}
  if(!is.null(jlabs)){colnames(dat2) <- jlabs}
  dat2
}

##### makeEqs: Generates a set of linear equations based on lagged matrices
makeEqs <- function(dat = NULL, varMods = NULL, mod = "min", type = "g", m = 1, 
                    scale = FALSE, addLast = TRUE, names = TRUE, datBack = FALSE){
  if(!is.null(varMods)){
    if(attr(varMods, "criterion") != "CV"){mod <- "min"}
    mod <- match.arg(mod, choices = c("min", "1se"))
    mod <- ifelse(mod == "min" | attr(varMods, "method") == "regsubsets", 1, 2)
    if(class(dat) == "list"){
      data <- dat$X
      if(any(grepl(":", colnames(data)))){data <- data[,!grepl(":", colnames(data))]}
    }
    newDat <- datMat(data = data, type = type, lags = 1, m = m, scale = scale)
    if(names == TRUE){
      colnames(newDat$X) <- colnames(dat$X)
      colnames(newDat$full)[!colnames(newDat$full) %in% colnames(dat$Y)] <- colnames(dat$X)
    }
    Ys <- colnames(newDat$Y)
    Xs <- colnames(newDat$X)
    eqs <- list()
    for(i in 1:length(Ys)){
      mains <- varMods[[i]][[mod]][!grepl(":", varMods[[i]][[mod]])]
      ints <- varMods[[i]][[mod]][grepl(":", varMods[[i]][[mod]])]
      if(m == 1){ints <- ints[grepl("X1", ints)]}
      predictors <- c(mains, ints)
      if(length(predictors) == 0){predictors <- 1}
      eqs[[i]] <- as.formula(paste0(Ys[i], " ~ ", paste(predictors, collapse = " + ")))
    }
    if(addLast == TRUE){
      for(j in 1:length(dat)){
        newDat[[j]] <- rbind(newDat[[j]], rep(NA, ncol(newDat[[j]])))
        newDat[[j]][nrow(newDat[[j]]),] <- dat[[j]][nrow(dat[[j]]),]
        rownames(newDat[[j]]) <- NULL
      }
    }
    newDat2 <- list(Y = newDat$Y, X = newDat$X, full = newDat$full, eqs = eqs)
    if(datBack == TRUE){return(newDat2)} else {return(newDat2$eqs)}
  } else {
    if(!is.null(dat)){
      x <- dat$X
      y <- dat$Y
    }
    eqs <- list()
    for(i in 1:ncol(y)){
      eqs[[i]] <- as.formula(paste(colnames(y)[i], "~", paste(colnames(x), collapse = " + ")))
    }
    return(eqs)
  }
}

##### matrixDist: compute similarity between two matrices
matrixDist <- function(mat1, mat2, ind = "cosine", similarity = TRUE, distMat = FALSE){
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
  if(ind == "cosine"){
    f1 <- norm(x = mat1, type = "F")
    f2 <- norm(x = mat2, type = "F")
    top <- sum(diag(t(mat1) %*% mat2))
    return(ifelse(similarity == TRUE, top/(f1 * f2), 1 - (top/(f1 * f2))))
  }
  if(ind == "correlation"){return(cor(as.vector(mat1), as.vector(mat2)))}
  if(ind == "ssd"){return(sum(d^2))}
  if(ind == "mse"){return(sum(d^2)/length(d))}
  if(ind == "rmse"){return(sqrt(sum(d^2)/length(d)))}
  if(ind == "mae"){return((sum(abs(d))/length(d)))}
  if(ind == "msd"){return(sum(d/length(d)))}
}

##### compareTrue: compute estimated and true parameter matrices
compareMats <- function(B, S, fits, d = 5, ind = "cosine", similarity = TRUE,
                        labs = c("fitOLS", "fit0", "fit1se")){
  x <- matrix(NA, ncol = 5, nrow = 3); s <- similarity
  x[,1] <- sapply(fits, function(z){matrixDist(B, getCoefs(z, "beta"), ind, s)})
  x[,2] <- sapply(fits, function(z){matrixDist(S, getCoefs(z, "sigma"), ind, s)})
  x[,3] <- sapply(fits, function(z){matrixDist(S, getCoefs(z, "res"), ind, s)})
  x[,4] <- sapply(fits, function(z){matrixDist(S, getCoefs(z, "cor"), ind, s)})
  k <- tryCatch({-solve(S)}, error = function(e){-corpcor::pseudoinverse(S)})
  diag(k) <- -diag(k)
  K <- t((1/sqrt(diag(k))) * k) * (1/sqrt(diag(k)))
  x[,5] <- sapply(fits, function(z){matrixDist(K, getCoefs(z, "pcor"), ind, s)})
  if(is.null(names(fits))){
    if(all(labs == c("fitOLS", "fit0", "fit1se"))){
      rownames(x) <- labs[rev(order(sapply(fits, function(z) z$rank)))]
    } else {rownames(x) <- labs}
  }
  colnames(x) <- c("Beta", "Sigma", "ResidCov", "ResidCor", "ResidPcor")
  round(x, digits = d)
}

################################################################################
################################################################################
##### datMat: Create a list of matrices for fitting lagged models
datMat <- function(data, type, m = NULL, lags = 1, scale = FALSE, 
                   center = TRUE, name = TRUE, full = TRUE){
  if(length(type) == 1 & ncol(data) > 1){type <- rep(type, ncol(data))}
  if(!is.null(m)){if(m != 1){m <- NULL}}
  if(!all(c(center, scale) == TRUE)){
    data <- scale(data, center, scale)
    scale <- FALSE
  }
  dat <- list(); Y <- matrix(NA, ncol = ncol(data), nrow = nrow(data) - lags)
  for(i in 1:ncol(data)){
    dat[[i]] <- setup(data = data, type = type, y = i, lags = lags, scale = scale)
    Y[,i] <- dat[[i]][,1]
  }
  if(!is.null(colnames(data))){colnames(Y) <- paste0(colnames(data), ".y")}
  rownames(Y) <- 1:nrow(Y)
  X <- interactionMatrix(data = dat[[1]], y = 1, type = type, moderators = m, lags = lags)
  if(name == TRUE & is.null(m)){colnames(X) <- colnames(data)}
  mats <- list(Y = Y, X = X)
  if(full == TRUE){
    mats <- list(Y = Y, X = X, full = data.frame(do.call(cbind, mats)))
    p <- ncol(mats$Y)
    if(!is.null(m)){
      colnames(mats$full)[((p * 2) + 1):ncol(mats$full)] <- colnames(X)[(p + 1):ncol(X)]
    }
  }
  mats
}

##### fitHierLASSO: performs variable selection using the hierarchical LASSO
fitHierLASSO <- function(data, yvar, type = "g", m = NULL, criterion = "CV", 
                         method = "glinternet", gamma = .5, nfolds = 10, 
                         nlam = 50, scale = TRUE, lags = NULL, useSE = TRUE, 
                         diag = FALSE, outMsgs = FALSE){
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
          if(is.null(lags)){m2 <- m} else {m2 <- mains[-m]}
          x2 <- c(colnames(x), paste0(colnames(x)[-m2], ":", colnames(x)[m2]))
          betas <- betas[which(names(betas) %in% x2)]
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
        if(m == 0){
          p <- length(mains)
        } else {
          p <- length(c(colnames(x), paste0(colnames(x)[-m], ":", colnames(x)[m])))
        }
      }
      ic_lambda <- -2 * LL_models + n_neighbors * ifelse(
        criterion == "AIC", 2, log(nrow(x))) + ifelse(
          criterion == "EBIC", 2 * gamma * n_neighbors * log(p), 0)
      betas <- allCoefs[[which.min(ic_lambda)]]
      lambda_min <- fit$lambda[which.min(ic_lambda) + 1]
      coefs <- Matrix(betas, sparse = TRUE)
      rownames(coefs) <- names(betas)
      fitobj <- list(fit = fit, fit0 = NA, crit = ic_lambda)
      if(ifelse(!is.null(m), ifelse(m == 0, FALSE, TRUE), TRUE) & method == "hiernet"){
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
varSelect <- function(dat, m = NULL, criterion = "EBIC", method = "glinternet",
                      seed = NULL, nlam = NULL, nfolds = 10, lags = NULL, 
                      type = "g", useSE = TRUE, verbose = TRUE, 
                      scale = FALSE, all = FALSE, gamma = .5){
  if(!is.null(m)){if(all(m == 0)){all <- TRUE}}
  if(class(dat) != "list"){
    if(is.null(m) | all == TRUE | !is.null(lags)){
      if(!is.null(lags)){dat <- datMat(data = dat, type = type, scale = scale)}
      if(is.null(lags)){dat <- list(Y = dat, X = dat)}
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
      stopifnot(all == FALSE)
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
  }
  if(!is.null(seed)){set.seed(seed)}
  p <- ncol(dat$Y); data <- dat$X; Y <- dat$Y
  method <- match.arg(
    tolower(method), 
    c("hiernet", "glinternet", "subset", "backward", "forward", "seqrep"))
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
      if(all == TRUE & !is.null(m)){if(all(m != 0)){
        m <- ifelse(i %in% m, list(NULL), ifelse(m < i, list(m), list(m - 1)))[[1]]}
      }
      tx <- Sys.time()
      hierMods[[i]] <- fitHierLASSO(data = data, yvar = i, type = type, m = m,
                                    nlam = nlam, nfolds = nfolds, lags = lags,
                                    method = method, useSE = useSE, scale = scale, 
                                    criterion = criterion, gamma = gamma)
      t[i] <- tx <- Sys.time() - tx
      names(t)[i] <- attr(tx, "units")
      if(all == TRUE & exists("m0", inherits = FALSE)){m <- m0}
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
    return(hierMods)
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

################################################################################
################################################################################
##### resample: perform either bootstrapping or multi-sample splits for variable selection
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
  if(!is.null(m)){if(all(is.na(m))){m <- NULL}}
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
		if(i == 1){cat("\n")}
        t2 <- Sys.time()
        cat("************ Sample Split: ", i, "/", niter, " ************\n", sep = "")
      }
      set.seed(seeds[i])
      vars[[i]] <- varSelect(dat = train[[i]], m = m, criterion = criterion, 
                             method = varMethod, nfolds = nfolds, gamma = gamma, 
                             nlam = nlam, type = type, seed = varSeed, all = !exogenous,
                             verbose = ifelse(sampMethod == "stability", ifelse(verbose, "pbar", FALSE), verbose)) 
      if(sampMethod == "stability"){
        vars1[[i]] <- varSelect(dat = test[[i]], m = m, criterion = criterion,
                                method = varMethod, nfolds = nfolds, gamma = gamma,
                                nlam = nlam, type = type, seed = varSeed,
                                all = !exogenous, verbose = ifelse(verbose, "pbar", FALSE)
        if(verbose){
          t3 <- Sys.time() - t2
          cat(paste0("Time: ", round(t3, 2), " ", attr(t3, "units")), "\n\n")
        }
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
  varMods0 <- lapply(vars, function(z) lapply(z, function(zz) zz[[lam]]))
  if(sampMethod == "stability"){
    if(is.logical(threshold)){threshold <- .6}
    #dfmax in glmnet: q <- ceiling(sqrt(EV * p * (2 * threshold - 1)))
    varMods1 <- lapply(vars1, function(z) lapply(z, '[[', 1))
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
      coefs1 <- summary(fit0$fitobj[[z]])$coef[-1, , drop = FALSE]
      coefs2 <- confint(fit0$fitobj[[z]], level = 1 - alpha)[-1, , drop = FALSE]
      outFreqs <- data.frame(predictor = factor(allNames[[z]]), select = freqs >= threshold,
                             lower = coefs2[, 1], b = coefs1[, 1], upper = coefs2[, 2], 
                             Pvalue = coefs1[, 4], freq = freqs, split1 = s1freqs, split2 = s2freqs)
      rownames(outFreqs) <- 1:nrow(outFreqs)
      return(outFreqs)
    })))
    preout$nlam <- nlam
    attributes(varFreqs) <- list(
      names = vs, type = type, criterion = criterion, rule = rule, 
      gamma = gamma, center = center, exogenous = exogenous, 
      moderators = m, threshold = FALSE)
    split1 <- lapply(seq_len(niter), function(z){list(vars = varMods0[[z]], varMods = vars[[z]])})
    split2 <- lapply(seq_len(niter), function(z){list(vars = varMods1[[z]], varMods = vars1[[z]])})
    names(split1) <- names(split2) <- paste0("iter", 1:niter)
    out <- list(call = preout, freqs = varFreqs, split1 = split1, split2 = split2)
    out$stability <- stability(out)
    if(!saveMods){out <- out[-c(3:4)]}
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
      if(!all(is.na(adjCIs0[[i]]))){
        if(any(adjCIs0[[i]]$lower == -Inf)){
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
      names = vs, type = type, criterion = criterion, rule = rule, 
      gamma = gamma, center = center, exogenous = exogenous, 
      moderators = m, threshold = FALSE)
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
    iter_p1 <- list(vars = varMods0[[z]], fit = fits[[z]])
    iter_p2 <- rev(list(sample = list(data = samps[[z]], inds = sampInd[[z]]), varMods = vars[[z]]))
    if(!saveMods){iter_p1 <- iter_p1[-2]; iter_p2 <- iter_p2[-1]}
    return(append(iter_p1, iter_p2))
  })
  names(boots) <- paste0("iter", 1:niter)
  out <- list(call = preout)
  if(saveMods){
    finalFits <- finalCoefs <- vector("list", 2)
    for(j in 1:2){
      varMods <- selectMod(obj = adjCIs0, select = c("select", "select_ci")[j])
      finalFits[[j]] <- fitNetwork(data = data, moderators = m, type = varMods, rule = rule, 
                                   threshold = threshold, gamma = gamma, center = center, 
                                   exogenous = exogenous, verbose = FALSE) #, ...)
      finalCoefs[[j]] <- getFitCIs(finalFits[[j]], allNames, alpha = alpha)
    }
    names(finalFits) <- c("fit0", "fit.5")
    names(finalCoefs) <- c("fitCIs", "fitCIs.5")
    out$fits <- append(finalFits, append(list(fit1 = fit0), finalCoefs))
    out$modLLs <- list(omnibus = modLL(out$fits[1:3]), 
                       modLL1 = modLL(out$fits[1:2], all = TRUE),
                       modLL2 = modLL(out$fits[2:3], all = TRUE))
  }
  out$adjCIs <- adjCIs0
  out$samples <- list(coefs = finalCoefs0, iters = boots)
  if(verbose){cat("\n"); print(Sys.time() - t1)}
  out
}


##### plotCoefs: plot coefficients from SURfit with confidence intervals
plotCoefs <- function(fit, true = FALSE, alpha = .05, plot = TRUE, col = "blue", 
                      flip = TRUE, data = NULL, select = FALSE, size = 1, 
                      labels = TRUE){
  if(class(fit) == "systemfit"){
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
      if(class(dat) == "list"){dat <- dat$full}
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
    if(!ifelse(select == TRUE, TRUE, FALSE)){
      dat$select <- ifelse(dat$lower <= 0 & dat$upper >= 0, FALSE, TRUE)
    }
    dat <- dat[dat$select, ]
  }
  if(plot == TRUE){
    plotCI <- function(dat, xlabs, true, flip = TRUE, size = 1, labels = TRUE){
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
      if(flip == TRUE){return(p + coord_flip())} else {return(p)}
    }
    plotCI(dat, xlabs = function(x){sub("[^_]*_", "", x)}, true = true, 
           flip = flip, size = size, labels = labels)
  } else {
    return(dat)
  }
}

##### selectMod: creates necessary input for fitNetwork when selecting variables
selectMod <- function(obj, data = NULL, fit = FALSE, select = "select", 
                      thresh = NULL, ascall = TRUE, type = "gaussian", ...){
  if(is.null(data)){ascall <- FALSE}
  cis <- c("adjCIs", "freqs")[which(c("adjCIs", "freqs") %in% names(obj))]
  if(length(cis) != 0){obj <- obj[[cis]]}
  if(any(grepl("adjCIs|freqs", names(obj)))){obj <- ifelse("adjCIs")}
  allNames <- lapply(lapply(obj, '[[', "predictor"), as.character)
  vs <- names(allNames) <- names(obj)
  if(is.numeric(select)){
    for(i in seq_along(obj)){obj[[i]]$select <- obj[[i]]$freq >= select}
    select <- "select"
  } else if(!is.null(thresh)){
    for(i in seq_along(obj)){obj[[i]]$select <- obj[[i]][[select]] >= thresh}
    select <- "select"
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
      length(type) == 1, ifelse(type == "gaussian", "g", "c"), 
      ifelse(type[i] %in% c("g", "gaussian"), "g", "c")
    )
  }
  if(!is.null(data)){
    args <- list(...)
    atts <- attributes(obj)
    atts[c("names", "criterion")] <- NULL
    atts$data <- data
    atts$type <- varMods
    if("err" %in% names(atts)){atts$err <- NULL}
    if(length(args) != 0){
      if(any(names(args) %in% names(atts))){
        dups <- names(atts)[which(names(atts) %in% names(args))]
        for(i in 1:length(dups)){atts[[dups[i]]] <- args[[dups[i]]]}
        args <- args[-which(names(args) %in% names(atts))]
      }
      if(length(args) != 0){atts <- append(atts, args)}
    }
    attr(varMods, "call") <- atts
  }
  if(fit & !is.null(data)){return(do.call(fitNetwork, atts))}
  if(ascall){return(attr(varMods, "call"))} else {return(varMods)}
}

##### getFitCIs: provides model coefficients with CIs 
getFitCIs <- function(fit, allNames, alpha = .05){
  if(all(c("mods", "call") %in% names(allNames))){
    allNames <- lapply(lapply(allNames$mods, '[[', "model"), function(z) rownames(z)[-1])
  }
  fits <- fit$fitobj
  fitCoefs <- suppressWarnings(suppressMessages(
    lapply(seq_along(fits), function(z){
      coefs1 <- summary(fits[[z]])$coef[-1, , drop = FALSE]
      if(nrow(coefs1) == 0){return(NA)}
      coefs2 <- confint(fits[[z]], level = 1 - alpha)[-1, , drop = FALSE]
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
  names(fitCoefs) <- names(fits)
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

##### fdr: False detection rates
fdr <- function(res, B, all = FALSE, ss = FALSE){
  if(ncol(B) == (2 * nrow(B))){B <- B[,-1]}
  B[B != 0] <- 1
  B <- t(B)
  m1 <- res$mod0$freqs * length(res$samps)
  m1 <- split(m1, rep(1:ncol(m1), each = nrow(m1)))
  for(i in 1:length(m1)){m1[[i]] <- cbind(m1[[i]], B[,i])}
  possiblePos <- as.vector(sapply(m1, function(z) sum(z[,2]))) * length(res$samps)
  possibleNeg <- as.vector(sapply(m1, function(z) sum(z[,2] == 0))) * length(res$samps)
  tp1 <- as.vector(sapply(m1, function(z) sum(z[,1] * z[,2])))
  fn1 <- possiblePos - tp1
  fp1 <- as.vector(sapply(m1, function(z) sum(z[z[,2] == 0, 1])))
  tn1 <- possibleNeg - fp1
  acc1 <- (tp1 + tn1)/(tp1 + fp1 + tn1 + fn1)
  prc1 <- tp1/(tp1 + fp1)
  fdr1 <- fp1/(tp1 + fp1)
  sens1 <- tp1/(tp1 + fn1)
  spec1 <- tn1/(tn1 + fp1)
  m2 <- res$mod1se$freqs * length(res$samps)
  m2 <- split(m2, rep(1:ncol(m2), each = nrow(m2)))
  for(i in 1:length(m2)){m2[[i]] <- cbind(m2[[i]], B[,i])}
  possiblePos <- as.vector(sapply(m2, function(z) sum(z[,2]))) * length(res$samps)
  possibleNeg <- as.vector(sapply(m2, function(z) sum(z[,2] == 0))) * length(res$samps)
  tp2 <- as.vector(sapply(m2, function(z) sum(z[,1] * z[,2])))
  fn2 <- possiblePos - tp2
  fp2 <- as.vector(sapply(m2, function(z) sum(z[z[,2] == 0, 1])))
  tn2 <- possibleNeg - fp2
  acc2 <- (tp2 + tn2)/(tp2 + fp2 + tn2 + fn2)
  prc2 <- tp2/(tp2 + fp2)
  fdr2 <- fp2/(tp2 + fp2)
  sens2 <- tp2/(tp2 + fn2)
  spec2 <- tn2/(tn2 + fp2)
  if(ss == TRUE){
    sss <- round(data.frame(Sens0 = sens1, Spec0 = spec1, Sens1se = sens2, Spec1se = spec2), 2)
    return(sss)
  }
  if(all == FALSE){
    round(data.frame(ACC0 = acc1, PRC0 = prc1, FDR0 = fdr1, Sens0 = sens1, 
                     Spec0 = spec1, ACC1se = acc2, PRC1se = prc2, FDR1se = fdr2, 
                     Sens1se = sens2, Spec1se = spec2), 2)
  } else {
    list(mod0 = data.frame(tp1, fp1, tn1, fn1, acc1, prc1, fdr1, sens1, spec1), 
         mod1se = data.frame(tp2, fp2, tn2, fn2, acc2, prc2, fdr2, sens2, spec2))
  }
}
makeFreqs <- function(res){
  niter <- length(res)
  allNames <- paste0("X", 1:4)
  allNames <- c(allNames, apply(combn(allNames, 2), 2, paste, collapse = ":"))[1:7]
  x <- list()
  for(i in 1:4){
    x[[i]] <- lapply(res, function(z) as.character(z[,i]))
    x[[i]] <- unname(unlist(x[[i]]))
    x[[i]] <- count(x[[i]][x[[i]] != "" & x[[i]] != "X2:X3" & x[[i]] != "X2:X4" & x[[i]] != "X3:X4"])
    x[[i]] <- x[[i]][match(allNames, x[[i]][,1]),]
  }
  m1 <- do.call(cbind, x)[,seq(2, 8, by = 2)]
  colnames(m1) <- paste0("X", 1:4, ".y")
  rownames(m1) <- allNames
  m1 <- as.matrix(m1)
  if(any(is.na(m1))){m1[is.na(m1)] <- 0}
  m1
}
getResids <- function(obj){
  r0 <- t(sapply(obj$fits, function(z) sapply(z$mod0$eq, function(zz) zz$df.residual)))
  r1 <- t(sapply(obj$fits, function(z) sapply(z$mod1se$eq, function(zz) zz$df.residual)))
  res0 <- lapply(obj$fits, function(z) residuals(z$mod0))
  res1 <- lapply(obj$fits, function(z) residuals(z$mod1se))
  colnames(r0) <- colnames(r1) <- colnames(res0[[1]])
  list(resid0 = res0, dfres0 = r0, resid1 = res1, dfres1 = r1)
}
mse <- function(x){
  rss0 <- lapply(x$residuals$resid0, function(z){
    apply(z, 2, function(zz) sum(zz^2))
  })
  rss0 <- do.call(rbind, rss0)
  mse0 <- matrix(0, ncol = ncol(x$residuals$dfres0), nrow = nrow(x$residuals$dfres0))
  for(i in 1:nrow(x$residuals$dfres0)){
    mse0[i,] <- rss0[i,]/x$residuals$dfres0[i,]
  }
  rss1 <- lapply(x$residuals$resid1, function(z){
    apply(z, 2, function(zz) sum(zz^2))
  })
  rss1 <- do.call(rbind, rss1)
  mse1 <- matrix(0, ncol = ncol(x$residuals$dfres1), nrow = nrow(x$residuals$dfres1))
  for(i in 1:nrow(x$residuals$dfres1)){
    mse1[i,] <- rss1[i,]/x$residuals$dfres1[i,]
  }
  list(mse0 = mse0, mse1 = mse1)
}

################################################################################
################################################################################

if(exists("clear")){
  if(clear == TRUE){
    message("Clearing .GlobalEnv")
    rm(clear)
    .modnets <- new.env()
    .modnets <- globalenv()
    suppressMessages(attach(.modnets))
    rm(list=ls())
  }
}
message("This is modnets 1.0.0")
