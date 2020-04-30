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
