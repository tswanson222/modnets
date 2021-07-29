#' Calculate prediction error from network models
#'
#' See the prediction error based on different statistics for either GGMs or SURs.
#' Also can compare and find the change values (such as R-squared change) between two networks
#'
#' @param object network
#' @param data dataset, or network for comparison
#' @param all logical
#' @param scale logical
#'
#' @return Table
#' @export
#'
#' @examples
#' 1 + 1
predictNet <- function(object, data = NULL, all = FALSE, scale = FALSE){
  # Use predictDelta if there are two networks
  if(is(data, 'ggm') | is(data, 'SURnet')){
    return(do.call(predictDelta, list(mod1 = object, mod0 = data, scale = scale)))
  }

  # SURnet vs. ggm
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

  # Calculate the type of error/predictive statistic
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

  # Compute errors
  errors <- lapply(1:ncol(data), function(z) calcError(Y[,z], yhat[,z], type[z], mods[[z]]$model))

  # Make everything into a pretty table
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

  # Decide what to return
  if(all == TRUE){
    output <- list(Y = Y, preds = yhat, probs = probObj, errors = errors)
    if(all(type == "g")){output$probs <- NULL}
  } else {
    output <- data.frame(Variable = errors[, 1], apply(errors[, -1], 2, unlist))
  }

  # RETURN
  return(output)
}

##### predictDelta: Calculate change in prediction for across nested models
predictDelta <- function(mod1, mod0, scale = FALSE){
  # SURnet vs. ggm
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
      } else if(length(type) == 1){
        type <- rep(substr(type, 1, 1), ncol(data))
      }
    } else {
      stopifnot(mod1$call$type == mod0$call$type)
    }
  }

  # Compute errors
  errors <- list(predictNet(mod0, all = FALSE, scale = scale), predictNet(mod1, all = FALSE, scale = scale))

  # See how many nodes are in each network, and whether they can be compared
  nodes <- unique(sapply(errors, nrow))
  if(length(nodes) > 1){stop('Cannot compare networks of different sizes')}

  # Collect variable types
  p <- ifelse(all(type == "c"), 4, ifelse(all(type == "g"), 5, 8))

  # Make matrix of comparisons
  delta <- data.frame(matrix(NA, nrow = length(nodes), ncol = p))
  for(i in seq_len(nodes)){
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

  # RETURN
  return(delta)
}

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
