##### nodewise: nodewise regression, with option to simulate/add moderators
nodewise <- function(data, mods = NULL, varMods = NULL, lambda = "min", center = TRUE, 
                     scale = FALSE, covariates = NULL, exogenous = TRUE, getEqs = FALSE){
  data <- dat <- dat0 <- data.frame(data)
  vs <- colnames(data)
  # Are there varMods? (type)
  if(!is.null(varMods)){
    # If varMods is list
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
      # If varMods is not list... 'type' now exists
    } else {
      type <- varMods
      if(is.null(mods)){
        varMods <- NULL
      } else if(length(mods) == 1){
        varMods <- NULL
        # If multiple moderators AND varMods is not a list..
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
    # If varMods is NULL...
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
  
  intMatrix2 <- function(data, mods = NULL, covariates = NULL, exogenous = TRUE){
    data <- data.frame(data)
    vars <- colnames(data)
    if(!is.null(mods)){
      mnames <- vars[mods]
      ynames <- switch(2 - exogenous, vars[-mods], vars)
      eqs <- lapply(ynames, function(z){
        intterms <- ifelse(length(mnames) == 1, mnames, 
                           paste0("(", paste(mnames, collapse = " + "), ")"))
        eq1 <- paste0(z, " ~ . ", ifelse(z %in% mnames, "^2", paste0("* ", intterms)))
        eq2 <- colnames(model.matrix(as.formula(eq1), data.frame(data)))[-1]
        return(paste0(z, " ~ ", paste(eq2, collapse = " + ")))
      })
    } else {
      eqs <- lapply(vars, function(z){
        eq1 <- as.formula(paste0(z, " ~ ."))
        eq2 <- colnames(model.matrix(eq1, data.frame(data)))[-1]
        return(paste0(z, " ~ ", paste(eq2, collapse = " + ")))
      })
    }
    if(!is.null(covariates)){
      eqs <- lapply(eqs, function(z) paste0(z, " + ", paste(names(covariates), collapse = " + ")))
    }
    return(lapply(eqs, as.formula))
  }
  
  # COVARIATES removing from data and making named list
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
  
  # MODS
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
  
  # CENTERING
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
  if(!is.null(varMods)){ # ONLY TRUE IF
    ints <- as.list(paste(names(varMods), "~", lapply(varMods, function(z) paste0(z[[lambda]], collapse = " + "))))
    if(!is.null(covariates) & ifelse("covariates" %in% names(attributes(varMods)), FALSE, TRUE)){
      ints <- lapply(ints, paste0, " + ", paste(names(covariates), collapse = " + "))
    }
  } else {
    ints <- intMatrix(data, mods, covariates)
    if(!is.null(mods) & !exogenous){ints <- append(ints, list(as.formula(paste0(names(mods), " ~ .^2"))))}
  }
  if(getEqs){return(lapply(ints, as.formula))}
  
  # CREATE MODELS
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
  
  # REJOINING DATA? AND TAKING IT AWAY
  if(!is.null(mods) & !exogenous){
    data0 <- data
    data <- dat0
  }
  
  # PROCESS MODELS
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
