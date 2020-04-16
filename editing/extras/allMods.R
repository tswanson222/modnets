allMods <- function(data, m, ebic = 0, cv = TRUE, fit = TRUE, zero = FALSE, ...){
  varMods1 <- varSelect(data, m, "AIC")
  varMods2 <- varSelect(data, m, "BIC")
  if(ebic > 0){
    varMods3 <- varSelect(data, m, "EBIC", gamma = ifelse(ebic == 1, .5, .25))
    if(ebic > 1){varMods4 <- varSelect(data, m, "EBIC", gamma = .5)}
  }
  if(cv != FALSE){
    if(cv != TRUE){seed <- cv} else {seed <- NULL}
    varMods5 <- varSelect(data, m, "CV", seed = seed)
  }
  varMods <- list(varMods1, varMods2)
  if(ebic > 0){varMods <- append(varMods, list(varMods3))}
  if(ebic > 1){varMods <- append(varMods, list(varMods4))}
  if(cv != FALSE){varMods <- append(varMods, list(varMods5))}
  names(varMods) <- paste0("varMods", 1:length(varMods))
  list2env(varMods, envir = .GlobalEnv)
  list2env(list(varMods = varMods), envir = .GlobalEnv)
  if(fit == FALSE){
    return(message("varMods in .GlobalEnv"))
  } else {
    if(zero){fit00 <- fitNetwork(data[, -m], ...)}
    fit0 <- fitNetwork(data, covariates = m, ...)
    fit1 <- fitNetwork(data, m, ...)
    fit2 <- fitNetwork(data, m, varMods1, ...)
    fit3 <- fitNetwork(data, m, varMods2, ...)
    if(ebic > 0){fit4 <- fitNetwork(data, m, varMods3, ...)}
    if(ebic > 1){fit5 <- fitNetwork(data, m, varMods4, ...)}
    if(cv != FALSE){
      fit6 <- fitNetwork(data, m, varMods5, ...)
      fit7 <- fitNetwork(data, m, varMods5, which.lam = "1se", ...)
    }
    fits <- list(fit0 = fit0, fit1 = fit1, fitAIC = fit2, fitBIC = fit3)
    if(zero){
      fits <- append(list(fit00), fits)
      names(fits)[1] <- "fit00"
    }
    bb <- which(names(fits) == "fitBIC")
    if(ebic > 0){
      fits <- append(fits, list(fit4))
      names(fits)[-c(1:bb)] <- "fitEBIC"
    }
    if(ebic > 1){
      fits <- append(fits, list(fit5))
      names(fits)[-c(1:bb)] <- c("fitEBIC.25", "fitEBIC.5")
    }
    if(cv != FALSE){
      bb <- length(fits)
      fits <- append(fits, list(fit6))
      fits <- append(fits, list(fit7))
      names(fits)[-c(1:bb)] <- c("fitCVmin", "fitCV1se")
    }
    list2env(fits, envir = .GlobalEnv)
    list2env(list(fits = fits), envir = .GlobalEnv)
    return(message("varMods and fits in .GlobalEnv"))
  }
}


allComparisons <- function(fits, num = FALSE, rev = FALSE){
  zzz <- t(combn(length(fits), 2))
  zzzz <- matrix(NA, ncol = length(fits), nrow = length(fits))
  zzzzz <- list()
  for(comp in 1:nrow(zzz)){
    zzzz[zzz[comp, 1], zzz[comp, 2]] <- ifelse(modLL(fits[[zzz[comp, 1]]], fits[[zzz[comp, 2]]])$decision[2] == "net1", 1, 0)
    zzzzz[[comp]] <- modLL(fits[[zzz[comp, 1]]], fits[[zzz[comp, 2]]], T)$decision
    if(any(zzzzz[[comp]] == "")){zzzzz[[comp]] <- zzzzz[[comp]][-which(zzzzz[[comp]] == "")]}
    if(length(zzzzz[[comp]]) == 0){
      zzzzz[[comp]] <- NA
    } else if(length(zzzzz[[comp]]) == 1){
      zzzzz[[comp]] <- ifelse(zzzzz[[comp]] == ifelse(rev, "net2", "net1"), 1, 0)
    } else if(!ifelse(rev, "net2", "net1") %in% zzzzz[[comp]]){
      zzzzz[[comp]] <- 0
    } else {
      if(num){zzzzz[[comp]] <- sum(zzzzz[[comp]] == ifelse(rev, "net2", "net1"))}
      if(!num){zzzzz[[comp]] <- round(sum(zzzzz[[comp]] == ifelse(rev, "net2", "net1"))/length(zzzzz[[comp]]), 2)}
    }
    zzzz[zzz[comp, 2], zzz[comp, 1]] <- zzzzz[[comp]]
  }
  if(any(is.na(zzzz))){zzzz[is.na(zzzz)] <- "-"}
  diag(zzzz) <- "-"
  zzzz <- data.frame(zzzz)
  dimnames(zzzz) <- rep(list(names(fits)), 2)
  zzzz
}


