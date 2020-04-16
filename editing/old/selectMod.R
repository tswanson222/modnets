selectVars <- function(obj, select = "select", type = "gaussian"){
  allNames <- lapply(lapply(obj, '[[', "predictor"), as.character)
  vs <- names(allNames) <- names(obj)
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
  return(varMods)
}


getFitCIs <- function(fit, allNames){
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

