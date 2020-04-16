out <- readRDS('~/Desktop/currentSims/four/outs_n200.RDS')
out0 <- evalFits(out1 = out, rmFailed = T, threshold = T)


rule <- 'or'
out <- setNames(lapply(dir(), readRDS), gsub('.*_|.RDS$', '', dir()))
out0 <- lapply(out, evalFits, rmFailed = TRUE, threshold = TRUE)
out00 <- setNames(lapply(seq_along(out0), function(i) setNames(lapply(
  split(out0[[i]]$x0$out2[[rule]], out0[[i]]$x0$out2[[rule]]$type), function(j){
    data.frame(do.call(cbind, split(j$value, j$measure)))[, 1:4]
  }), names(out3[[1]][[1]]))), names(out))


out1 <- lapply(out, function(z) z[-which(sapply(z, length) != 8)])
all(sapply(seq_along(out1), function(z) identical(out0[[z]][[1]][[1]], out1[[z]])))
out2 <- lapply(out1, function(z){
  lapply(z, lapply, netInts, rule = rule, threshold = TRUE, avg = F, empty = FALSE)
})
out3 <- lapply(out2, function(z) lapply(seq_along(z), function(i){lapply(
  z[[i]][-1], function(j){as.matrix(unlist(CompareNetworks(
    true = z[[i]]$x, est = j, directed = FALSE)))})}))
out4 <- lapply(out3, function(z) setNames(lapply(seq_len(unique(sapply(z, length))), function(i){
  data.frame(t(do.call(cbind, lapply(z, '[[', i))))[, 1:4]
}), names(out3[[1]][[1]])))


sapply(seq_along(out4), function(i) sapply(seq_along(out4[[i]]), function(j) identical(out4[[i]][[j]], out00[[i]][[j]])))


sapply(seq_along(out4), function(z) identical(out4[[z]], out00[[z]]))


################################################################################
################################################################################
bias <- function(x,y) mean(abs(x-y),na.rm=TRUE)

MSE <- function(x,y) mean((x-y)^2,na.rm=TRUE)


cor0 <- function(x,y){
  if (all(is.na(x)) || all(is.na(y))){
    return(NA)
  }
  
  if (sd(x,na.rm=TRUE)==0 | sd(y,na.rm=TRUE) == 0){
    return(0)
  }
  
  return(cor(x,y,use="pairwise.complete.obs"))
}


# Function to compute sensitivity and specificity:
CompareNetworks <- function(true,est, directed = TRUE){
  if (is.matrix(true) & is.matrix(est)){
    if (directed){
      real <- c(true)
      est <- c(est)
    } else {
      real <- true[upper.tri(true,diag=FALSE)]
      est <- est[upper.tri(est,diag=FALSE)]
    }        
  } else {
    real <- true
  }
  
  # True positives:
  TruePos <- sum(est != 0 &  real != 0)
  
  # False pos:
  FalsePos <- sum(est != 0 & real == 0)
  
  # True Neg:
  TrueNeg <- sum(est == 0 & real == 0)
  
  # False Neg:
  FalseNeg <- sum(est == 0 & real != 0)
  
  out <- list()
  
  ### Sensitivity:
  out$Sensitivity <- TruePos / (TruePos + FalseNeg)
  
  # Specificity:
  out$Specificity <- TrueNeg / (TrueNeg + FalsePos)
  
  # Precision:
  out$Precision <- TruePos / (TruePos + FalsePos)
  
  # Accuracy:
  out$Accuracy <- (TruePos + TrueNeg) / (TruePos + TrueNeg + FalsePos + FalseNeg)
  
  # Correlation:
  out$Correlation <- cor0(est,real)
  
  out$MAE <- bias(est,real)
  
  out$MSE <- MSE(est,real)
  
  return(out)
}

################################################################################
################################################################################
combineSamples <- function(x){
  o23 <- c('out2', 'out3')
  N <- gsub('[^0-9]', '', names(x))
  x0 <- sapply(x, length)
  if(all(x0 == 2)){
    out0 <- setNames(lapply(o23, function(z) lapply(lapply(x, '[[', 'x0'), '[[', z)), o23)
    out1 <- setNames(lapply(o23, function(z) lapply(lapply(x, '[[', 'x1'), '[[', z)), o23)
    
  }
}

