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
                 simult = simulVars, call = obj$call)
  return(output)
}


plotStability <- function(obj, pp, s = 3, thresh = .5, color = "black"){
  if("stability" %in% names(obj)){obj <- obj$stability}
  node <- names(obj[[s]])[pp]
  p <- length(obj[[s]]) - ifelse(s == 3, 0, ifelse(!obj$call$exogenous, 2, 1))
  stopifnot(pp <= p + 1)
  n <- nrow(obj[[s]][[pp]])
  k <- ncol(obj[[s]][[pp]])
  plot(0, type = "n", ylim = c(0, 1), xlim = c(1, n + 1), axes = F,
       main = node, xlab = expression(lambda), ylab = "Selection Probability") 
  axis(1); axis(2)
  if(any(obj$call$moderators == 0)){
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

