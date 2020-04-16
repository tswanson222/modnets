################################################################################
##### makeNet: make Ising network 
makeNet <- function(connectivity, external = NULL, niter = 10000, 
                    situWeights = c(0, -.1, .2, .6), parameters = NULL, 
                    k = 9, p = 0, sitTime = c(300, 1000), plot = TRUE){
  if(length(connectivity) == 1){connectivity <- rep(connectivity, niter)}
  if(!is.null(external)){if(length(external) == 1){external <- rep(external, niter)}}
  if(is.null(parameters)){
    A <- (1/c(7.3, 16, 9.3, 5, 7.5, 5.4, 5.1, 11.1, 9.1)) * 12
    B <- (1/c(18.7, 11.3, 15, 26, 19.6, 24.3, 15.4, 15.5, 13.8)) * 28
  }
  if(p == 0){
    edges <- cbind(c(rep(0, 4), 1, 1:4, 6, 6:8), 
                   c(1, 6, 8, 7, 2, 3, 4, 5, 6, 8, 2, 2, 5))
    edges <- edges + 1
    W <- matrix(0, 9, 9)
    for(i in 1:nrow(edges)){W[edges[i, 1], edges[i, 2]] <- 1}
    W <- W + t(W)
  } else {
    W <- matrix(0, k, k)
    n <- (k * (k - 1))/2
    W[upper.tri(W)] <- if(p < 1){rbinom(n, 1, p)} else {sample(c(rep(0, n - p), rep(1, p)))}
    W <- W + t(W)
  }
  useSitus <- TRUE
  if(is.null(external)){external <- 0} else {useSitus <- FALSE}
  activation <- totalAct <- chanceAct <- matrix(NA, nrow = niter, ncol = k)
  activation[1, ] <- totalAct[1, ] <-  0
  probAct <- function(a, b, Act, shock = 1){(1/(1 + exp(a * (b - Act)))) * shock}
  chanceAct[1, ] <- sapply(1:9, function(z) probAct(A[z], B[z], 0))
  deltaSit <- sample(sitTime[1]:sitTime[2], 1)
  situ <- sample(1:length(situWeights), 1)
  for(i in 2:niter){
    active <- activation[i - 1, ]
    for(j in 1:k){
      totalAct[i, j] <- sum(W[j, -j] * active[-j] * connectivity[i - 1]) + external[i - 1]
      chanceAct[i, j] <- probAct(A[j], B[j], totalAct[i, j])
      activation[i, j] <- ifelse(active[j] == 1, 
                                 ifelse(runif(1) > chanceAct[i, j], 0, 1), 
                                 ifelse(runif(1) < chanceAct[i, j], 1, 0))
    }
    if(i < deltaSit & useSitus){
      external[i] <- rnorm(1, situWeights[situ], sd = .001)
    } else if(useSitus){
      deltaSit <- deltaSit + sample(sitTime[1]:sitTime[2], 1)
      situ <- sample(1:length(situWeights), 1)
      external[i] <- rnorm(1, situWeights[situ], sd = .001)
    }
  }
  data <- data.frame(activation)
  out <- data.frame(external, data)
  if(plot){plot.ts(rowSums(out[,-1]), ylim = c(0, k))}
  out$active <- rowSums(out[, -1])
  out$connect <- connectivity
  out
}

hysteresis <- function(connectivity = 1, niter = 10000, omega = .005, 
                       bin = .05, amp = 3, eps = .02){
  f <- (2 * pi)/omega
  external <- amp * sin(omega * 1:niter) + rnorm(niter, sd = eps)
  out <- makeNet(connectivity = connectivity, niter = niter, 
                 external = external, plot = FALSE)
  stress_bins <- seq(min(out$external), max(out$external), by = bin)
  x1 <- round(f * .25)
  x2 <- round(f * .5)
  fn <- floor(niter/f)
  direction <- c(rep(1, x1), rep(rep(c(0, 1), each = x2), fn), rep(0, x2), rep(1, x1))
  out$direction <- direction[1:niter]
  out$id = 1:niter
  out <- out[order(out$external), ]
  out$bins <- cut(out$external, length(stress_bins), labels = FALSE)
  out$connectivity <- connectivity
  out0 = out[out$direction == 0, ]
  out1 = out[out$direction == 1, ]
  d0 = sapply(split(out0$active, out0$bins), mean)
  d1 = sapply(split(out1$active, out1$bins), mean)
  plot(stress_bins, d0, type = "b", col = "black", lwd = 2, xlab = "Stress", 
       ylab = "Depression", axes = F); axis(1); axis(2); 
  lines(stress_bins, d1, col = "grey50", lwd = 2)
  points(stress_bins, d1, col = "grey50", lwd = 2)
  legend("topleft", legend = c("Decreasing", "Increasing"), bty = "n",
         col = c("black", "grey50"), lwd = c(2, 2), pch = c(1, 1))
  out <- out[order(out$id), ]
  return(out)
}



w <- function(k, p = .5, plot = TRUE, e = NULL){
  n <- (k * (k - 1))/2
  W <- matrix(0, k, k)
  if(is.null(e)){
    W[upper.tri(W)] <- rbinom(n, 1, p)
  } else {
    W[upper.tri(W)] <- sample(c(rep(0, n - e), rep(1, e)))
  }
  W <- W + t(W)
  if(!plot){return(W)}
  plotNet(W, layout = "circle", title = paste("Total Edges =", sum(W)/2))
}


vs <- c("Depressed mood", "Psychomotor change", "Appetite change", 
        "Sleep disturbance", "Anhedonia", "Fatigue", "Suicidal ideation",
        "Worthlessness", "Trouble concentrating")
selves <- c("romantic", "student", "work", "friend")
