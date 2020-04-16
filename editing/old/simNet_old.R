##### simNet: cross-sectional simulations (OLD)
simNet <- function(N = 100, p = 5, m = FALSE, m2 = .1, b1 = NULL, b2 = NULL, 
                   sparsity = .5, intercepts = NULL, nIter = 250, msym = TRUE, 
                   onlyDat = FALSE, pbar = TRUE, div = 1000, gibbs = TRUE,
                   ordinal = FALSE, nLevels = 5, mord = FALSE, 
                   mbinary = FALSE, minOrd = 3){
  m2_range <- c(.1, .3)
  if(is.null(intercepts)){intercepts <- rep(0, p)}
  if(mord & FALSE){
    mfun <- function(prob = NULL, x = 1:nLevels, size = 1){
      if(!is.null(prob)){
        if(length(prob) == 1){
          prob <- abs(rnorm(n = length(x), mean = prob, sd = 1))
          prob <- prob/sum(prob)
        }
      }
      sample(x = x, size = size, prob = prob)
    }
  } else if(mbinary & FALSE){
    mfun <- function(prob = NULL, x = 0:1, size = 1){
      sample(x = x, size = size, prob = prob)
    }
  } else {
    mfun <- function(mean = 0, n = 1, sd = 1){rnorm(n = n, mean = mean, sd = sd)}
  }
  if(isTRUE(m)){m <- mfun()}
  m <- ifelse(identical(m, FALSE), 0, m)
  if(is.null(b1)){b1 <- simPcor(p, sparsity, parRange = c(.1, .4))}
  if(is.null(b2)){
    b2 <- matrix(sample(c(0, 1), p^2, TRUE, prob = c(1 - m2, m2)), p, p)
    b2 <- b2 * sample(c(-1, 1), p^2, TRUE, prob = c(1 - .5, .5))
    b2 <- matrix(b2 * runif(p^2, min(m2_range), max(m2_range)), ncol = p)
    if(all(m == 0)){b2 <- diag(0, p)}
    if(msym){b2 <- simPcor(p, sparsity = 1 - m2, parRange = m2_range)}
    if(max(abs(b2)) > max(m2_range)){b2 <- b2/(max(abs(b2))/max(m2_range))}
  }
  diag(b1) <- diag(b2) <- 0
  if(!all(m == 0) | isTRUE(gibbs)){
    data <- matrix(NA, nrow = N, ncol = p)
    M <- c()
    if(pbar){pb <- txtProgressBar(max = N, style = 3)}
    for(case in 1:N){
      sampling <- matrix(NA, nrow = nIter, ncol = p)
      if(!ordinal | TRUE){
        sampling[1, ] <- rnorm(p)
        m0 <- c(ifelse(m == 0, 0, mfun()), numeric(nIter - 1))
        for(iter in 2:nIter){
          m0[iter] <- ifelse(m == 0, 0, mfun(m))
          for(v in 1:p){
            v_mu <- 0
            v_ps <- which(b1[v, ] != 0)
            if(length(v_ps) > 0){
              for(vp in v_ps){
                v_it <- iter - as.numeric(vp >= v)
                v_mu <- c(v_mu, sampling[v_it, vp] * b1[v, vp])
              }
            }
            v_ps2 <- which(b2[v, ] != 0)
            if(length(v_ps2) > 0){
              for(vp in v_ps2){
                v_it <- iter - as.numeric(vp >= v)
                v_mu <- c(v_mu, ((sampling[v_it, vp] * m0[v_it]) * b2[v, vp]))
              }
            }
            v_mu <- intercepts[v] + sum(v_mu)
            sampling[iter, v] <- rnorm(1, v_mu, 1)
            if(any(abs(sampling[!is.na(sampling)]) > div)){stop('Sampler diverged')}
          }
        }
      } else { ## NOT WORKING
        sampling[1, ] <- sample(1:nLevels, 1)
        m0 <- c(ifelse(m == 0, 0, mfun()), numeric(nIter - 1)) # FIX
        for(iter in 2:nIter){
          m0[iter] <- ifelse(m == 0, 0, mfun(m)) # FIX
          for(v in 1:p){
            potentials <- rep(NA, nLevels)
            for(cat in 1:nLevels){
              v_mu <- 0
              v_ps <- which(b1[v, ] != 0)
              if(length(v_ps) > 0){
                for(vp in v_ps){
                  v_it <- iter - as.numeric(vp >= v)
                  v_mu <- c(v_mu, sampling[v_it, vp] * b1[v, vp])
                }
              }
              v_ps2 <- which(b2[v, ] != 0)
              if(length(v_ps2) > 0){
                for(vp in v_ps2){
                  v_it <- iter - as.numeric(vp >= v)
                  v_mu <- c(v_mu, ((sampling[v_it, vp] * m0[v_it]) * b2[v, vp]))
                }
              }
              potentials[cat] <- intercepts[v] + sum(v_mu)
            }
            probs <- exp(potentials)/sum(exp(potentials))
            sampling[iter, v] <- sample(1:nLevels, 1, prob = probs)
          }
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
    data <- mvtnorm::rmvnorm(N, sigma = Sigma)
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
    if(mord){
      mord <- c()
      while(length(unique(mord)) < minOrd){
        mord <- as.numeric(cut(M, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
      }
    }
    data <- cbind(data, M)
  }
  out <- append(list(data = data.frame(data)), out)
  if(onlyDat){out <- data.frame(data)}
  if(!all(m == 0)){attributes(out)[c('m', 'm2')] <- list(m, m2)}
  class(out) <- c(ifelse(onlyDat, 'data.frame', 'list'), 'mgmSim')
  return(out)
}
