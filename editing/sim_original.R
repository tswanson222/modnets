################################################################################
##### simNets: simulate datasets and fit models for different sample sizes
simNets <- function(n, niter = 10, tt = 50, p = 3, m = "random", ind = "cor",
                    FUN = "mlGVAR", argsDat = NULL, argsFit = NULL, 
                    verbose = TRUE){
  nit <- paste0("iter", 1:niter)
  A1 <- list(nTime = tt, nPerson = n, nNode = p, m = m)
  if(!is.null(argsDat)){A1 <- append(A1, argsDat)}
  A2 <- list(data = NULL, m = p + 1, verbose = FALSE)
  if(is.null(m)){A2$m <- NULL}
  if(!is.null(argsFit)){A2 <- append(A2, argsFit)}
  A1 <- A1[intersect(formalArgs("mlGVARsim"), names(A1))]
  FUN <- match.arg(FUN, c("mlGVAR", "lmerVAR"))
  FUNargs <- setdiff(formalArgs(FUN), "...")
  A2 <- A2[intersect(FUNargs, names(A2))]
  outs <- suppressWarnings(lapply(seq_along(n), function(z){
    A1$nPerson <- n[z]
    message("Simulating datasets")
    pb1 <- txtProgressBar(max = niter, style = 3)
    x <- setNames(lapply(seq_len(niter), function(i){
      ii <- do.call(mlGVARsim, A1)
      setTxtProgressBar(pb1, i)
      return(ii)
    }), nit)
    close(pb1)
    message("Fitting models")
    pb2 <- txtProgressBar(max = niter, style = 3)
    fits <- setNames(lapply(seq_along(x), function(i){
      A2$data <- switch(2 - is.null(m), x[[i]]$data, x[[i]]$mm$data)
      k <- do.call(match.fun(FUN), A2)
      setTxtProgressBar(pb2, i)
      return(k)
    }), nit)
    close(pb2)
    n1 <- c("sensitivity", "specificity", "precision", "accuracy")
    n2 <- c("temporal", "contemporaneous", "between")
    y1 <- lapply(seq_along(x), function(j){
      compareMods(fits[[j]], x[[j]], ind = ind)})
    yy1 <- do.call(rbind, lapply(y1, '[[', 1))
    y2 <- lapply(seq_along(x), function(j){
      compareMods(fits[[j]], x[[j]], threshold = TRUE, ind = ind)})
    yy2 <- do.call(rbind, lapply(y2, '[[', 1))
    rownames(yy1) <- rownames(yy2) <- 1:niter
    y3 <- lapply(seq_along(x), function(j) performance(fits[[j]], x[[j]]))
    yy3 <- do.call(abind::abind, c(y3, along = 3))
    yy3 <- setNames(lapply(1:3, function(m){
      mmm <- sapply(1:4, function(mm) yy3[m, mm, ])
      colnames(mmm) <- n1
      rownames(mmm) <- 1:niter
      return(mmm)
    }), n2)
    y4 <- lapply(seq_along(x), function(j){
      performance(fits[[j]], x[[j]], threshold = TRUE)})
    yy4 <- do.call(abind::abind, c(y4, along = 3))
    yy4 <- setNames(lapply(1:3, function(m){
      mmm <- sapply(1:4, function(mm) yy4[m, mm, ])
      colnames(mmm) <- n1
      rownames(mmm) <- 1:niter
      return(mmm)
    }), n2)
    out1 <- list(regular = yy1, thresholded = yy2)
    out2 <- list(regular = yy3, thresholded = yy4)
    mains <- list(similarity = out1, performance = out2)
    allNets <- c("kappa", "pcc", "beta", "pdc", "between")
    nets0 <- setNames(lapply(x, function(f1){
      setNames(lapply(allNets, function(f2){
        net(f1, f2)}), allNets)}), nit)
    nets1 <- setNames(lapply(fits, function(f1){
      setNames(lapply(allNets, function(f2){
        net(f1, f2)}), allNets)}), nit)
    nets2 <- setNames(lapply(fits, function(f1){
      setNames(lapply(allNets, function(f2){
        net(f1, f2, threshold = TRUE)}), allNets)}), nit)
    nets <- list(nets0 = nets0, nets1 = nets1, nets2 = nets2)
    if(!is.null(m)){
      for(N in 1:niter){
        nets0[[N]]$ints <- netInts(x[[N]])
        nets1[[N]]$ints <- netInts(fits[[N]])
        nets2[[N]]$ints <- netInts(fits[[N]], threshold = TRUE)
      }
      nets <- list(nets0 = nets0, nets1 = nets1, nets2 = nets2)
      y5 <- t(sapply(seq_along(x), function(w){
        w1 <- matrixDist(netInts(fits[[w]]), netInts(x[[w]]), ind = ind)
        w2 <- matrixDist(netInts(fits[[w]], threshold = T), netInts(x[[w]]), ind = ind)
        return(c(regular = w1, thresholded = w2))
      }))
      rownames(y5) <- 1:niter
      y6 <- lapply(seq_along(x), function(q){
        q1 <- performance(netInts(fits[[q]]), netInts(x[[q]]))
        q2 <- performance(netInts(fits[[q]], threshold = T), netInts(x[[q]]))
        return(list(regular = q1, thresholded = q2))
      })
      out4 <- list(regular = do.call(rbind, lapply(y6, '[[', 1)), 
                   thresholded = do.call(rbind, lapply(y6, '[[', 2)))
      rownames(out4[[1]]) <- rownames(out4[[2]]) <- 1:niter
      ints <- list(similarity = y5, performance = out4)
    }
    output <- list(x = x, nets = nets, mains = mains)
    if(!is.null(m)){output$ints <- ints}
    return(output)
  }))
  return(outs)
}



