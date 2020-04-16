##### modLL: log-likelihood of whole network, & LRT comparing two networks
modLL <- function(net0, net1 = NULL, nodemods = FALSE, all = FALSE, 
                  ftest = FALSE, d = 4, alpha = .05, omnibus = NULL){
  if("fits" %in% names(net0)){
    net1 <- net0$fits$fit1; net0 <- net0$fits$fit0; all <- TRUE
  } else if(!"fitobj" %in% names(net0)){omnibus <- net0}
  if(class(nodemods) == "character"){ftest <- switch(match.arg(nodemods, c("lrt", "ftest")), lrt = FALSE, ftest = TRUE)}
  mll <- function(object){
    k <- length(object$fitobj)
    res <- as.matrix(do.call(cbind, lapply(object$fitobj, resid)))
    n <- nrow(res)
    sigma <- (t(res) %*% res)/n
    inv <- solve(sigma)
    ll <- sum(((-k/2) * log(2 * pi)) - (.5 * log(det(sigma))) - (.5 * diag(res %*% inv %*% t(res))))
    df <- sum(sapply(lapply(object$mods, '[[', "model"), nrow)) + 1
    return(c(LL = ll, df = df))
  }
  if(is.null(net1) & is.null(omnibus)){omnibus <- TRUE}
  if(!is.null(omnibus) & !all){
    if(class(omnibus) != "list"){
      k1 <- length(net0$fitobj)
      n1 <- nrow(net0$data) * k1
      ll1 <- mll(net0)
      aic1 <- (2 * ll1[2]) - (2 * ll1[1])
      bic1 <- (ll1[2] * log(n1)) - (2 * ll1[1])
      out1 <- c(ll1, AIC = unname(aic1), BIC = unname(bic1))
      if(is.null(net1)){return(out1)}
      k2 <- length(net1$fitobj)
      n2 <- nrow(net1$data) * k2
      ll2 <- mll(net1)
      aic2 <- (2 * ll2[2]) - (2 * ll2[1])
      bic2 <- (ll2[2] * log(n2)) - (2 * ll2[1])
      out2 <- c(ll2, AIC = unname(aic2), BIC = unname(bic2))
      out <- rbind(out1, out2)
      rownames(out) <- c("net0", "net1")
      return(out)
    } else {
      nets <- length(omnibus)
      k <- sapply(omnibus, function(z) length(z$fitobj))
      n0 <- sapply(omnibus, function(z) nrow(z$data))
      n <- n0 * k
      out <- data.frame(t(sapply(omnibus, mll)))
      out$AIC <- sapply(1:nets, function(z) (2 * out[z, 2]) - (2 * out[z, 1]))
      out$BIC <- sapply(1:nets, function(z) (log(n[z]) * out[z, 2]) - (2 * out[z, 1]))
      rownames(out) <- ifelse(!is.null(names(omnibus)), list(names(omnibus)), list(paste0("net", 0:(nets - 1))))[[1]]
      return(out)
    }
  } else {
    if(is.null(net1) & !is.null(omnibus)){net1 <- net0[[2]]; net0 <- net0[[1]]}
    fitnames <- ifelse(!is.null(omnibus) & all == TRUE, 
                       list(names(omnibus)[1:2]), list(paste0("net", 0:1)))[[1]]
    k1 <- length(net0$fitobj); n1 <- nrow(net0$data) * k1
    k2 <- length(net1$fitobj); n2 <- nrow(net1$data) * k2
    ll1 <- unname(mll(net0)); ll2 <- unname(mll(net1))
    abic <- function(LL, n = NULL, ab = FALSE){
      if(ab){
        aic <- 2 * LL[2] - 2 * LL[1]
        bic <- log(n) * LL[2] - 2 * LL[1]
        return(c(aic, bic))
      }
      ifelse(is.null(n), 2, log(n)) * LL[2] - 2 * LL[1]
    }
    lls <- data.frame(matrix("", ncol = 8, nrow = 2))
    maxdf <- which.max(c(ll1[2], ll2[2]))
    lls[, 1:2] <- do.call(rbind, ifelse(maxdf == 1, list(list(ll2, ll1)), list(list(ll1, ll2)))[[1]])
    rownames(lls) <- ifelse(maxdf == 1, list(rev(fitnames)), list(fitnames))[[1]]
    lldiff <- abs(ll1[1] - ll2[1]) * 2
    dfdiff <- abs(ll1[2] - ll2[2])
    pval <- pchisq(lldiff, dfdiff, lower.tail = FALSE)
    decision <- ifelse(pval <= alpha, rownames(lls)[which.max(lls[, 1])], 
                       rownames(lls)[which.min(lls[, 2])])
    out <- c(lldiff, dfdiff, pval)
    if(!is.null(d)){out <- round(out, d)}
    lls[2, 3:6] <- c(out, decision)
    lls[ , 7:8] <- t(do.call(
      cbind, ifelse(
        maxdf == 1, list(list(abic(ll2, n2, T), abic(ll1, n1, T))), 
        list(list(abic(ll1, n1, T), abic(ll2, n2, T))))[[1]]))
    colnames(lls) <- c("LL", "df", "LL_diff2", "Df_diff", "pval", "decision", "AIC", "BIC")
    if(all(c(nodemods, all, ftest) == FALSE)){return(lls)}
    netLRT <- function(net0, net1, d = NULL, all = FALSE, alpha = .05, ftest = FALSE){
      ss1 <- unlist(sapply(net0$mods, '[', "deviance"))
      ss2 <- unlist(sapply(net1$mods, '[', "deviance"))
      ll1 <- unlist(sapply(net0$mods, '[', "LL_model"))
      ll2 <- unlist(sapply(net1$mods, '[', "LL_model"))
      df1 <- unlist(sapply(net0$fitobj, '[', "df.residual"))
      df2 <- unlist(sapply(net1$fitobj, '[', "df.residual"))
      lldiff <- abs(ll1 - ll2) * 2
      dfdiff <- abs(df1 - df2)
      ssdiff <- abs(ss1 - ss2)
      if(ftest){
        dfs <- apply(cbind(df1, df2), 1, which.min)
        fullmods <- sapply(1:length(dfs), function(z){
          unname(cbind(ss1, ss2)[z, ][dfs[z]]/cbind(df1, df2)[z, ][dfs[z]])
        })
        dffull <- sapply(1:length(dfs), function(z) cbind(df1, df2)[z, ][dfs[z]])
        fchange <- (ssdiff/dfdiff)/fullmods
        if(any(is.na(fchange))){fchange[is.na(fchange)] <- 0}
        ps <- suppressWarnings(pf(fchange, dfdiff, dffull, lower.tail = FALSE))
        if(any(is.na(ps))){ps[is.na(ps)] <- 1}
        names(fchange) <- names(dfdiff) <- names(ps) <- names(dffull) <- names(net0$mods)
        decision <- c()
        for(i in 1:length(ps)){
          if(ps[i] <= alpha){
            decision[i] <- paste0("net", ifelse(ss1[i] < ss2[i], 0, 1))
          } else if(ps[i] == 1){
            decision[i] <- ""
          } else {
            decision[i] <- paste0("net", ifelse(df1[i] > df2[i], 0, 1))
          }
        }
        out <- list(Fchange = fchange, Df_diff = dfdiff, Df_res = dffull, pval = ps, decision = decision)
      } else {
        ps <- pchisq(q = lldiff, df = dfdiff, lower.tail = FALSE)
        names(lldiff) <- names(dfdiff) <- names(ps) <- names(net0$mods)
        decision <- c()
        for(i in 1:length(ps)){
          if(ps[i] <= alpha){
            decision[i] <- paste0("net", ifelse(ll1[i] > ll2[i], 0, 1))
          } else if(ps[i] == 1){
            decision[i] <- ""
          } else {
            decision[i] <- paste0("net", ifelse(df1[i] > df2[i], 0, 1))
          }
        }
        out <- list(LL_diff2 = lldiff, Df_diff = dfdiff, pval = ps, decision = decision)
      }
      if(!is.null(d)){out$pval <- round(out$pval, d)}
      out <- do.call(cbind.data.frame, out)
      if(all){
        aic1 <- unlist(sapply(net0$mods, '[', "AIC"))
        aic2 <- unlist(sapply(net1$mods, '[', "AIC"))
        bic1 <- unlist(sapply(net0$mods, '[', "BIC"))
        bic2 <- unlist(sapply(net1$mods, '[', "BIC"))
        out1 <- list(RSS = ss1, LL = ll1, df = df1, AIC = aic1, BIC = bic1)
        out2 <- list(RSS = ss2, LL = ll2, df = df2, AIC = aic2, BIC = bic2)
        out1 <- do.call(cbind.data.frame, out1)
        out2 <- do.call(cbind.data.frame, out2)
        rownames(out1) <- rownames(out2) <- rownames(out)
        out3 <- list(net0 = out1, net1 = out2, LRT = out)
        names(out3)[3] <- ifelse(ftest, "Ftest", "LRT")
        return(out3)
      } else {
        return(out)
      }
    }
    if(all){
      univariate <- netLRT(net0, net1, d, all, alpha, ftest = FALSE)
      univariate$Ftest <- netLRT(net0, net1, d, all = FALSE, alpha, ftest = TRUE)
      return(append(univariate, list(omnibus = lls)))
    } else {
      return(netLRT(net0, net1, d, all, alpha, ftest))
    }
  }
}
