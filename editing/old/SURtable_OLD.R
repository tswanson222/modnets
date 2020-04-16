##### SURtable_OLD: print fit indices for 3 SUR models
SURtable_OLD <- function(fits, ind = "RMSE", verbose = TRUE, d = 3,
                         labs = c("fitOLS", "fit0", "fit1se")){
  stopifnot(class(fits) == "list")
  if(all(sapply(lapply(fits, names), function(z) "SURfit" %in% z))){
    lrdat <- t(sapply(fits, SURll))[, 1:2]
    colnames(lrdat)[1] <- "logLik"
    fits <- lapply(fits, function(z) z$SURfit)
  }
  if(ind %in% c("RSS", "rss")){ind <- "SSR"}
  ind <- match.arg(ind, choices = c("RMSE", "MSE", "SSR", "R2", "adjR2", "omnibus", "lrtest"))
  measure <- switch(ind, "RMSE" = "sigma", "MSE" = "sigma", "R2" = "r.squared", "lrtest" = "lrtest",
                    "adjR2" = "adj.r.squared", "SSR" = "ssr", "omnibus" = "omnibus")
  if(is.null(names(fits))){
    if(all(labs == c("fitOLS", "fit0", "fit1se"))){
      names(fits) <- labs[rev(order(sapply(fits, function(z) z$rank)))]
    } else {
      names(fits) <- labs
    }
  }
  p <- length(fits[[1]]$eq)
  if(measure == "omnibus"){
    dat <- matrix(NA, ncol = 8, nrow = 3)
    dat[,1] <- sapply(fits, logLik)
    dat[,2] <- sapply(fits, function(z) summary(z)$df[1] + (p * (p + 1))/2)
    dat[,3] <- sapply(fits, AIC)
    dat[,4] <- sapply(fits, BIC)
    dat[,5] <- sapply(fits, function(z) sum(colSums(residuals(z)^2)))
    dat[,6] <- sapply(fits, function(z) summary(z)$detResidCov)
    dat[,7] <- sapply(fits, function(z) summary(z)$ols.r.squared)
    dat[,8] <- sapply(fits, function(z) summary(z)$mcelroy.r.squared)
    colnames(dat) <- c("LL", "df", "AIC", "BIC", "SSR", "detSigma", "OLS.R2", "McElroy.R2")
    rownames(dat) <- names(fits)
  } else if(measure == "lrtest"){
    lrdat <- lrdat[order(lrdat[,"df"], decreasing = TRUE),]
    x <- t(combn(1:3, 2))
    dat <- data.frame(LL1 = rownames(lrdat)[x[,1]], LL2 = rownames(lrdat)[x[,2]], 
                      Chisq = NA, Df = NA, p.value = NA, sig = NA)
    for(i in 1:3){
      dat[i,3] <- (2 * (lrdat[x[i,1], 1] - lrdat[x[i,2], 1]))
      dat[i,4] <- abs(lrdat[x[i,1], 2] - lrdat[x[i,2], 2])
      dat[i,5] <- pchisq(dat[i,3], dat[i,4], lower.tail = FALSE)
      dat[i,6] <- ifelse(dat[i,5] < .001, "***", 
                         ifelse(dat[i,5] < .01, "**", 
                                ifelse(dat[i,5] < .05, "*", 
                                       ifelse(dat[i,5] < .1, ".", ""))))
    }
    dat[,c(3,5)] <- round(dat[,c(3,5)], digits = d)
    RMSEA <- function(X2, df, N){
      rmsea <- sqrt(max(c(((X2/N)/df) - (1/N), 0)))
      lower.l <- function(lambda){(pchisq(X2, df = df, ncp = lambda) - .95)}
      lambda.l <- tryCatch({uniroot(f = lower.l, lower = 0, upper = X2)$root}, error = function(e){0})
      rmsea.lower <- sqrt(lambda.l/(N * df))
      upper.l <- function(lambda){(pchisq(X2, df = df, ncp = lambda) - .05)}
      lambda.u <- tryCatch({uniroot(f = upper.l, lower = 0, upper = max(N, X2 * 4))$root}, error = function(e){1})
      rmsea.upper <- sqrt(lambda.u/(N * df))
      rmsea.pvalue <- 1 - pchisq(X2, df = df, ncp = (N * df * (.05^2)))
      out <- c(lower = rmsea.lower, RMSEA = rmsea, upper = rmsea.upper, p.value = rmsea.pvalue)
      out
    }
    N <- fits[[1]]$df.residual + fits[[1]]$rank
    rmsea <- round(t(mapply(RMSEA, dat$Chisq, dat$Df, N)), 3)
    rownames(rmsea) <- rownames(dat)
    lrdat <- cbind(lrdat, AIC = (-2 * lrdat[,1]) + (2 * lrdat[,2]), 
                   BIC = (-2 * lrdat[,1]) + (lrdat[,2] * log(N)))
    dat <- list(LRtest = dat, RMSEA = rmsea, LLs = lrdat)
    if(verbose == TRUE){
      cat("\nIf difference is NOT significant, 'LL2' is preferred\nIf difference IS significant, 'LL1' is preferred\n\n")}
  } else {
    fs <- lapply(fits, function(z) lapply(z$eq, summary))
    dat <- sapply(fs, function(z) sapply(z, function(y) y[[measure]]))
    if(ind == "MSE"){dat <- dat^2}
    rownames(dat) <- sapply(1:p, function(z) as.character(fits[[1]]$eq[[z]]$terms[[2]]))
    colnames(dat) <- paste0(colnames(dat), ".", ind)
  }
  if(measure != "lrtest"){dat <- round(dat, digits = d)}
  dat
}
