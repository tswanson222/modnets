################################################################################
##### Variable selection --> Model fitting procedure outline
source("functions.R")
settings("gauss1")
if(is.null(colnames(data))){colnames(data) <- paste0("V", 1:ncol(data))}
dat <- datMat(data, type, m = 1)
q <- makeEqs(dat = dat)
fit <- systemfit(formula = q, method = "SUR", data = dat$full)
coefs <- names(coef(fit)) # Fitting the model just to get coef names

##### Variable Selection
# 1) LASSO
net1 <- fitNetwork(data, type, lags = 1, seed = 333)
adj <- as.vector(t(net1$adjMat))
N <- list()
for(i in 1:ncol(dat$Y)){N[[i]] <- paste0("eq", i, "_", colnames(data))}
coefs <- data.frame(b = unlist(N), value = adj)

# 2) Best subset selection (exhaustive search)
regMods <- list()
ind <- "cp"
for(i in 1:ncol(dat$Y)){
  regMods[[i]] <- summary(regsubsets(dat$X, dat$Y[,i], nvmax = ncol(dat$X)))
  regMods[[i]] <- ifelse(regMods[[i]]$which, 1, 0)[which.min(regMods[[i]][[ind]]), -1]
}
adj <- as.vector(t(do.call(rbind, regMods)))
N <- list()
for(i in 1:ncol(dat$Y)){N[[i]] <- paste0("eq", i, "_", names(regMods[[1]]))}
coefs <- data.frame(b = unlist(N), value = adj)

##### Fit restricted SUR model
fit2 <- systemfit(formula = q, method = "SUR", data = dat$full,
                  restrict.matrix = coefs[coefs$value == 0, "b"])

# Can fit best subset selection with different indices, although Cp appears to be popular
# Forward/backward selection (stepwise) seem to be the only other routes with unbiased estimators
# 'selectiveInference' can be used for LASSO
# Principal components regression ('pls::pcr') using 10-fold CV
# Partial least squares ('pls::plsr') using 10-fold CV
# Problem is referred to as 'inference after variable selection'

### NEED THIS FOR CLUSTER SCRIPTS THAT USE MODNETS!
.keepn <- new.env()
.keepn$n <- n
attach(.keepn)
clear <- TRUE

################################################################################
###################### Example inputs for the SURsampler #######################
################################################################################
modnets()

##### Generate sample covariance matrices using the Wishart distribution
# MCMCpack::rwish | MCMCpack::dwish | MCMCpack::riwish | MCMCpack::diwish
# tangentially, matrixcalc::is.positive.definite is a nice handy function
### 1) rwish(df, S) will return a symmetric, positive-definite matrix 
# that resembles (S*df)
################### Three Predictors, Two Interaction Terms ####################
# Step 1: Inputs and data generation
n <- 1000
S <- matrix(c(1, .4, .3, .4, 1, .6, .3, .6, 1), 3, 3)
beta <- matrix(c(.01, .014, 0, .005, .07, .3, .01, .09, .1, .23, .3, .001), 
               ncol = 4, nrow = 3, byrow = TRUE)
beta2 <- matrix(c(.01, .05, .1, .1, .01, .005), ncol = 2, byrow = T)
B <- cbind(beta, beta2)
dat <- SURsampler(B = B, S = S, n = n, seed = 666, time = T)

# Step 2: Fit unconstrained model
fit <- SURfit(dat)
net <- SURnet(fit, dat)
results(B, fit); results(S, fit)

# betas & residualCors identical to OLS; residualCovs differ
m <- list(); res <- list()
for(i in 1:length(q)){
  m[[i]] <- coef(lm(q[[i]], data.frame(sampDat$full)))
  res[[i]] <- resid(lm(q[[i]], data.frame(sampDat$full)))
}
m <- do.call(rbind, m); res <- do.call(cbind, res)
sparsify(m, 2); sparsify(beta, 2)
sparsify(cov(res), 1); S

# Step 3: Variable selection and constrained models
vars <- varSelect(dat = dat, nlam = 10)
fit1 <- SURfit(dat = dat, varMods = vars, mod = "min")
fit2 <- SURfit(dat, vars, "1se")
net1 <- SURnet(fit1, dat)
net2 <- SURnet(fit2, dat)
results(B, fit2, labels = fit); results(S, fit2)
results(B, fit1, labels = fit); results(S, fit1)

# Step 4: Model comparison
SURtable(fit, fit1, fit2, "RMSE")
SURtable(fit, fit1, fit2, "adjR2")
SURtable(fit, fit1, fit2, "omni")
lrtest(fit2, fit1, fit)

################################################################################
################### Four Predictors, Three Interaction Terms ###################
################################################################################
modnets()
n <- 500
B <- list()
B[[1]] <- c(0, .3, 0, .25, 0, 0, .07, 0)
B[[2]] <- c(0, .4, 0, 0, .15, 0, 0, .09)
B[[3]] <- c(0, .18, .27, 0, 0, .08, 0, 0)
B[[4]] <- c(0, .5, 0, .31, 0, 0, .1, 0)
B <- do.call(rbind, B)
S <- matrix(c(1, .3, .7, .2, .3, 1, 0, .06, .7, 0, 1, .4, .2, .06, .4, 1), ncol = 4, byrow = T)
dat <- SURsampler(S = S, B = B, n = n, seed = 666, time = T)
vars <- varSelect(dat = dat, seed = 666, nlam = 20)
vars1 <- varSelect(dat = dat, method = "glinternet", nlam = 50)
vars2 <- varSelect(dat = dat, method = "glinternet", m = 1, nlam = 50)
vars3 <- varSelect(dat = dat, method = "glinternet", nlam = 50, useSE = F)
vars4 <- varSelect(dat = dat, method = "glinternet", nlam = 50, m = 1, useSE = F)

res1 <- resample(dat, 50, "glinternet", "split", nfolds = 10, m = 1)
res2 <- resample(dat, 50, "glinternet", "bootstrap", nfolds = 10)

fit <- SURfit(dat = dat)
fit1 <- SURfit(dat = dat, varMods = vars, mod = "min")
fit2 <- SURfit(dat = dat, varMods = vars, mod = "1se")

net <- SURnet(fit, dat)
net1 <- SURnet(fit = fit1, dat = dat)
net2 <- SURnet(fit = fit2, dat = dat)

results(B, fit); results(S, fit)
results(B, fit1); results(S, fit1)
results(B, fit2); results(S, fit2)

fits <- list(fit, fit1, fit2)
SURtable(fits, "RMSE")
SURtable(fits, "adjR2")
SURtable(fits, "omni")
SURtable(fits, "lrt")


################################################################################
modnets()
n <- 500
B <- list()
B[[1]] <- c(0, .4, 0, .25, 0, 0, .07, 0)
B[[2]] <- c(0, .3, 0, 0, .15, 0, 0, .09)
B[[3]] <- c(0, .35, .27, 0, 0, .08, 0, 0)
B[[4]] <- c(0, .5, 0, .31, 0, 0, .1, 0)
B <- do.call(rbind, B)
S <- matrix(c(1, .3, .7, .2, .3, 1, 0, .06, .7, 0, 1, .4, .2, .06, .4, 1), ncol = 4, byrow = T)
dat <- SURsampler(S = S, B = B, n = n, seed = 666, time = T)
res1 <- resample(dat, 100, "glinternet", "split", 666, m = 1)
res2 <- resample(dat, 100, "glinternet", "bootstrap", 666, m = 1)

##################### Two Predictors, One Interaction Term #####################
##### THIS MODEL DOESN'T WORK!!
modnets()
n <- 900
S <- matrix(c(1, .4, .4, 1), 2, 2)
beta <- matrix(c(0, .14, .27, 0, .07, .3), ncol = 3, nrow = 2, byrow = T)
beta2 <- matrix(c(.1, .05), ncol = 1)
B <- cbind(beta, beta2)
dat <- SURsampler(S = S, B = B, n = n, seed = 666)
fit <- systemfit(formula = makeEqs(dat), method = "SUR", data = dat$full)
results(B, fit); results(S, fit)
################################################################################
################################################################################
##### 3D plotting
invisible(suppressMessages(lapply(c("plotly", "cusp"), require, character.only = TRUE)))
a <- seq(-4, 4, by = .1)
b <- seq(-4, 4, by = .1)
ab <- expand.grid(a, b)

z0 <- array()
for(i in 1:nrow(ab)){z0[i] <- dcusp(ab[i, 1], 0, ab[i, 2])}
z0 <- matrix(z0, nrow = length(a), ncol = length(b))

a0 <- plot_ly(x = a, y = b, z = z0) %>% add_surface()
a0p <- persp(x = a, y = b, z = z0, theta = 150)

################################################################################
################################################################################
##### Breaking down SUR in 'systemfit'
data <- data[1:100,]
dat <- datMat(data, type)
q <- makeEqs(dat = dat)
X <- cbind(1, dat$X)
Y <- dat$Y

beta <- solve(t(X) %*% X) %*% t(X) %*% Y
S <- 1/(nrow(X) - 1) * (t(Y - X %*% beta) %*% (Y - X %*% beta))

I <- diag(nrow = nrow(X))
omega <- solve(S) %x% I
e <- matrix(0, ncol = 6, nrow = nrow(X))
x1 <- rbind(X, e, e, e, e)
x2 <- rbind(e, X, e, e, e)
x3 <- rbind(e, e, X, e, e)
x4 <- rbind(e, e, e, X, e)
x5 <- rbind(e, e, e, e, X)
XX <- cbind(x1, x2, x3, x4, x5)

beta_1 <- solve(t(XX) %*% omega %*% XX) %*% t(XX) %*% omega %*% as.vector(Y)
  
control <- systemfit.control()
modelFrame <- systemfit:::.prepareData.systemfit(dat$full)
modelFrameEq <- list(); evalModelFrameEq <- list()
termsEq <- list(); yVecEq <- list(); xMatEq <- list()
nEq <- length(q)
nCoefEq <- numeric(nEq)
for(i in 1:length(q)){
  modelFrameEq[[i]] <- modelFrame
  modelFrameEq[[i]]$formula <- q[[i]]
  evalModelFrameEq[[i]] <- eval(modelFrameEq[[i]])
  termsEq[[i]] <- attr(evalModelFrameEq[[i]], "terms")
  yVecEq[[i]] <- model.extract(evalModelFrameEq[[i]], "response")
  xMatEq[[i]] <- model.matrix(termsEq[[i]], evalModelFrameEq[[i]])
  nCoefEq[i] <- ncol(xMatEq[[i]])
}

nCoefLiEq <- nCoefEq
nObsWithNa <- length(yVecEq[[1]])
validObsEq <- matrix(NA, nrow = nObsWithNa, ncol = nEq)
for(i in 1:nEq){validObsEq[, i] <- !is.na(yVecEq[[i]]) & rowSums(is.na(xMatEq[[i]])) == 0}
yVecAll <- matrix(0, 0, 1)
for(i in 1:length(q)){yVecAll <- c(yVecAll, yVecEq[[i]])}
xMatAll <- systemfit:::.stackMatList(xMatEq, way = "diag", useMatrix = TRUE)


coef <- solve(crossprod(xMatAll), crossprod(xMatAll, yVecAll), tol = control$solvetol)

if(method %in% c("SUR")){
  coefOld <- coef
  coefDiff <- coef
  iter <- 0
  while((sum(coefDiff^2)/sum(coefOld^2))^0.5 > control$tol & iter < control$maxiter){
    iter <- iter + 1
    coefOld <- coef
    resids <- yVecAll - xMatAll %*% coef
    rcov <- systemfit:::.calcResidCov(resids, methodResidCov = control$methodResidCov, 
                                      validObsEq = validObsEq, nCoefEq = nCoefLiEq, 
                                      xEq = xMatEq, centered = control$centerResiduals, 
                                      useMatrix = control$useMatrix, solvetol = control$solvetol)
    coef <- systemfit:::.calcGLS(xMat = xMatAll, yVec = yVecAll, 
                                 R.restr = R.restr, q.restr = q.restr, sigma = rcov, 
                                 validObsEq = validObsEq, useMatrix = control$useMatrix, 
                                 solvetol = control$solvetol)
    coefDiff <- coef - coefOld
  }
  coefCov <- .calcGLS(xMat = xMatAll, R.restr = R.restr, 
                      q.restr = q.restr, sigma = rcov, validObsEq = validObsEq, 
                      useMatrix = control$useMatrix, solvetol = control$solvetol)
  resids <- yVecAll - xMatAll %*% coef
}

omega <- systemfit:::.calcXtOmegaInv(xMat = xMatAll, sigma = rcov, validObsEq = validObsEq, useMatrix = T, solvetol = control$solvetol)

################################################################################
################################################################################
##### Log-likelihood function?
ll <- function(dat){
  S <- var(dat)
  n <- nrow(dat)
  p <- ncol(dat)
  one <- (-(n*p)/2) * log(2*pi)
  two <- (n/2) * log(abs(S))
  three <- (n/2) * sum(diag(solve(S) %*% S))
  logLik <- one - two - three
  logLik
}

################################################################################
################################################################################
##### Maximum likelihood confidence intervals for partial correlations
mle_CI <- function(X, alpha){
  X <- as.matrix(X)
  # X: data frame
  if(!require("qgraph")) install.packages("qgraph")
  if(!require("Matrix")) install.packages("Matrix")
  # number of observations (rows)
  n <- nrow(X)
  # number of variables (columns)
  p <- ncol(X)
  ## compute maximum likelihood estimator
  ## for covariance matrix
  mle_cov <- crossprod(scale(X, scale = F))/n
  ## compute maximum likelihood estimator of precision matrix
  ## (inverse covariance matrix)
  mle_inv <- solve(mle_cov)
  ## standardize and revese sign = partial correaltions
  par_cors <- as.matrix(qgraph::wi2net(mle_inv))
  mle_parcors <- mle_ci_helper(alpha = alpha, par_cors = par_cors, n = n, s = p - 1)
  mle_inv <- mle_parcors$sig_mat * mle_inv
  list(mle_parcors = mle_parcors, mle_inv = mle_inv)
}

mle_ci_helper <- function(alpha, par_cors, s, n){
  # n: sample size
  # s: p - 1 (controlled for)
  # alpha: confidence level
  # par_cors: partial correlations
  mat <- matrix(0,nrow = s + 1, ncol = s + 1)
  CI_ls <- list()
  par_cor <- par_cors[upper.tri(par_cors)]
  cov <- list()
  for(i in 1:length(par_cor)){
    # critical value
    z_crit <- qnorm(1 - alpha/2)
    # standard error
    se <- sqrt(1/((n - s - 3)))
    # z transformation
    z <- log((1 + par_cor[i])/(1 - par_cor[i]))/2
    # z lower bound
    Z_L <- z - z_crit * se
    # Z upper bound
    Z_U <- z + z_crit * se
    rho_L <- (exp(2*Z_L) - 1)/(exp(2*Z_L) + 1)
    rho_U <- (exp(2*Z_U) - 1)/(exp(2*Z_U) + 1)
    CI <- c(rho_L, rho_U)
    CI_ls[[i]] <- CI
    cov[[i]] <- ifelse(CI[1] < 0 & CI[2] > 0, 0, 1)
  }
  ci_dat <- do.call(rbind.data.frame , CI_ls)
  colnames(ci_dat) <- c("low", "up")
  ci_dat$pcor <- unlist(par_cor)
  diag(mat) <- 1
  mat[upper.tri(mat)] <- unlist(cov)
  mat <- as.matrix(Matrix::forceSymmetric(mat))
  list(sig_mat = mat, par_cors = par_cors,
       par_sig = mat * par_cors,
       cis = ci_dat, cov_prob = unlist(cov))
}


# 95 % CI
est_mle_95 <- mle_CI(X, alpha = 1 - 0.95)
# sparsified partial correlation matrix
est_mle_95$mle_parcors$par_sig
# 99 % CI
est_mle_99 <- mle_CI(X, alpha = 1 - 0.99)
# sparsified partial correlation matrix
est_mle_99$mle_parcors$par_sig

################################################################################
################################################################################
##### GraphicalVAR
trevAR <- function(data, nLambda = 50, verbose = TRUE, gamma = 0.5, scale = TRUE, 
                   lambda_beta, lambda_kappa, maxit.in = 100, maxit.out = 100, 
                   deleteMissings = TRUE, penalize.diagonal = TRUE, lambda_min_kappa = 0.05, 
                   lambda_min_beta = lambda_min_kappa, mimic = c("current", "0.1.2", "0.1.4", "0.1.5", "0.2"), 
                   vars, beepvar, dayvar, idvar, lags = 1, centerWithin = TRUE, 
                   likelihood = c("unpenalized", "penalized"), manual = FALSE){
  if(OLS == TRUE){lambda_beta <- lambda_kappa <- 0}
  if(manual == TRUE){
    nLambda = 50
    verbose = TRUE
    gamma = .5
    scale = T
    lambda_min_beta <- lambda_min_kappa <- .05
    mimic = "current"
    vars = colnames(data)
    centerWithin = T
    likelihood = "unpenalized"
    manual = T
    penalize.diagonal = T
    maxit.in <- maxit.out <- 100
    deleteMissings = T
    lags = 1
    i = 1
    if(!is.null(Sacha)){if(Sacha == TRUE){data <- Data; vars <- Vars}}
  } else {
    mimic <- match.arg(mimic)
  }
  if(mimic == "0.1.2"){
    if(lambda_min_beta != lambda_min_kappa){warning("mimic = 0.1.2 only uses lambda_min_kappa, not lambda_min_beta")}
    if(lambda_min_kappa != 0.01){warning("Set lambda_min_kappa = 0.01 to mimic 0.1.2 default behavior")}
  }
  if(is.list(data) && !is.data.frame(data)){
    if(!("data_c" %in% names(data) & "data_l" %in% names(data))){stop("'data_c' and 'data_l' must be contained in 'data'")}
    data_c <- data$data_c
    data_l <- data$data_l
  } else {
    if(mimic == "0.1.5"){
      if(is.data.frame(data)){data <- as.matrix(data)}
      stopifnot(is.matrix(data))
      data <- scale(data, TRUE, scale)
      data_c <- data[-1, , drop = FALSE]
      data_l <- cbind(1, data[-nrow(data), , drop = FALSE])
    } else {
      if(manual == FALSE){
        data <- graphicalVAR:::tsData(as.data.frame(data), vars = vars, 
                                      beepvar = beepvar, dayvar = dayvar, idvar = idvar, 
                                      scale = scale, centerWithin = centerWithin, lags = lags)
      } else {
        data <- graphicalVAR:::tsData(as.data.frame(data), vars = vars, 
                                      scale = scale, centerWithin = centerWithin, lags = lags)
      }
      data_c <- data$data_c
      data_l <- data$data_l
    }
  }
  data_c <- as.matrix(data_c)
  data_l <- as.matrix(data_l)
  Nvar <- ncol(data_c)
  Ntime <- nrow(data_c)
  if(any(is.na(data_c)) || any(is.na(data_l))){
    if(deleteMissings){
      warnings("Data with missings deleted")
      missing <- rowSums(is.na(data_c)) > 0 | rowSums(is.na(data_l)) > 0
      data_c <- data_c[!missing, ]
      data_l <- data_l[!missing, ]
    } else {
      stop("Missing data not supported")
    }
  }
  if(manual == TRUE){
    if(!exists("lambda_beta") | !exists("lambda_kappa")){
      if(mimic == "0.1.2"){
        lams <- graphicalVAR:::SparseTSCGM_lambdas(data_l, data_c, nLambda, lambda.min.ratio = lambda_min_kappa)
      } else {
        lams <- graphicalVAR:::generate_lambdas(data_l, data_c, nLambda, 
                                                nLambda, lambda_min_kappa = lambda_min_kappa, 
                                                lambda_min_beta = lambda_min_beta, penalize.diagonal = penalize.diagonal, 
                                                version0.1.4 = mimic == "0.1.4")
      }
      if(!exists("lambda_beta")){lambda_beta <- lams$lambda_beta}
      if(!exists("lambda_kappa")){lambda_kappa <- lams$lambda_kappa}
    }
  } else {
    if(missing(lambda_beta) | missing(lambda_kappa)){
      if(mimic == "0.1.2"){
        lams <- SparseTSCGM_lambdas(data_l, data_c, nLambda, lambda.min.ratio = lambda_min_kappa)
      } else {
        lams <- generate_lambdas(data_l, data_c, nLambda, 
                                 nLambda, lambda_min_kappa = lambda_min_kappa, 
                                 lambda_min_beta = lambda_min_beta, penalize.diagonal = penalize.diagonal, 
                                 version0.1.4 = mimic == "0.1.4")
      }
      if(missing(lambda_beta)){lambda_beta <- lams$lambda_beta}
      if(missing(lambda_kappa)){lambda_kappa <- lams$lambda_kappa}
    }
  }
  Nlambda_beta <- length(lambda_beta)
  Nlambda_kappa <- length(lambda_kappa)
  lambdas <- expand.grid(kappa = lambda_kappa, beta = lambda_beta)
  Estimates <- vector("list", nrow(lambdas))
  if(verbose){pb <- txtProgressBar(0, nrow(lambdas), style = 3)}
  for(i in seq_len(nrow(lambdas))){
    if(lambdas$beta[i] == 0 & lambdas$kappa[i] == 0){
      X <- data_l
      Y <- data_c
      nY <- ncol(Y)
      nX <- ncol(X)
      n <- nrow(X)
      beta <- t(Y) %*% X %*% solve(t(X) %*% X)
      S <- 1/(nrow(Y) - 1) * (t(Y) %*% Y - t(Y) %*% X %*% t(beta) - beta %*% t(X) %*% Y + beta %*% t(X) %*% X %*% t(beta))
      S <- (S + t(S))/2
      if(any(eigen(S)$value < 0)){stop("Residual covariances not postive definite")}
      kappa <- solve(S)
      kappa <- (kappa + t(kappa))/2
      lik1 = determinant(kappa)$modulus[1]
      lik2 <- sum(diag(kappa %*% S))
      pdO = sum(sum(kappa[upper.tri(kappa, diag = FALSE)] != 0))
      pdB = sum(sum(beta != 0))
      LLk <- (n/2) * (lik1 - lik2)
      LLk0 <- (n/2) * (-lik2)
      EBIC <- -2 * LLk + (log(n)) * (pdO + pdB) + (pdO + pdB) * 4 * gamma * log(2 * nY)
      Estimates[[i]] <- list(beta = beta, kappa = kappa, EBIC = EBIC)
    } else {
      tryres <- try(Rothmana(data_l, data_c, lambdas$beta[i], 
                             lambdas$kappa[i], gamma = gamma, maxit.in = maxit.in, 
                             maxit.out = maxit.out, penalize.diagonal = penalize.diagonal, 
                             mimic = mimic, likelihood = likelihood))
      if(is(tryres, "try-error")){
        Estimates[[i]] <- list(beta = matrix(NA, Nvar, Nvar + 1), 
                               kappa = matrix(NA, Nvar, Nvar), 
                               EBIC = Inf, error = tryres)
      } else {
        Estimates[[i]] <- tryres
      }
    }
    if(verbose){setTxtProgressBar(pb, i)}
  }
  if(verbose){close(pb)}
  lambdas$ebic <- sapply(Estimates, "[[", "EBIC")
  if(all(lambdas$ebic == Inf)){stop("No model estimated without error")}
  min <- which.min(lambdas$ebic)
  Results <- Estimates[[min]]
  if(length(lambda_beta) > 1){
    if(lambdas$beta[[min]] == min(lambda_beta)){message("Minimal tuning parameter for beta selected.")}
  }
  if(length(lambda_kappa) > 1){
    if(lambdas$kappa[[min]] == min(lambda_kappa)){message("Minimal tuning parameter for kappa selected.")}
  }
  colnames(Results$beta) <- colnames(data_l)
  rownames(Results$beta) <- colnames(data_c)
  Results$PCC <- graphicalVAR:::computePCC(Results$kappa)
  if(1 %in% lags){
    Results$PDC <- graphicalVAR:::computePDC(Results$beta[, c("1", paste0(data$vars, "_lag1"))], Results$kappa)
    if(length(lags) > 1){warning("Partial directed correlations only computed for lag 1 network.")}
  }
  Results$path <- lambdas
  Results$labels <- colnames(data_c)
  if(is.null(Results$labels)){Results$labels <- paste0("V", seq_len(ncol(data_c)))}
  rownames(Results$beta) <- colnames(Results$kappa) <- rownames(Results$kappa) <- colnames(Results$PCC) <- rownames(Results$PCC) <- colnames(Results$PDC) <- rownames(Results$PDC) <- Results$labels
  Results$gamma <- gamma
  Results$allResults <- Estimates
  Results$N <- nrow(data_c)
  Results$data <- data
  class(Results) <- "graphicalVAR"
  return(Results)
}

################################################################################
##### Another function for the graphicalVAR
Rothmana <- function(X, Y, lambda_beta, lambda_kappa, convergence = 0.0001, 
                     gamma = 0.5, maxit.in = 100, maxit.out = 100, penalize.diagonal, 
                     interceptColumn = 1, mimic = "current", likelihood = c("unpenalized", "penalized")){
  likelihood <- match.arg(likelihood)
  nY <- ncol(Y)
  nX <- ncol(X)
  if(missing(penalize.diagonal)){
    if(mimic == "0.1.2"){
      penalize.diagonal <- nY != nX
    } else {
      penalize.diagonal <- (nY != nX - 1) & (nY != nX)
    }
  }
  lambda_mat <- matrix(lambda_beta, nX, nY)
  if(!penalize.diagonal){
    if(nY == nX){
      add <- 0
    } else if(nY == nX - 1){
      add <- 1
    } else {
      stop("Beta is not P x P or P x P+1, cannot detect diagonal.")
    }
    for(i in 1:min(c(nY, nX))){lambda_mat[i + add, i] <- 0}
  }
  if(!is.null(interceptColumn) && !is.na(interceptColumn)){lambda_mat[interceptColumn, ] <- 0}
  n <- nrow(X)
  beta_ridge <- graphicalVAR:::beta_ridge_C(X, Y, lambda_beta)
  beta <- matrix(0, nX, nY)
  it <- 0
  repeat{
    it <- it + 1
    kappa <- graphicalVAR:::Kappa(beta, X, Y, lambda_kappa)
    beta_old <- beta
    beta <- graphicalVAR:::Beta_C(kappa, beta, X, Y, lambda_beta, lambda_mat, convergence, maxit.in)
    if(sum(abs(beta - beta_old)) < (convergence * sum(abs(beta_ridge)))){break}
    if(it > maxit.out){warning("Model did NOT converge in outer loop"); break}
  }
  ZeroIndex <- which(kappa == 0, arr.ind = TRUE)
  WS <- (t(Y) %*% Y - t(Y) %*% X %*% beta - t(beta) %*% t(X) %*% Y + t(beta) %*% t(X) %*% X %*% beta)/(nrow(X))
  if(any(eigen(WS, only.values = TRUE)$values < -sqrt(.Machine$double.eps))){stop("Residual covariance matrix is not non-negative definite")}
  if(likelihood == "unpenalized"){
    if(nrow(ZeroIndex) == 0){
      out4 <- suppressWarnings(glasso(WS, rho = 0, trace = FALSE))
    } else {
      out4 <- suppressWarnings(glasso(WS, rho = 0, zero = ZeroIndex, trace = FALSE))
    }
    lik1 <- determinant(out4$wi)$modulus[1]
    lik2 <- sum(diag(out4$wi %*% WS))
  } else {
    lik1 <- determinant(kappa)$modulus[1]
    lik2 <- sum(diag(kappa %*% WS))
  }
  pdO = sum(sum(kappa[upper.tri(kappa, diag = FALSE)] != 0))
  if(mimic == "0.1.2"){
    pdB = sum(sum(beta != 0))
  } else {
    pdB = sum(sum(beta[lambda_mat != 0] != 0))
  }
  LLk <- (n/2) * (lik1 - lik2)
  LLk0 <- (n/2) * (-lik2)
  EBIC <- -2 * LLk + (log(n)) * (pdO + pdB) + (pdO + pdB) * 4 * gamma * log(2 * nY)
  return(list(beta = t(beta), kappa = kappa, EBIC = EBIC))
}

#####
trevMLGVAR <- function(data, vars, beepvar, dayvar, idvar, scale = TRUE, 
                       centerWithin = TRUE, gamma = 0.5, verbose = TRUE, 
                       subjectNetworks = TRUE, lambda_min_kappa_fixed = 0.001, 
                       lambda_min_beta_fixed = 0.001, lambda_min_kappa = 0.05, 
                       lambda_min_beta = lambda_min_kappa, lambda_min_glasso = 0.01, ...){
  if(missing(idvar)){stop("'idvar' must be assigned")}
  dataPrepped <- tsData(data, vars = vars, beepvar = beepvar, 
                        dayvar = dayvar, idvar = idvar, scale = scale, 
                        centerWithin = centerWithin)
  if(verbose){message("Estimating fixed networks")}
  ResFixed <- graphicalVAR(dataPrepped, lambda_min_kappa = lambda_min_kappa_fixed, 
                           lambda_min_beta = lambda_min_beta_fixed, 
                           gamma = gamma, ...)
  if(verbose){message("Estimating between-subjects network")}
  meansData <- dataPrepped$data_means
  meansData <- meansData[, names(meansData) != idvar]
  meansData <- meansData[rowMeans(is.na(meansData)) != 1, ]
  ResBetween <- qgraph::EBICglasso(cov(meansData), nrow(meansData), 
                                   gamma, returnAllResults = TRUE, 
                                   lambda.min.ratio = lambda_min_glasso)
  IDs <- unique(dataPrepped$data[[idvar]])
  idResults <- list()
  if(!identical(subjectNetworks, FALSE)){
    if(isTRUE(subjectNetworks)){subjectNetworks <- IDs}
    if(verbose){
      message("Estimating subject-specific networks")
      pb <- txtProgressBar(max = length(subjectNetworks), style = 3)
    }
    for(i in seq_along(subjectNetworks)){
      capture.output({
        idResults[[i]] <- try(suppressWarnings(graphicalVAR(dataPrepped$data[dataPrepped$data[[idvar]] == subjectNetworks[i], ], 
                                                            vars = dataPrepped$vars, beepvar = dataPrepped$beepvar, 
                                                            dayvar = dataPrepped$dayvar, idvar = dataPrepped$idvar, 
                                                            scale = scale, lambda_min_kappa = lambda_min_kappa, 
                                                            lambda_min_beta = lambda_min_beta, gamma = gamma, 
                                                            centerWithin = centerWithin, ..., verbose = FALSE)))
      })
      if(verbose){setTxtProgressBar(pb, i)}
      if(is(idResults[[i]], "try-error")){idResults[[i]] <- list()}
    }
    if(verbose){close(pb)}
  } else {
    idResults <- lapply(seq_along(IDs), function(x) list())
  }
  Results <- list(fixedPCC = ResFixed$PCC, fixedPDC = ResFixed$PDC, 
                  fixedResults = ResFixed, betweenNet = ResBetween$optnet, 
                  betweenResults = ResBetween, ids = IDs, 
                  subjectPCC = lapply(idResults, "[[", "PCC"), 
                  subjectPDC = lapply(idResults, "[[", "PDC"), 
                  subjecResults = idResults)
  class(Results) <- "mlGraphicalVAR"
  return(Results)
}

#####
trevTSdata <- function(data, vars, beepvar, dayvar, idvar, lags = 1, 
                       scale = TRUE, centerWithin = TRUE, deleteMissings = TRUE){
  . <- NULL
  data <- as.data.frame(data)
  if(missing(idvar)){
    idvar <- "ID"
    data[[idvar]] <- 1
  }
  if(missing(dayvar)){
    dayvar <- "DAY"
    data[[dayvar]] <- 1
  }
  if(missing(beepvar)){
    beepvar <- "BEEP"
    data <- data %>% dplyr::group_by_(dayvar, idvar) %>% 
      dplyr::mutate_(BEEP = ~seq_len(n()))
  }
  if(missing(vars)){
    vars <- names(data[!names(data) %in% c(idvar, dayvar, beepvar)])
  }
  data <- data[, c(vars, idvar, dayvar, beepvar)]
  data[, vars] <- scale(data[, vars], TRUE, scale)
  MeansData <- data %>% dplyr::group_by_(idvar) %>% 
    dplyr::summarise_at(funs(mean(., na.rm = TRUE)), .vars = vars)
  if(centerWithin){
    if(length(unique(data[[idvar]])) > 1){
      data <- data %>% dplyr::group_by_(idvar) %>% 
        dplyr::mutate_at(funs(scale(., center = TRUE, scale = FALSE)), .vars = vars)
    }
  }
  augData <- data
  beepsPerDay <- eval(substitute(dplyr::summarize_(data %>% group_by_(idvar, dayvar), 
                                                   first = ~min(beepvar, na.rm = TRUE), 
                                                   last = ~max(beepvar, na.rm = TRUE)), 
                                 list(beepvar = as.name(beepvar))))
  allBeeps <- expand.grid(unique(data[[idvar]]), unique(data[[dayvar]]), 
                          seq(min(data[[beepvar]], na.rm = TRUE), 
                              max(data[[beepvar]], na.rm = TRUE)))
  names(allBeeps) <- c(idvar, dayvar, beepvar)
  allBeeps <- eval(substitute({
    allBeeps %>% dplyr::left_join(beepsPerDay, by = c(idvar, dayvar)) %>% 
      dplyr::group_by_(idvar, dayvar) %>% 
      dplyr::filter_(~BEEP >= first, ~BEEP <= last) %>% 
      dplyr::arrange_(idvar, dayvar, beepvar)
  }, list(BEEP = as.name(beepvar))))
  augData <- augData %>% dplyr::right_join(allBeeps, by = c(idvar, dayvar, beepvar))
  data_c <- augData %>% ungroup %>% dplyr::select_(.dots = vars)
  data_l <- do.call(cbind, lapply(lags, function(l){
    data_lagged <- augData %>% dplyr::group_by_(idvar, dayvar) %>% 
      dplyr::mutate_at(funs(shift), .vars = vars) %>% ungroup %>% 
      dplyr::select_(.dots = vars)
    names(data_lagged) <- paste0(vars, "_lag", l)
    data_lagged
  }))
  if(deleteMissings){
    isNA <- rowSums(is.na(data_c)) > 0 | rowSums(is.na(data_l)) > 0
    data_c <- data_c[!isNA, ]
    data_l <- data_l[!isNA, ]
  }
  Results <- list(data = augData, data_c = data_c[, vars], 
                  data_l = cbind(1, data_l), data_means = MeansData, vars = vars, 
                  idvar = idvar, dayvar = dayvar, beepvar = beepvar, lags = lags)
  class(Results) <- "tsData"
  return(Results)
}

shift <- function(x, n){c(rep(NA, n), x)[1:length(x)]}

################################################################################
################################################################################
##### To fit your own models with SEM...
# No periods in colnames---change these to underscores
nVar <- ncol(dat$Y)

# Lambda (identity):
Lambda <- diag(nVar * 2)

# Theta (zeros):
Theta <- matrix(0, nVar * 2, nVar * 2)

# Beta (block): This uses the estimated beta matrix without intercepts
gvarBeta <- ifelse(beta[,-1] == 0, 0, NA)
O <- matrix(0, nVar, nVar)
Beta <- rbind(cbind(O, O), cbind(gvarBeta, O))

# Latent network (exo block for t-1, cont network for t):
gvarKappa <- ifelse(kappa == 0, 0, NA)
Omega_psi <- rbind(cbind(matrix(NA, nVar, nVar), O), cbind(O, gvarKappa))
diag(Omega_psi) <- 0

# Free latent scale:
delta_psi <- diag(NA, nVar * 2)

# Fit model:
lvfit <- lvnet(dat$full, lambda = Lambda, theta = Theta, 
               omega_psi = Omega_psi, beta = Beta, delta_psi = delta_psi, 
               fitFunction = "ML", exogenous = 1:nVar)

# Compare to latent variable model:
# Lambda (equality constrained factor loadings):
Lambda <- rbind(cbind(c(1, paste0("l", 2:nVar)), 0), 
                cbind(0, c(1, paste0("l", 2:nVar))))

# Equality constrained residuals:
Theta <- diag(nVar * 2)
for(i in 1:nVar){Theta[i, nVar + i] <- Theta[nVar + i, i] <- NA}
diag(Theta) <- c(paste0("t", 1:nVar), paste0("t", 1:nVar))

# Beta (regression t-1 -> t1):
Beta <- matrix(c(0, 0, NA, 0), 2, 2, byrow = TRUE)

# Psi (exo var and residual)
Psi <- diag(NA, 2)

# Equality constrained residuals:
Theta[1:nVar,1:nVar] <- NA_character_

# Fit model:
lvfit2 <- lvnet(augData, lambda = Lambda, theta = Theta, beta = Beta, 
                psi = Psi, fitFunction = "ML", exogenous = 1:nVar)

# Compare the models:
tab <- lvnetCompare(lvfit, lvfit2)
