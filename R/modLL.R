#' Omnibus log-likelihood of a network and LRT comparing two networks.
#'
#' Computes log-likelihood, AIC, and BIC for a whole network.
#' Also compares two networks using a likelihoo ratio test.
#'
#' @param net0 GGM
#' @param net1 GGM or null
#' @param nodes logical
#' @param lrt logical
#' @param all logical
#' @param d numeric
#' @param alpha numeric
#' @param orderBy character
#' @param decreasing logical
#'
#' @return A named vector
#' @export
#'
#' @examples
#' 1 + 1
modLL <- function(net0, net1 = NULL, nodes = FALSE, lrt = NULL, all = FALSE,
                  d = 4, alpha = .05, orderBy = NULL, decreasing = TRUE){
  # Log-likelihood for whole network
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

  # Log-likelihood and information criteria for single regression
  uni_mll <- function(object){
    getInd <- function(ind, x){
      ind <- c("deviance", "LL_model", "df.residual", "AIC", "BIC")[ind]
      X <- ifelse(ind == "df.residual", list(x$fitobj), list(x$mods))[[1]]
      return(unname(sapply(X, '[[', ind)))
    }
    out <- do.call(cbind.data.frame, lapply(1:5, getInd, x = object))
    colnames(out) <- c("RSS", "LL", "df", "AIC", "BIC")
    rownames(out) <- names(object$mods)
    out
  }

  # Combines mll with information criteria for full network
  omni_mll <- function(object){
    k <- length(object$fitobj)
    n <- nrow(object$data) * k
    ll <- mll(object = object)
    aic <- (2 * ll[2]) - (2 * ll[1])
    bic <- (ll[2] * log(n)) - (2 * ll[1])
    out <- c(ll, AIC = unname(aic), BIC = unname(bic))
    out
  }

  # Likelihood ratio test for comparing two networks
  mod_lrt <- function(object, d = 4, alpha = .05, N = NULL){
    if(is.list(object)){
      if(length(object) > 2){object <- object[1:2]}
      nn <- names(object)
      ll0 <- object[[1]]$LL; df0 <- object[[1]]$df
      ll1 <- object[[2]]$LL; df1 <- object[[2]]$df
      omnibus <- FALSE
    } else {
      if(nrow(object) > 2){object <- object[1:2, ]}
      nn <- rownames(object)
      ll0 <- object[1, 1]; df0 <- object[1, 2]
      ll1 <- object[2, 1]; df1 <- object[2, 2]
      omnibus <- TRUE
    }
    lldiff <- abs(ll0 - ll1) * 2
    dfdiff <- abs(df0 - df1)
    ps <- pchisq(q = lldiff, df = dfdiff, lower.tail = FALSE)
    decision <- c()
    for(i in seq_along(ps)){
      if(ps[i] <= alpha){
        decision[i] <- ifelse(ll0[i] > ll1[i], nn[1], nn[2])
      } else if(ps[i] == 1){
        decision[i] <- "- "
      } else if(ps[i] > alpha){
        if(omnibus){
          decision[i] <- ifelse(df0 < df1, nn[1], nn[2])
        } else {
          decision[i] <- ifelse(df0[i] > df1[i], nn[1], nn[2])
        }
      }
    }
    if(!omnibus){
      if(!is.null(d)){ps <- round(ps, d)}
      out <- data.frame(LL_diff2 = lldiff, Df_diff = dfdiff, pval = ps, decision = decision)
      rownames(out) <- rownames(object[[1]])
    } else {
      RMSEA <- function(X2, df, N){
        rmsea <- sqrt(max(c(((X2/N)/df) - (1/N), 0)))
        lower.l <- function(lambda){(pchisq(q = X2, df = df, ncp = lambda) - .95)}
        lambda.l <- tryCatch({uniroot(f = lower.l, lower = 0, upper = X2)$root}, error = function(e){0})
        rmsea.lower <- sqrt(lambda.l/(N * df))
        upper.l <- function(lambda){(pchisq(q = X2, df = df, ncp = lambda) - .05)}
        lambda.u <- tryCatch({uniroot(f = upper.l, lower = 0, upper = max(N, X2 * 4))$root}, error = function(e){1})
        rmsea.upper <- sqrt(lambda.u/(N * df))
        rmsea.pvalue <- 1 - pchisq(q = X2, df = df, ncp = (N * df * (.05^2)))
        return(c(lower = rmsea.lower, RMSEA = rmsea, upper = rmsea.upper, p.value = rmsea.pvalue))
      }
      rmsea <- RMSEA(lldiff, dfdiff, N)
      if(!is.null(d)){
        ps <- round(ps, d)
        lldiff <- round(lldiff, d)
        rmsea <- round(rmsea, d)
      }
      if(object[2, 2] < object[1, 2]){object <- object[order(object[, 2]), ]}
      out0 <- data.frame(LL_diff2 = c("", lldiff), Df_diff = c("", dfdiff),
                         pval = c("", ps), decision = c("", decision))
      out <- data.frame(object[, 1:2], out0, object[, 3:4])
      attr(out, "RMSEA") <- rmsea
    }
    return(out)
  }

  # Helps simplify interface; you can put two networks in a list for net0
  nn <- paste0("net", 0:1)
  if(length(net0) == 2 & ifelse(
    is.null(net1), TRUE, ifelse(is.logical(net1), net1, FALSE))){
    if(is.logical(net1)){nodes <- net1}
    nn <- ifelse(!is.null(names(net0)), list(names(net0)), list(nn))[[1]]
    net1 <- net0[[2]]; net0 <- net0[[1]]
  }

  # Take the between-subjects network if mlGVAR
  if(isTRUE(attr(net0, "mlGVAR"))){net0 <- net0$betweenNet}

  # If we have a ggm...
  if("ggm" %in% names(attributes(net0))){
    omni0 <- switch(2 - isTRUE('modLL' %in% names(net0)), net0$modLL$omnibus, omni_mll(net0))
    uni0 <- switch(2 - isTRUE('modLL' %in% names(net0)), net0$modLL$nodes, uni_mll(net0))
    if(is.logical(net1)){nodes <- net1; net1 <- NULL}
    if(!is.null(net1)){
      if(isTRUE(attr(net1, "mlGVAR"))){net1 <- net1$betweenNet}
      if(is.null(lrt)){lrt <- TRUE}
      omni1 <- switch(2 - isTRUE('modLL' %in% names(net1)), net1$modLL$omnibus, omni_mll(net1))
      uni1 <- switch(2 - isTRUE('modLL' %in% names(net1)), net1$modLL$nodes, uni_mll(net1))
      omni0 <- rbind(omni0, omni1)
      uni0 <- list(uni0, uni1)
      rownames(omni0) <- names(uni0) <- nn
    }
  } else {
    if(all(sapply(net0, function(z) isTRUE(attr(z, "mlGVAR"))))){
      net0 <- lapply(net0, '[[', "betweenNet")
    }
    stopifnot(all(sapply(net0, function(z) isTRUE(attr(z, "ggm")))))
    yLL <- all(sapply(net0, function(z) 'modLL' %in% names(z)))
    if(is.null(lrt)){lrt <- FALSE}
    if(is.logical(net1)){nodes <- net1}
    omni0 <- t(switch(2 - yLL, sapply(lapply(net0, '[[', 'modLL'), '[[', 'omnibus'), sapply(net0, omni_mll)))
    if(!is.null(orderBy)){
      orderBy <- switch(match.arg(tolower(as.character(
        orderBy)), c("true", "ll", "loglik", "aic", "bic", "df", "rank")),
        ll =, loglik = "LL", aic = "AIC", bic = "BIC", "df")
      net0 <- net0[order(omni0[, orderBy], decreasing = decreasing)]
      #net0 <- net0[order(sapply(net0, function(z){
      #  omni_mll(z)[orderBy]}), decreasing = decreasing)]
    }
    nn <- ifelse(!is.null(names(net0)), list(names(net0)),
                 list(paste0("net", 0:(length(net0) - 1))))[[1]]
    uni0 <- switch(2 - yLL, lapply(lapply(net0, '[[', 'modLL'), '[[', 'nodes'), lapply(net0, uni_mll))
    #omni0 <- do.call(rbind, lapply(net0, omni_mll))
    #uni0 <- lapply(net0, uni_mll)
    rownames(omni0) <- names(uni0) <- nn
  }

  # Collect output
  out <- ifelse(all, list(list(nodes = uni0, omnibus = omni0)),
                ifelse(nodes, list(uni0), list(omni0)))[[1]]

  # Do you want to conduct an LRT?
  if(ifelse(!is.null(net1), lrt, FALSE)){
    lrt0 <- mod_lrt(object = omni0, d = d, alpha = alpha, N = prod(dim(net1$dat$Y)))
    lrt1 <- mod_lrt(object = uni0, d = d, alpha = alpha)
    out <- list(nodes = uni0, LRT = lrt1, omnibus = lrt0)
    if(!all){out <- ifelse(nodes, list(lrt1), list(lrt0))[[1]]}
  }

  # RETURN
  return(out)
}

#' Moderated network table
#'
#' Fits LRT to a list of network models to compare them all against each other.
#'
#' @param fits list of network models
#' @param nodes logical
#' @param orderBy character
#' @param d numeric
#' @param alpha numeric
#' @param decreasing logical
#' @param names logical
#' @param rmsea logical
#'
#' @return A list of tables
#' @export
#'
#' @examples
#' 1 + 1
modTable <- function(fits, nodes = FALSE, orderBy = TRUE, d = 4, alpha = .05,
                     decreasing = TRUE, names = NULL, rmsea = FALSE){
  n <- length(fits)
  stopifnot(is.list(fits) & n > 2)
  if(!is.null(names)){names(fits) <- names}
  if(is.null(names(fits))){names(fits) <- paste0("fit", 1:n)}
  if(all(sapply(fits, function(z) isTRUE(attr(z, "mlGVAR"))))){
    fits <- lapply(fits, '[[', "betweenNet")
  }
  stopifnot(all(sapply(fits, function(z) isTRUE(attr(z, "ggm")))))
  orderBy <- ifelse(is.null(orderBy), FALSE, ifelse(identical(tolower(orderBy), 'lrt'), TRUE, orderBy))
  if(!is.logical(orderBy)){
    orderBy <- switch(match.arg(tolower(as.character(
      orderBy)), c("ll", "loglik", "aic", "bic", "df", "rank")),
      ll =, loglik = "LL", aic = "AIC", bic = "BIC", "df")
    fits <- fits[order(sapply(fits, function(z){
      modLL(z)[orderBy]}), decreasing = decreasing)]
  }
  tt <- combn(n, 2)
  lls0 <- modLL(fits)
  lrts0 <- lapply(seq_len(ncol(tt)), function(i){
    modLL(fits[tt[, i]], d = d, alpha = alpha)})
  out1 <- t(sapply(lrts0, rownames))
  out2 <- t(sapply(lrts0, function(z) as.numeric(z[2, 3:5])))
  out3 <- sapply(lrts0, function(z) z[2, 6])
  out4 <- cbind.data.frame(out1, out2, out3)
  colnames(out4) <- c("net0", "net1", "Chisq", "Df", "pval", "decision")
  rmsea0 <- data.frame(out4[, 1:2], do.call(rbind, lapply(lrts0, attr, "RMSEA")))
  select <- table(out4$decision)
  if(any(names(select) == "- ")){select <- select[names(select) != "- "]}
  lls0 <- cbind(lls0, LRT = numeric(nrow(lls0)))
  lls0[match(names(select), rownames(lls0)), "LRT"] <- unname(select)
  if(isTRUE(orderBy)){
    lls0 <- lls0[order(lls0[, 'LRT'], decreasing = TRUE), ]
  }
  out <- list(LRT = out4, omnibus = lls0, RMSEA = rmsea0)
  if(!rmsea){out$RMSEA <- NULL}
  if(nodes != FALSE){
    lls1 <- modLL(fits, nodes = TRUE)
    stopifnot(length(unique(lapply(lls1, rownames))) == 1)
    nodenames <- unique(lapply(lls1, rownames))[[1]]
    lls2 <- lapply(nodenames, function(fit){
      nodemod <- matrix(unlist(sapply(lls1, function(node){
        node[fit, -1]})), nrow = length(lls1), ncol = 4, byrow = TRUE)
      dimnames(nodemod) <- list(names(fits), c("LL", "df", "AIC", "BIC"))
      return(nodemod)
    })
    names(lls2) <- nodenames
    lrts1 <- lapply(seq_len(ncol(tt)), function(i){
      modLL(fits[tt[, i]], nodes = TRUE, d = d, alpha = alpha)})
    netLRTs <- lapply(colnames(lrts1[[1]]), function(nn){
      n1 <- data.frame(out4[, 1:2], "|", do.call(rbind, lapply(lrts1, '[[', nn)))
      colnames(n1)[-c(1:2)] <- c("|", nodenames)
      return(n1)
    })
    names(netLRTs) <- colnames(lrts1[[1]])
    deci <- netLRTs$decision[, -c(1:3)]
    deci2 <- do.call(cbind.data.frame, lapply(1:ncol(deci), function(z){
      z <- table(deci[, z])
      if(any(names(z) == "- ")){z <- z[names(z) != "- "]}
      zz <- setNames(numeric(length(fits)), names(fits))
      zz[match(names(z), names(zz))] <- unname(z)
      return(zz)
    }))
    colnames(deci2) <- colnames(deci)
    out <- list(nodes = lls2, LRT = netLRTs, counts = deci2)
    if(is.character(nodes)){
      out$decision <- netLRTs$decision
      out$LRT <- NULL
    }
  }
  attr(out, "alpha") <- alpha
  return(out)
}
