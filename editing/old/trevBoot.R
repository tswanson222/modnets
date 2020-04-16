bootMod <- function(data, nboots, threshold = FALSE, ci = .95, verbose = TRUE, ...){
  ci <- (1 - ci)/2
  n <- nrow(data)
  inds <- data.frame(replicate(nboots, sample(1:n, n, replace = TRUE)))
  if(verbose){pb <- txtProgressBar(min = 0, max = nboots + 1, style = 3)}
  fit0 <- fitNetwork(data = data, threshold = threshold, ...)
  if(verbose){setTxtProgressBar(pb, 1)}
  fits <- lapply(1:nboots, function(z){
    fit <- fitNetwork(data = data[inds[[z]], ], threshold = threshold, ...)
    if(verbose){setTxtProgressBar(pb, z + 1)}
    return(fit)
  })
  adj0 <- fit0$adjMat[lower.tri(fit0$adjMat)]
  str0 <- colSums(abs(fit0$adjMat))
  ei0 <- colSums(fit0$adjMat)
  adj1 <- do.call(cbind, lapply(lapply(fits, '[[', "adjMat"), function(z) z[lower.tri(z)]))
  colnames(adj1) <- colnames(inds) <- paste0("boot", 1:nboots)
  adj2 <- as.data.frame(t(apply(adj1, 1, function(z){
    ints <- quantile(z, c(ci, .5, 1 - ci))
    return(c(ints, ifelse(ints[1] <= 0 & ints[3] >= 0, FALSE, TRUE)))
  })))
  colnames(adj2) <- c("lower", "median", "upper", "sig")
  adj2$sig <- as.logical(adj2$sig)
  nodes <- names(fits[[1]]$fitobj)
  v <- length(nodes)
  node1 <- rep(nodes[-v], (v - 1):1)
  node2 <- matrix(nodes, v, v)[lower.tri(matrix(nodes, v, v))]
  edges <- paste0(node1, "--", node2)
  e <- length(edges)
  edge1 <- rep(edges[-e], (e - 1):1)
  edge2 <- matrix(edges, e, e)[lower.tri(matrix(edges, e, e))]
  adjmeans <- data.frame(mean = rowMeans(adj1), sd = apply(adj1, 1, sd))
  adj2 <- data.frame(edge = edges, node1 = node1, node2 = node2, adjmeans, adj2,
                     sample = adj0, sample_lower = adj0 + (adjmeans$sd * qnorm(ci)),
                     sample_upper = adj0 + (adjmeans$sd * qnorm(1 - ci)))
  rownames(adj1) <- edges
  strengths <- EIs <- list()
  for(i in 1:length(nodes)){
    strengths[[i]] <- adj1[which(node1 == nodes[i] | node2 == nodes[i]), -c(1:2)]
    EIs[[i]] <- colSums(strengths[[i]])
    strengths[[i]] <- colSums(abs(strengths[[i]]))
  }
  strengths <- do.call(cbind.data.frame, strengths)
  EIs <- do.call(cbind.data.frame, EIs)
  colnames(strengths) <- colnames(EIs) <- nodes
  estDiffs <- function(ests, ci = .025){
    comps <- t(combn(ncol(ests), 2))
    diffOut <- as.data.frame(t(sapply(1:nrow(comps), function(z){
      ints <- quantile(ests[, comps[z, 1]] - ests[, comps[z, 2]], c(ci, 1 - ci))
      return(c(ints, ifelse(ints[1] <= 0 & ints[2] >= 0, FALSE, TRUE)))
    })))
    colnames(diffOut) <- c("lower", "upper", "sig")
    diffOut$sig <- as.logical(diffOut$sig)
    return(diffOut)
  }
  sdiffs <- data.frame(node1 = node1, node2 = node2, estDiffs(strengths, ci))
  eidiffs <- data.frame(node1 = node1, node2 = node2, estDiffs(EIs, ci))
  smeans <- data.frame(node = nodes, mean = colMeans(strengths), sd = apply(strengths, 2, sd))
  smeans <- data.frame(smeans, sample = str0, sample_lower = str0 + (smeans$sd * qnorm(ci)),
                       sample_upper = str0 + (smeans$sd * qnorm(1 - ci)))
  eimeans <- data.frame(node = nodes, mean = colMeans(EIs), sd = apply(EIs, 2, sd))
  eimeans <- data.frame(eimeans, sample = ei0, sample_lower = ei0 + (eimeans$sd * qnorm(ci)),
                        sample_upper = ei0 + (eimeans$sd * qnorm(1 - ci)))
  rownames(smeans) <- rownames(eimeans) <- 1:v
  edge_diffs <- data.frame(edge1 = edge1, edge2 = edge2, estDiffs(t(adj1), ci))
  out <- list(edges = list(edge_cis = adj2, edge_diffs = edge_diffs),
              strength = list(str_means = smeans, str_diffs = sdiffs),
              EI = list(ei_means = eimeans, ei_diffs = eidiffs),
              boots = list(edges = adj1, strengths = strengths, EIs = EIs, inds = inds))
  attr(out, "call") <- match.call()
  attr(out, "ci") <- 1 - ci * 2
  return(out)
}
