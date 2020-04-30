### ======================================================================== ###
### ======================================================================== ###
##### centTable: mimics centralityTable
centTable <- function(Wmats, scale = TRUE, which.net = "temporal", labels = NULL, 
                      relative = FALSE, weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(Wmats, "mlGVAR"))){
    Wmats <- switch(which.net, between = Wmats$betweenNet, Wmats$fixedNets)}
  if("SURnet" %in% c(names(Wmats), names(attributes(Wmats)))){
    if("SURnet" %in% names(Wmats)){Wmats <- Wmats$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p")[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc"))
    Wmats <- Wmats[[ifelse(which.net == "contemporaneous", "contemporaneous", "temporal")]]
    if(which.net == "pdc"){Wmats <- Wmats$PDC}
  }
  if("adjMat" %in% names(Wmats)){Wmats <- t(Wmats$adjMat)}
  if(any(grepl("lag", dimnames(Wmats)))){dimnames(Wmats) <- lapply(dimnames(Wmats), function(z) gsub("[.]lag1.*|[.]y$", "", z))}
  if(!is.list(Wmats)){Wmats <- list(Wmats)}
  if(any(sapply(Wmats, ncol) == 1)){stop("Not supported for single-node graphs")}
  names(Wmats) <- qgraph:::fixnames(Wmats, "graph ")
  centOut <- lapply(Wmats, centAuto, which.net = which.net, weighted = weighted, signed = signed)
  for(g in seq_along(centOut)){
    if(!is(centOut[[g]], "centrality_auto")){
      names(centOut[[g]]) <- qgraph:::fixnames(centOut[[g]], "type ")
      for(t in seq_along(centOut[[g]])){
        if(!is.null(labels)){
          centOut[[g]][[t]][["node.centrality"]][["node"]] <- labels
        } else if(!is.null(rownames(centOut[[g]][[t]][["node.centrality"]]))){
          centOut[[g]][[t]][["node.centrality"]][["node"]] <- rownames(centOut[[g]][[t]][["node.centrality"]])
        } else {
          centOut[[g]][[t]][["node.centrality"]][["node"]] <- paste("Node", seq_len(nrow(centOut[[g]][[t]][["node.centrality"]])))
        }
        centOut[[g]][[t]]$node.centrality$graph <- names(centOut)[g]
        centOut[[g]][[t]]$node.centrality$type <- names(centOut[[g]])[t]
      }
    } else {
      centOut[[g]]$node.centrality$graph <- names(centOut)[g]
      if(!is.null(labels)){
        centOut[[g]][["node.centrality"]][["node"]] <- labels
      } else if(!is.null(rownames(centOut[[g]][["node.centrality"]]))){
        centOut[[g]][["node.centrality"]][["node"]] <- rownames(centOut[[g]][["node.centrality"]])
      } else {
        centOut[[g]][["node.centrality"]][["node"]] <- paste("Node", seq_len(nrow(centOut[[g]][["node.centrality"]])))
      }
    }
  }
  isList <- sapply(centOut, function(x) !"centrality_auto" %in% class(x))
  if(any(isList)){
    for(l in which(isList)){centOut <- c(centOut, centOut[[l]])}
    centOut <- centOut[-which(isList)]
  }
  for(i in seq_along(centOut)){
    if(relative | scale){
      if(relative & scale){warning("Using 'relative' and 'scale' together is not recommended")}
      for(j in which(sapply(centOut[[i]][["node.centrality"]], mode) == "numeric")){
        if(scale){centOut[[i]][["node.centrality"]][, j] <- qgraph:::scale2(centOut[[i]][["node.centrality"]][, j])}
        if(relative){
          mx <- max(abs(centOut[[i]][["node.centrality"]][, j]), na.rm = TRUE)
          if(mx != 0){centOut[[i]][["node.centrality"]][, j] <- centOut[[i]][["node.centrality"]][, j]/mx}
        }
        attributes(centOut[[i]][["node.centrality"]][, j]) <- NULL
      }
    }
  }
  wideCent <- plyr::rbind.fill(lapply(centOut, "[[", "node.centrality"))
  if(is.null(wideCent$type)){wideCent$type <- NA}
  longCent <- reshape2::melt(wideCent, variable.name = "measure", id.var = c("graph", "type", "node"))
  if(any(is.nan(longCent$value))){warning("NaN detected in centrality measures. Try relative = FALSE")}
  return(longCent)
}

##### centAuto: mimics centrality_auto
centAuto <- function(x, which.net = "temporal", weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(x, "mlGVAR"))){
    x <- switch(which.net, between = x$betweenNet, x$fixedNets)}
  if("SURnet" %in% c(names(x), names(attributes(x)))){
    if("SURnet" %in% names(x)){x <- x$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p")[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc"))
    x <- x[[ifelse(which.net == "contemporaneous", "contemporaneous", "temporal")]]
    if(which.net == "pdc"){x <- x$PDC}
  }
  if("adjMat" %in% names(x)){x <- t(x$adjMat)}
  if(any(grepl("lag", dimnames(x)))){dimnames(x) <- lapply(dimnames(x), function(z) gsub("[.]lag1.*|[.]y$", "", z))}
  if(is.list(x)){return(lapply(x, centAuto, which.net = which.net, weighted = weighted, signed = signed))}
  if(!weighted){x <- sign(x)}
  if(!signed){x <- abs(x)}
  if(!is.matrix(x)){stop("The input network must be an adjacency or weights matrix")}
  diag(x) <- 0
  directed.gr <- ifelse(isSymmetric.matrix(object = x, tol = 0.000000000001), FALSE, TRUE)
  weighted.gr <- ifelse(all(qgraph::mat2vec(x) %in% c(0, 1)), FALSE, TRUE)
  net_qg <- qgraph::qgraph(x, diag = FALSE, labels = colnames(x), DoNotPlot = TRUE, minimum = 0)
  centr <- qgraph::centrality(net_qg)
  if(directed.gr & !weighted.gr){
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness, Closeness = centr$Closeness, 
      InDegree = centr$InDegree, OutDegree = centr$OutDegree, 
      OutExpectedInfluence = centr$OutExpectedInfluence, 
      InExpectedInfluence = centr$InExpectedInfluence))
  }
  if(directed.gr & weighted.gr){
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness, Closeness = centr$Closeness, 
      InStrength = centr$InDegree, OutStrength = centr$OutDegree, 
      OutExpectedInfluence = centr$OutExpectedInfluence, 
      InExpectedInfluence = centr$InExpectedInfluence))
  }
  if(!directed.gr & !weighted.gr){
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness/2, Closeness = centr$Closeness, 
      Degree = centr$OutDegree, ExpectedInfluence = centr$OutExpectedInfluence))
  }
  if(!directed.gr & weighted.gr){ 
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness/2, Closeness = centr$Closeness, 
      Strength = centr$OutDegree, ExpectedInfluence = centr$OutExpectedInfluence))
  }
  row.names(centr1) <- colnames(x)
  log <- capture.output({
    graph <- igraph::graph.adjacency(
      adjmatrix = 1 * (x != 0), 
      mode = ifelse(directed.gr, "directed", "undirected"))
    comps <- igraph::components(graph)
    largcomp <- comps$membership == which.max(comps$csize)
  })
  if(sum(largcomp) < ncol(x) & sum(largcomp) > 1){
    x2 <- x[largcomp, largcomp]
    clos <- qgraph::centrality(qgraph::qgraph(
      x2, diag = FALSE, labels = colnames(x)[largcomp], 
      DoNotPlot = TRUE, minimum = 0))$Closeness
    centr1$Closeness[largcomp] <- clos
    centr1$Closeness[!largcomp] <- NA
  }
  net_ig_abs <- igraph::graph.adjacency(
    adjmatrix = abs(1/x), mode = ifelse(directed.gr, "directed", "undirected"), 
    weighted = ifelse(weighted.gr, list(TRUE), list(NULL))[[1]], diag = FALSE)
  edgebet <- igraph::edge.betweenness(graph = net_ig_abs, directed = directed.gr)
  el <- data.frame(igraph::get.edgelist(graph = net_ig_abs), stringsAsFactors = FALSE)
  edgebet <- merge(el, edgebet, by = 0)
  edgebet$Row.names <- NULL
  names(edgebet) <- c("from", "to", "edgebetweenness")
  edgebet <- edgebet[order(edgebet$edgebetweenness, decreasing = TRUE), ]
  ShortestPathLengths <- centr$ShortestPathLengths
  rownames(ShortestPathLengths) <- colnames(ShortestPathLengths) <- colnames(x)
  Res <- list(node.centrality = centr1, edge.betweenness.centrality = edgebet, 
              ShortestPathLengths = ShortestPathLengths)
  class(Res) <- c("list", "centrality_auto")
  return(Res)
}

##### centPlot: mimics centralityPlot
centPlot <- function(Wmats, scale = c("z-scores", "raw", "raw0", "relative"), 
                     which.net = "temporal", include = "all", labels = NULL, 
                     orderBy = NULL, decreasing = FALSE, plot = TRUE, 
                     verbose = TRUE, weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(Wmats, "mlGVAR"))){
    Wmats <- switch(which.net, between = Wmats$betweenNet, Wmats$fixedNets)}
  if(is.logical(scale)){scale <- ifelse(scale, "z-scores", "raw")}
  invisible(suppressMessages(require(ggplot2)))
  measure <- value <- node <- type <- NULL
  scale <- match.arg(scale)
  include0 <- c("Degree", "Strength", "OutDegree", "InDegree", "OutStrength", 
                "InStrength", "Closeness", "Betweenness", "ExpectedInfluence", 
                "OutExpectedInfluence", "InExpectedInfluence")
  if(all(tolower(include) == "all")){
    include <- include0
  } else if(isTRUE(attr(Wmats, "SURnet")) & which.net != "contemporaneous"){
    include0 <- include0[!grepl("Degree|^S|^E", include0)]
    include <- include0[grep(paste(tolower(include), collapse = "|"), tolower(include0))]
  } else if(which.net %in% c("between", "contemporaneous")){
    include0 <- include0[!grepl("Degree|^Out|^In", include0)]
    include <- include0[grep(paste(tolower(include), collapse = "|"), tolower(include0))]
  }
  include <- match.arg(include, c("Degree", "Strength", "OutDegree", "InDegree", "OutStrength", 
                                  "InStrength", "Closeness", "Betweenness", "ExpectedInfluence", 
                                  "OutExpectedInfluence", "InExpectedInfluence"), several.ok = TRUE)
  if(scale == "z-scores" & verbose & plot){message("Note: z-scores are shown on x-axis.")}
  if(scale == "relative" & verbose & plot){message("Note: relative centrality indices are shown on x-axis.")}
  Long <- centTable(Wmats = Wmats, scale = (scale == "z-scores"), labels = labels, 
                    which.net = which.net, relative = (scale == "relative"), 
                    weighted = weighted, signed = signed)
  Long <- subset(Long, measure %in% include)
  Long$measure <- factor(Long$measure)
  if(ifelse(is.null(orderBy), FALSE, ifelse(orderBy == "default", TRUE, FALSE))){
    nodeLevels <- unique(gtools::mixedsort(
      as.character(Long$node), decreasing = decreasing))
  } else if(!is.null(orderBy)){
    nodeLevels <- names(sort(tapply(
      Long$value[Long$measure == orderBy], 
      Long$node[Long$measure == orderBy], mean), 
      decreasing = decreasing))
  } else {
    nodeLevels <- rev(unique(as.character(Long$node)))
  }
  Long$node <- factor(as.character(Long$node), levels = nodeLevels)
  Long <- Long[gtools::mixedorder(Long$node), ]
  if(length(unique(Long$type)) > 1){
    g <- ggplot(Long, aes(x = value, y = node, group = type, colour = type))
  } else {
    g <- ggplot(Long, aes(x = value, y = node, group = type))
  }
  g <- g + geom_path() + xlab("") + ylab("") + geom_point() + theme_bw()
  if(length(unique(Long$graph)) > 1){
    g <- g + facet_grid(graph ~ measure, scales = "free")
  } else {
    g <- g + facet_grid(~measure, scales = "free")
  }
  if(scale == "raw0"){g <- g + xlim(0, NA)}
  if(plot){plot(g)} else {invisible(g)}
}

 
### ======================================================================== ###
### ======================================================================== ###
##### clustTable: mimics clusteringTable
clustTable <- function(Wmats, scale = TRUE, labels = NULL, 
                       relative = FALSE, signed = TRUE){
  if("SURnet" %in% c(names(Wmats), names(attributes(Wmats)))){
    if("SURnet" %in% names(Wmats)){Wmats <- Wmats$SURnet}
    Wmats <- Wmats$contemporaneous$adjMat
  } else if("adjMat" %in% names(Wmats)){
    Wmats <- Wmats$adjMat
  }
  if(any(grepl("lag", dimnames(Wmats)))){
    dimnames(Wmats) <- lapply(dimnames(Wmats), function(z) gsub("[.]lag1.*|[.]y$", "", z))
  }
  if(!is.list(Wmats)){Wmats <- list(Wmats)}
  if(any(sapply(Wmats, ncol) == 1)){stop("Not supported for single-node graphs")}
  syms <- sapply(Wmats, isSymmetric)
  if(any(!syms)){
    if(all(!syms)){stop("No symmetrical graphs detected")}
    warning(paste0(sum(!syms), " Nonsymmetrical graph", ifelse(sum(!syms) > 1, "s ", " "), "removed"))
    Wmats <- Wmats[-which(!syms)]
  }
  names(Wmats) <- qgraph:::fixnames(Wmats, "graph ")
  clustOut <- lapply(Wmats, clustAuto)
  for(g in seq_along(clustOut)){
    if(!is(clustOut[[g]], "clustcoef_auto")){
      names(clustOut[[g]]) <- qgraph:::fixnames(clustOut[[g]], "type ")
      for(t in seq_along(clustOut[[g]])){
        if(!is.null(labels)){
          clustOut[[g]][[t]][["node"]] <- labels
        } else if(!is.null(rownames(clustOut[[g]][[t]]))){
          clustOut[[g]][[t]][["node"]] <- rownames(clustOut[[g]][[t]])
        } else {
          clustOut[[g]][[t]][["node"]] <- paste("Node", seq_len(nrow(clustOut[[g]][[t]])))
        }
        clustOut[[g]][[t]]$graph <- names(clustOut)[g]
        clustOut[[g]][[t]]$type <- names(clustOut[[g]])[t]
      }
    } else {
      clustOut[[g]]$graph <- names(clustOut)[g]
      if(!is.null(labels)){
        clustOut[[g]][["node"]] <- labels
      } else if(!is.null(rownames(clustOut[[g]]))){
        clustOut[[g]][["node"]] <- rownames(clustOut[[g]])
      } else {
        clustOut[[g]][["node"]] <- paste("Node", seq_len(nrow(clustOut[[g]])))
      }
    }
  }
  isList <- sapply(clustOut, function(x) !"clustcoef_auto" %in% class(x))
  if(any(isList)){
    for(l in which(isList)){clustOut <- c(clustOut, clustOut[[l]])}
    clustOut <- clustOut[-which(isList)]
  }
  for(i in seq_along(clustOut)){
    if(any(grepl("signed_", names(clustOut[[i]])))){
      clustOut[[i]] <- clustOut[[i]][, sapply(clustOut[[i]], mode) != "numeric" | grepl("signed_", names(clustOut[[i]])) == signed]
      names(clustOut[[i]]) <- gsub("signed_", "", names(clustOut[[i]]))
    }
    names(clustOut[[i]]) <- gsub("clust", "", names(clustOut[[i]]))
    if(relative | scale){
      if(relative & scale){warning("Using 'relative' and 'scale' together is not recommended")}
      for(j in which(sapply(clustOut[[i]], mode) == "numeric")){
        if(scale){clustOut[[i]][, j] <- qgraph:::scale2(clustOut[[i]][, j])}
        if(relative){
          mx <- max(abs(clustOut[[i]][, j]), na.rm = TRUE)
          if(mx != 0){clustOut[[i]][, j] <- clustOut[[i]][, j]/mx}
        }
        attributes(clustOut[[i]][, j]) <- NULL
      }
    }
  }
  WideCent <- plyr::rbind.fill(clustOut)
  if(is.null(WideCent$type)){WideCent$type <- NA}
  LongCent <- reshape2::melt(WideCent, variable.name = "measure", id.var = c("graph", "type", "node"))
  return(LongCent)
}

##### clustAuto: mimics clustcoef_auto
clustAuto <- function(x, thresholdWS = 0, thresholdON = 0){
  if("SURnet" %in% c(names(x), names(attributes(x)))){
    if("SURnet" %in% names(x)){x <- x$SURnet}
    x <- x$contemporaneous$adjMat
  } else if("adjMat" %in% names(x)){
    x <- x$adjMat
  }
  if(any(grepl("lag", dimnames(x)))){
    dimnames(x) <- lapply(dimnames(x), function(z) gsub("[.]lag1.*|[.]y$", "", z))
  }
  if(is.list(x)){
    return(lapply(x, clustAuto, thresholdWS = thresholdWS, thresholdON = thresholdWS))
  }
  dim = dim(x)
  if(is.null(dim) || length(dim) != 2){stop("adjacency is not two-dimensional")}
  if(!is.numeric(x)){stop("adjacency is not numeric")}
  if(dim[1] != dim[2]){stop("adjacency is not square")}
  if(max(abs(x - t(x)), na.rm = TRUE) > 0.000000000001){stop("adjacency is not symmetric")}
  if(min(x, na.rm = TRUE) < -1 || max(x, na.rm = TRUE) > 1){x <- x/max(abs(x))}
  weighted.gr <- ifelse(all(abs(x) %in% c(0, 1)), FALSE, TRUE)
  signed.gr <- ifelse(all(x >= 0), FALSE, TRUE)
  net_ig <- igraph::graph.adjacency(
    adjmatrix = abs(x), mode = "undirected", 
    weighted = ifelse(weighted.gr, list(TRUE), list(NULL))[[1]], diag = FALSE)
  cb <- igraph::transitivity(net_ig, type = "barrat", isolates = "zero")
  cw <- qgraph:::clustWS(x, thresholdWS)
  cz <- qgraph:::clustZhang(x)
  co <- qgraph:::clustOnnela(x, thresholdON)
  if(!signed.gr & !weighted.gr){output <- cbind(clustWS = cw[, 1])}
  if(!signed.gr & weighted.gr){
    output <- cbind(clustWS = cw[, 1], clustZhang = cz[, 1], 
                    clustOnnela = co[, 1], clustBarrat = cb)
  }
  if(signed.gr & !weighted.gr){
    output <- cbind(clustWS = cw[, 1], signed_clustWS = cw[, 2])
  }
  if(signed.gr & weighted.gr){
    output <- cbind(clustWS = cw[, 1], signed_clustWS = cw[, 2], 
                    clustZhang = cz[, 1], signed_clustZhang = cz[, 2], 
                    clustOnnela = co[, 1], signed_clustOnnela = co[, 2], 
                    clustBarrat = cb)
  }
  output[is.na(output)] <- 0
  Res <- data.frame(output)
  class(Res) <- c("data.frame", "clustcoef_auto")
  rownames(Res) <- colnames(x)
  return(Res)
}

##### clustPlot: mimics clusteringPlot
clustPlot <- function(Wmats, scale = c("z-scores", "raw", "raw0", "relative"), 
                      include = "all", labels = NULL, orderBy = NULL, 
                      decreasing = FALSE, plot = TRUE, signed = TRUE, 
                      verbose = TRUE){
  if(is.logical(scale)){scale <- ifelse(scale, "z-scores", "raw")}
  invisible(suppressMessages(require(ggplot2)))
  measure <- value <- node <- type <- NULL
  scale <- match.arg(scale)
  if(scale == "z-scores" & verbose & plot){message("Note: z-scores are shown on x-axis.")}
  if(scale == "relative" & verbose & plot){message("Note: relative centrality indices are shown on x-axis.")}
  Long <- clustTable(Wmats = Wmats, scale = (scale == "z-scores"), labels = labels, 
                     relative = (scale == "relative"), signed = signed)
  Long$value[!is.finite(Long$value)] <- 0
  if(all(include == "all")){include <- c("WS", "Zhang", "Onnela", "Barrat")}
  include <- match.arg(include, c("WS", "Zhang", "Onnela", "Barrat"), several.ok = TRUE)
  Long <- subset(Long, measure %in% include)
  Long$measure <- factor(Long$measure)
  if(ifelse(is.null(orderBy), FALSE, ifelse(orderBy == "default", TRUE, FALSE))){
    nodeLevels <- unique(gtools::mixedsort(
      as.character(Long$node), decreasing = decreasing))
  } else if(!is.null(orderBy)){
    nodeLevels <- names(sort(tapply(
      Long$value[Long$measure == orderBy], 
      Long$node[Long$measure == orderBy], mean), 
      decreasing = decreasing))
  } else {
    nodeLevels <- rev(unique(as.character(Long$node)))
  }
  Long$node <- factor(as.character(Long$node), levels = nodeLevels)
  Long <- Long[gtools::mixedorder(Long$node), ]
  if(length(unique(Long$type)) > 1){
    g <- ggplot(Long, aes(x = value, y = node, group = type, colour = type))
  } else {
    g <- ggplot(Long, aes(x = value, y = node, group = type))
  }
  g <- g + geom_path() + xlab("") + ylab("") + geom_point() + theme_bw()
  if(length(unique(Long$graph)) > 1){
    g <- g + facet_grid(graph ~ measure, scales = "free")
  } else {
    g <- g + facet_grid(~measure, scales = "free")
  }
  if(scale == "raw0"){g <- g + xlim(0, NA)}
  if(plot){plot(g)} else {invisible(g)}
}


### ======================================================================== ###
### ======================================================================== ###
##### plotCentrality: plot centrality (and clustering) for multiple networks
plotCentrality <- function(Wmats, which.net = "temporal", scale = TRUE, 
                           labels = NULL, plot = TRUE, centrality = "all", 
                           clustering = "Zhang"){
  if(any(c("ggm", "SURnet", "mlGVAR") %in% names(attributes(Wmats)))){Wmats <- list(net1 = Wmats)}
  if(all(sapply(Wmats, function(z) isTRUE(attr(z, "mlGVAR"))))){
    Wmats <- lapply(Wmats, function(z) switch(
      which.net, between = z$betweenNet, z$fixedNets))
  }
  if(any(grepl("ggm", lapply(Wmats, function(z) names(attributes(z)))))){which.net <- "contemporaneous"}
  if(length(unique(lapply(Wmats, checkInclude))) != 1){stop("All networks must be of the same type")}
  if(is.null(names(Wmats))){names(Wmats) <- paste0("net", seq_along(Wmats))}
  which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc"))
  c0 <- c01 <- do.call(rbind, lapply(seq_along(Wmats), function(z){
    cbind.data.frame(
      centTable(Wmats[[z]], scale = scale, which.net = which.net, 
                labels = labels), group = names(Wmats)[z])
  }))
  if(all(centrality != "all")){
    include0 <- checkInclude(Wmats[[1]], which.net = which.net)
    include0 <- include0[!grepl(ifelse(
      which.net != "contemporaneous", "Degree|^S|^E", 
      "Degree|^Out|^In"), include0)]
    centrality <- include0[grep(paste(tolower(
      centrality), collapse = "|"), tolower(include0))]
    c0 <- c01 <- subset(c0, measure %in% centrality)
  }
  if(which.net == "contemporaneous" & clustering != FALSE){
    c1 <- do.call(rbind, lapply(seq_along(Wmats), function(z){
      z1 <- clustTable(Wmats[[z]], scale = scale, labels = labels)
      z1 <- z1[z1$measure == ifelse(is.logical(clustering), "Zhang", clustering), ]
      z1$measure <- "Clust. coef."
      z1$node <- as.character(z1$node)
      rownames(z1) <- 1:nrow(z1)
      return(cbind.data.frame(z1, group = names(Wmats)[z]))
    }))
    c01 <- rbind(c0, c1)
  }
  c01 <- c01[order(c01$node), ]
  c01 <- c01[order(c01$group), ]
  rownames(c01) <- 1:nrow(c01)
  c01$node <- stringr::str_sub(c01$node, 1, 6)
  if(!plot){
    list2env(list(c0 = c0, c1 = c1, c01 = c01), .GlobalEnv)
  } else {
    invisible(suppressMessages(require(ggplot2)))
    g1 <- ggplot(c01, aes(x = value, y = node, group = group, color = group, shape = group)) +
      geom_path(alpha = 1, size = 1) + geom_point(size = 2) + xlab("") + ylab("") + theme_bw() + 
      facet_grid(. ~ measure, scales = "free") + scale_x_continuous(breaks = c(-1, 0, 1)) + 
      theme(axis.line.x = element_line(colour = "black"),
            axis.ticks.x = element_line(colour = "black"),
            axis.ticks.y = element_line(colour = "white", size = 0),
            axis.text.y = element_text(colour = "black"),
            axis.text.x = element_text(angle = 45, colour = "black"))
    g1
  }
}

##### checkInclude
checkInclude <- function(x, which.net = "temporal"){
  if("SURnet" %in% c(names(x), names(attributes(x)))){
    if("SURnet" %in% names(x)){x <- x$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p")[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc"))
    x <- x[[ifelse(which.net == "contemporaneous", "contemporaneous", "temporal")]]
    if(which.net == "pdc"){x <- x$PDC}
  }
  if("adjMat" %in% names(x)){x <- t(x$adjMat)}
  directed <- !isTRUE(all.equal(x, t(x), check.attributes = FALSE))
  weighted <- any(!x %in% c(0, 1))
  include <- c("Degree", "Strength", "OutDegree", "InDegree", "OutStrength", 
               "InStrength", "Closeness", "Betweenness", "ExpectedInfluence", 
               "OutExpectedInfluence", "InExpectedInfluence")
  if(directed){
    include <- c("OutStrength", "InStrength", "Closeness", "Betweenness",
                 "OutExpectedInfluence", "InExpectedInfluence")
  } else {
    include <- c("Strength", "Closeness", "Betweenness", "ExpectedInfluence")
  }
  if(!weighted){include <- gsub("Strength", "Degree", include)}
  return(include)
}


##### expinf: one- or two-step expected influence
expinf <- function(x, ei = 'both', type = 'pairwise', directed = FALSE, ...){
  #if(directed){stop('Not implemented for directed graphs yet')}
  if(ei %in% c('pairwise', 'interactions')){type <- ei; ei <- 'both'}
  type <- match.arg(tolower(type), c('pairwise', 'interactions'))
  if(!is(x, 'matrix')){
    FUN <- match.fun(switch(type, pairwise = 'net', interactions = 'netInts'))
    if(type == 'interactions'){formals(FUN)$avg <- !directed}
    x <- FUN(x, ...)
  }
  if(!directed){diag(x) <- 0}
  p <- ncol(x)
  ei1 <- rowSums(x)
  ei2 <- ei1 + rowSums(x * matrix(rep(ei1, p), p, p, T))
  out <- setNames(list(ei1, ei2), paste0('step', 1:2))
  if(ifelse(is.numeric(ei), ei %in% 1:2, FALSE)){out <- out[[ei]]}
  return(out)
}
