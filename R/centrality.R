#' Makes table of centrality values
#'
#' Lots of ways to do centrality
#'
#' @param Wmats Model(s)
#' @param scale logical
#' @param which.net character
#' @param labels character
#' @param relative logical
#' @param weighted logical
#' @param signed logical
#'
#' @return Table
#' @export
#'
#' @examples
#' 1 + 1
centTable <- function(Wmats, scale = TRUE, which.net = "temporal", labels = NULL,
                      relative = FALSE, weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(Wmats, "mlGVAR"))){
    Wmats <- switch(which.net, between = Wmats$betweenNet, Wmats$fixedNets)}
  if("SURnet" %in% c(names(Wmats), names(attributes(Wmats)))){
    if("SURnet" %in% names(Wmats)){Wmats <- Wmats$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p", 'i')[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc", "interactions"))
    Wmats <- Wmats[[ifelse(which.net == "contemporaneous", "contemporaneous", ifelse(which.net == 'interactions', 'interactions', 'temporal'))]]
    if(which.net == "pdc"){Wmats <- Wmats$PDC}
    if(which.net == 'interactions'){names(Wmats)[1] <- 'adjMat'}
  } else if(startsWith(tolower(which.net), 'i')){stop('Interaction centrality only supported for temporal networks.')}
  if("adjMat" %in% names(Wmats)){Wmats <- t(Wmats$adjMat)}
  if(any(grepl("lag", dimnames(Wmats)))){dimnames(Wmats) <- lapply(dimnames(Wmats), function(z) gsub("[.]lag1.*|[.]y$", "", z))}
  if(!is.list(Wmats)){Wmats <- list(Wmats)}
  if(any(sapply(Wmats, ncol) == 1)){stop("Not supported for single-node graphs")}
  #names(Wmats) <- qgraph:::fixnames(Wmats, "graph ")
  names(Wmats) <- fnames(Wmats, 'graph ')
  centOut <- lapply(Wmats, centAuto, which.net = which.net, weighted = weighted, signed = signed)
  for(g in seq_along(centOut)){
    if(!is(centOut[[g]], "centrality_auto")){
      #names(centOut[[g]]) <- qgraph:::fixnames(centOut[[g]], "type ")
      names(centOut[[g]]) <- fnames(centOut[[g]], 'type ')
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
        if(scale){
          #centOut[[i]][["node.centrality"]][, j] <- qgraph:::scale2(centOut[[i]][["node.centrality"]][, j])
          centOut[[i]][["node.centrality"]][, j] <- scaleNA(centOut[[i]][["node.centrality"]][, j])
        }
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

#' Mimics centrality_auto
#'
#' Yep yep
#'
#' @param x Not sure
#' @param which.net character
#' @param weighted logical
#' @param signed logical
#'
#' @return A list of stuff
#' @export
#'
#' @examples
#' 1 + 1
centAuto <- function(x, which.net = "temporal", weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(x, "mlGVAR"))){
    x <- switch(which.net, between = x$betweenNet, x$fixedNets)}
  if("SURnet" %in% c(names(x), names(attributes(x)))){
    if("SURnet" %in% names(x)){x <- x$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p", "i")[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc", "interactions"))
    x <- x[[ifelse(which.net == "contemporaneous", "contemporaneous", ifelse(which.net == 'interactions', 'interactions', "temporal"))]]
    if(which.net == "pdc"){x <- x$PDC}
    if(which.net == 'interactions'){names(x)[1] <- 'adjMat'}
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

#' Mimics centralityPlot from qgraph
#'
#' Indeed
#'
#' @param Wmats Models
#' @param scale character
#' @param which.net character
#' @param include character
#' @param labels character
#' @param orderBy character
#' @param decreasing logical
#' @param plot logical
#' @param verbose logical
#' @param weighted logical
#' @param signed logical
#'
#' @return Plot
#' @export
#'
#' @examples
#' 1 + 1
centPlot <- function(Wmats, scale = c("z-scores", "raw", "raw0", "relative"),
                     which.net = "temporal", include = "all", labels = NULL,
                     orderBy = NULL, decreasing = FALSE, plot = TRUE,
                     verbose = TRUE, weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(Wmats, "mlGVAR"))){
    Wmats <- switch(which.net, between = Wmats$betweenNet, Wmats$fixedNets)}
  if(is.logical(scale)){scale <- ifelse(scale, "z-scores", "raw")}
  #invisible(suppressMessages(require(ggplot2)))
  measure <- value <- node <- type <- NULL
  scale <- tryCatch({match.arg(scale)}, error = function(e){scale})
  #scale <- match.arg(scale)
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
