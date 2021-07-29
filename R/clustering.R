#' Mimics clusteringTable from qgraph
#'
#' But with moderators!
#'
#' @param Wmats Models
#' @param scale logical
#' @param labels character
#' @param relative logical
#' @param signed logical
#'
#' @return Table
#' @export
#'
#' @examples
#' 1 + 1
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
#' Mimics clustcoef_auto
#'
#' @param x something
#' @param thresholdWS something
#' @param thresholdON something
#'
#' @return A list of stuff
#' @export
#'
#' @examples
#' 1 + 1
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
#' Mimics clusteringPlot
#'
#' @param Wmats Models
#' @param scale character
#' @param include character
#' @param labels character
#' @param orderBy character
#' @param decreasing logical
#' @param plot logical
#' @param signed logical
#' @param verbose logical
#'
#' @return Plot
#' @export
#'
#' @examples
#' 1 + 1
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
