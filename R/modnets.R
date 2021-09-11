#' modnets: A package for fitting moderated network models
#'
#' The modnets package is all about fitting and analyzing moderated network
#' models. Currently in development -- more updates and details to come.
#'
#' @docType package
#' @name modnets
#'
#' @importFrom grDevices rgb terrain.colors
#' @importFrom graphics axis hist layout legend lines plot.new text
#' @importFrom methods formalArgs is
#' @import stats
#' @import ggplot2
#' @import qgraph
#' @importFrom utils capture.output combn data sessionInfo setTxtProgressBar txtProgressBar stack

utils::globalVariables(
  c('B', 'Y', 'b', 'group', 'id1', 'id2', 'lower', 'node', 'ord', 'pred', 'prop0',
    'sig', 'subN', 'upper', 'value', 'x', 'xmax', 'xmin', 'y', 'ymax', 'ymin',
    'Type', 'measure')
)
