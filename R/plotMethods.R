#' @describeIn plotNet For ggms!
#' @export
plot.ggm <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                     predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                     scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                     plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                     binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.SURnet <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                        predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                        scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                        plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                        binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.mlGVAR <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                        predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                        scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                        plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                        binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.lmerVAR <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                         predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                         scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                         plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                         binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.mgmSim <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                        predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                        scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                        plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                        binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.mlVARsim <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                          predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                          scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                          plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                          binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.simMLgvar <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                           predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                           scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                           plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                           binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotPower
#' @export
plot.mnetPower <- function(x, by = 'type', yvar = 'default', yadd = NULL, hline = .8,
                           xlab = 'Number of cases', title = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotPower, args)
}

#' @rdname plotBoot
#' @export
plot.bootNet <- function(x, type = 'edges', net = 'temporal', plot = 'all', cor = .7,
                         order = 'mean', ci = .95, pairwise = TRUE, interactions = TRUE,
                         labels = NULL, title = NULL, cis = 'quantile', true = NULL,
                         errbars = FALSE, vline = FALSE, threshold = FALSE,
                         difference = FALSE, color = FALSE, text = FALSE,
                         textPos = 'value', multi = NULL, directedDiag = FALSE, ...){
  args <- as.list(match.call())[-1]
  do.call(plotBoot, args)
}

#' Plot method for output of resample function
#'
#' Description
#'
#' @param x output from resample function
#' @param what network, bootstrap, or coefs
#' @param ... Additional arguments.
#'
#' @return A plot
#' @export
#'
#' @examples
#' 1 + 1
plot.resample <- function(x, what = 'network', ...){
  args <- tryCatch({list(...)}, error = function(e){list()})
  if(isTRUE(what) | is(what, 'numeric')){
    args$threshold <- ifelse(!'threshold' %in% names(args), what, args$threshold)
    what <- 'network'
  }
  what <- match.arg(tolower(what), c('network', 'bootstrap', 'coefs'))
  FUN <- switch(what, network = plotNet, bootstrap = plotBoot, coefs = plotCoefs)
  if(what == 'bootstrap' & !is.null(x$call$moderators)){
    if('fit0' %in% names(x)){
      if(!any(grepl(':', unlist(x$fit0$call$type)))){
        x$call <- replace(x$call, 'moderators', list(moderators = NULL))
      }
    }
  }
  args$XX <- x
  names(args)[which(names(args) == 'XX')] <- switch(what, network = 'x', bootstrap = 'x', coefs = 'fit')
  do.call(FUN, args)
}
