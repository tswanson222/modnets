# @describeIn plotNet For ggms!
#' @rdname plotNet
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
#' Allows one to plot results from the \code{resample()} function based on a few
#' different options.
#'
#' For the \code{what} argument, the options correspond with calls to the
#' following functions: \itemize{ \item{\code{"network": \link{plotNet}()}}
#' \item{\code{"bootstrap": \link{plotBoot}()}} \item{\code{"coefs":
#' \link{plotCoefs}()}} \item{\code{"stability": \link{plotStability}()}}
#' \item{\code{"pvals": \link{plotPvals}()}} }
#'
#' @param x Output from the \code{resample()} function.
#' @param what Can be one of three options for all \code{resample()} outputs:
#'   \code{what = "network"} will plot the final network model selected from
#'   resampling. \code{what = "bootstrap"} will run \code{bootNet()} based on
#'   the final model to create bootstrapped estimates of confidence bands around
#'   each edge estimate. \code{what = "coefs"} will plot the confidence
#'   intervals based on the model parameters in the final network. Additionally,
#'   if the object was fit with \code{sampMethod = "stability"}, a stability
#'   plot can be created with the \code{"stability"} option. Otherwise, if
#'   \code{sampMethod = "bootstrap"} or \code{sampMethod = "split"}, a plot of
#'   the empirical distribution function of p-values can be displayed with the
#'   \code{"pvals"} option.
#' @param ... Additional arguments.
#'
#' @return Plots various aspects of output from the \code{\link{resample}()}
#'   function.
#' @export
#'
#' @examples
#' \dontrun{
#' x <- resample(data)
#'
#' plot(x, 'network')
#' plot(x, 'bootstrap')
#' plot(x, 'coefs')
#' }
plot.resample <- function(x, what = 'network', ...){
  args <- tryCatch({list(...)}, error = function(e){list()})
  if(isTRUE(what) | is(what, 'numeric')){
    args$threshold <- ifelse(!'threshold' %in% names(args), what, args$threshold)
    what <- 'network'
  }
  what <- match.arg(tolower(what), c('network', 'bootstrap', 'coefs', 'stability', 'pvals'))
  if(what == 'stability'){stopifnot('stability' %in% names(x))}
  if(what == 'pvals'){stopifnot(x$call$sampMethod != 'stability')}
  FUN <- switch(what, network = plotNet, bootstrap = plotBoot, coefs = plotCoefs,
                stability = plotStability, pvals = plotPvals)
  if(what == 'bootstrap' & !is.null(x$call$moderators)){
    if('fit0' %in% names(x)){
      if(!any(grepl(':', unlist(x$fit0$call$type)))){
        x$call <- replace(x$call, 'moderators', list(moderators = NULL))
      }
    }
  }
  args$XX <- x
  names(args)[which(names(args) == 'XX')] <- ifelse(what == 'coefs', 'fit', 'x')
  #names(args)[which(names(args) == 'XX')] <- switch(what, network = 'x', bootstrap = 'x', coefs = 'fit')
  do.call(FUN, args)
}
