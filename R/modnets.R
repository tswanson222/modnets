#' modnets: Modeling Moderated Networks
#'
#' Methods for modeling and plotting various types of moderated networks,
#' including tools for model selection. Model selection tools can be employed
#' for any type of moderated network, and include methods based on the LASSO as
#' well as resampling techniques such as bootstrapping, multi-sample splitting,
#' and stability selection. The primary model types include: \itemize{
#' \item{Cross-section moderated networks} \item{Temporal (idiographic)
#' moderated networks} \item{Multi-level moderated networks} }
#'
#' @section Details:
#'
#'   \tabular{ll}{ Package: \tab modnets \cr Title: \tab Modeling Moderated
#'   Networks \cr Version: \tab 0.9.0 \cr Author: \tab Trevor Swanson \cr
#'   Maintainer: \tab <trevorswanson222@gmail.com> \cr URL: \tab
#'   \url{https://github.com/tswanson222/modnets} \cr BugReports: \tab
#'   \url{https://github.com/tswanson222/modnets/issues} \cr License: \tab GPL-3
#'   \cr Imports: \tab abind, corpcor, ggplot2, glinternet, glmnet, gridExtra,
#'   gtools, igraph, interactionTest, leaps, lme4, lmerTest, Matrix, methods,
#'   mvtnorm, parallel, pbapply, plyr, psych, qgraph, reshape2, systemfit \cr
#'   Suggests: \tab sn, arm, hierNet, psychTools \cr Date: \tab 2021-09-26 }
#'
#' @author Trevor Swanson
#'
#'   \strong{Maintainer:} Trevor Swanson <trevorswanson222@gmail.com>
#'
#' @references
#'
#' Swanson, T. J. (2020). \emph{Modeling moderators in psychological networks}
#' (Publication No. 28000912) \[Doctoral dissertation, University of Kansas\].
#' ProQuest Dissertations Publishing.
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
#' @importFrom utils capture.output combn data sessionInfo setTxtProgressBar
#'   txtProgressBar stack

utils::globalVariables(
  c('B', 'Y', 'b', 'group', 'id1', 'id2', 'lower', 'node', 'ord', 'pred', 'prop0',
    'sig', 'subN', 'upper', 'value', 'x', 'xmax', 'xmin', 'y', 'ymax', 'ymin',
    'Type', 'measure')
)
