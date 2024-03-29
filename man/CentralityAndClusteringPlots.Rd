% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centralityClustering.R
\name{CentralityAndClusteringPlots}
\alias{CentralityAndClusteringPlots}
\alias{centPlot}
\alias{clustPlot}
\alias{plotCentrality}
\title{Plots for node centrality values or clustering coefficients}
\usage{
centPlot(
  Wmats,
  scale = c("z-scores", "raw", "raw0", "relative"),
  which.net = "temporal",
  include = "all",
  labels = NULL,
  orderBy = NULL,
  decreasing = FALSE,
  plot = TRUE,
  verbose = TRUE,
  weighted = TRUE,
  signed = TRUE
)

clustPlot(
  Wmats,
  scale = c("z-scores", "raw", "raw0", "relative"),
  include = "all",
  labels = NULL,
  orderBy = NULL,
  decreasing = FALSE,
  plot = TRUE,
  signed = TRUE,
  verbose = TRUE
)

plotCentrality(
  Wmats,
  which.net = "temporal",
  scale = TRUE,
  labels = NULL,
  plot = TRUE,
  centrality = "all",
  clustering = "Zhang"
)
}
\arguments{
\item{Wmats}{Output from one of the primary \code{modnets} functions.}

\item{scale}{If \code{"z-scores"}, then standardized values will be plotted.
If \code{"relative"}, then values will be scaled relative to the largest
value on each measure. \code{"raw"} can be used to plot raw values.}

\item{which.net}{Only applies to SUR networks, as well as those fit with the
\code{\link{mlGVAR}} function. Character string to indicate which type of
network to compute centrality values for. Options are \code{"temporal"} for
the temporal network, \code{"contemporaneous"} for the contemporaneous
network, \code{"PDC"} for the partial directed correlation network, and
\code{"interactions"} for the temporal interaction network.}

\item{include}{Character vector of which centrality measures to plot.
\code{"Betweenness"} and \code{"Closeness"} are available for all types of
network. \code{"Strength"} and \code{"ExpectedInfluence"} are only
available for GGMs. And \code{"InStrength", "OutStrength",
  "InExpectedInfluence", "OutExpectedInfluence"} are only available for SUR
networks. Defaults to \code{"all"}}

\item{labels}{Character vector listing the node names. If \code{NULL}, then
the names specified by the model are used.}

\item{orderBy}{Character string specifying which measure to order values by.}

\item{decreasing}{Logical. Only relevant if \code{orderBy} is specified.
Determines whether values are organized from highest to lowest, or vice
versa.}

\item{plot}{Logical. Determines whether to plot the output or not.}

\item{verbose}{Logical. Determines whether to return a message about the plot
(messages are only shown if values are scaled).}

\item{weighted}{See \code{\link{centTable}} or \code{\link{clustTable}}.}

\item{signed}{See \code{\link{centTable}} or \code{\link{clustTable}}.}

\item{centrality}{Character vector of centrality measures to plot. Defaults
to \code{"all"}.}

\item{clustering}{Character vector of clustering measures to plot. Defaults
to \code{"Zhang"}.}
}
\value{
A plot of centrality values or clustering coefficients for several
measures.
}
\description{
Mimics the \code{\link[qgraph:centralityPlot]{qgraph::centralityPlot}} and
\code{\link[qgraph:clusteringPlot]{qgraph::clusteringPlot}} functions. The
purpose of revising this function was to make it compatible with outputs from
the \code{modnets} package.
}
\details{
The only utility of the \code{\link{plotCentrality}} function is as an easy
way to combine centrality measures and clustering coefficients into a single
plot.
}
\examples{
x <- fitNetwork(ggmDat)

centPlot(x)
clustPlot(x)
plotCentrality(x)
}
\seealso{
\code{\link{centTable}, \link{clustTable}, \link{centAuto},
  \link{clustAuto}, \link[qgraph:centralityPlot]{qgraph::centralityPlot},
  \link[qgraph:clusteringPlot]{qgraph::clusteringPlot}}
}
