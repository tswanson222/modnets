% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mnetPowerSim.R
\name{mnetPowerSim}
\alias{mnetPowerSim}
\title{Power simulator for cross-sectional and idiographic networks}
\usage{
mnetPowerSim(
  niter = 10,
  N = 100,
  p = 5,
  m = FALSE,
  m1 = 0,
  m2 = 0.1,
  sparsity = 0.5,
  lags = NULL,
  trueNet = NULL,
  threshold = TRUE,
  rule = "OR",
  avg = TRUE,
  maxiter = 100,
  saveFits = TRUE,
  saveData = FALSE,
  intercepts = NULL,
  mbinary = FALSE,
  select = NULL,
  vargs = list(),
  type = "g",
  gibbs = TRUE,
  ordinal = FALSE,
  mord = FALSE,
  nLevels = 5,
  minOrd = 3,
  div = 1000,
  modType = "none",
  m1_range = NULL,
  m2_range = c(0.1, 0.3),
  time = TRUE,
  skewErr = FALSE,
  nCores = 1,
  cluster = "mclapply",
  fixedPar = NULL,
  V2 = 1,
  ...
)
}
\arguments{
\item{niter}{Number of iterations/samples to take for each combination of
parameters.}

\item{N}{Numeric value, or vector of sample sizes to generate data with.}

\item{p}{Numeric value, or vector of network sizes.}

\item{m}{If a value is provided then a moderated network will be simulated.
See \code{\link{simNet}} for details.}

\item{m1}{Functions similarly to \code{m2}, except that this argument refers
to the number/probability of main effects of the moderator on any given
node.}

\item{m2}{Numeric. If \code{m2 >= 1}, then this will determine the number of
interaction effects between the moderator and some node in the network. If
a value between 0 and 1 is provided, then this determines the probability
of any given edge being moderated by the moderator.}

\item{sparsity}{Numeric value between 0 and 1. Determines the sparsity of
sampled network matrices.}

\item{lags}{Determines whether the network should be a temporal network or
not. If simulating a temporal network, set to \code{TRUE} or 1.}

\item{trueNet}{The adjacency matrix of the data-generating network model, or
a list containing the adjacency matrix as the first element, and the
interaction matrix as the second element.}

\item{threshold}{See corresponding argument in \code{\link{fitNetwork}}.
Automatically set to \code{TRUE} if \code{select} is not \code{NULL}.}

\item{rule}{Only applies to GGMs (including between-subjects networks) when a
threshold is supplied. The \code{"AND"} rule will only preserve edges when
both corresponding coefficients have p-values below the threshold, while
the \code{"OR"} rule will preserve an edge so long as one of the two
coefficients have a p-value below the supplied threshold.}

\item{avg}{See corresponding argument of \code{\link{netInts}}}

\item{maxiter}{If a model fails to be fit, this determines the maximum number
of iterations to re-try it before giving up. Will also simulate new
datasets at each iteration.}

\item{saveFits}{Logical. Determines whether to save the models fit to each
dataset at each iteration.}

\item{saveData}{Logical. Determines whether to save the datasets generated at
each iteration.}

\item{intercepts}{A vector of means for sampling node values.}

\item{mbinary}{Logical. Determines whether the moderator should be a binary
variable.}

\item{select}{Identifies a variable selection function -- either
\code{\link{varSelect}} or \code{\link{resample}} -- to use for introducing
variable selection at each iteration. The usefulness of this is to mimic a
real-world situation, wherein the researcher may be interested in seeing
how well datasets of different sizes afford models that approximate a true
model after employing iterated variable selection. If \code{TRUE} then this
defaults to \code{"varSelect"}. Highly recommended to use the \code{vargs}
argument to supply necessary information about the parameters of the
variable selection process, such as \code{sampMethod}, \code{criterion},
etc.}

\item{vargs}{A named list of arguments relevant to the variable selection
procedure specified by the \code{select} argument.}

\item{type}{Can supply a variable selection object, such as the output from
either \code{\link{varSelect}} or \code{\link{modSelect}}, can be supplied
to choose a specific constrained model to fit on all iterations. This is
essentially an alternative to \code{select}, in that \code{select} performs
variable selection at each iteration, whereas this argument defines a
constrained model that is applied at every iteration.}

\item{gibbs}{If \code{TRUE}, then Gibbs sampling will be used. Otherwise,
data are generated from the \code{\link[mvtnorm:rmvnorm]{mvtnorm::rmvnorm}}
function based on the partial correlation matrix that is created.}

\item{ordinal}{Logical. Determines whether to generate ordinal values or not.}

\item{mord}{Logical. Determines whether the moderator variable should be
simulated as ordinal.}

\item{nLevels}{Number of levels for the ordinal variables. Only relevant if
\code{ordinal} is not \code{FALSE}.}

\item{minOrd}{The minimum number of unique values allowed for each variable.}

\item{div}{A value to use as a sign that the sampler diverged. Can be
increased based on expected range of values. If a datapoint is larger than
\code{div}, then the sampler will stop.}

\item{modType}{Determines the type of moderation to employ, such as
\code{"none", "full", "partial"}. See \code{\link{simNet}} for details.}

\item{m1_range}{Numeric vector of length 2. The range of values for moderator
main effect coefficients.}

\item{m2_range}{Numeric vector of length 2. The range of values for moderator
interaction effect coefficients.}

\item{time}{If \code{TRUE} then the time it takes to simulate the data is
printed to screen at the end of the sampling.}

\item{skewErr}{The skewness parameter for the \code{alpha} argument in the
\code{\link[sn:rmsn]{sn::rmsn}} function. Only relevant when \code{gibbs =
  FALSE} and no moderator is specified.}

\item{nCores}{Numeric value indicating the number of CPU cores to use for the
resampling. If \code{TRUE}, then the
\code{\link[parallel:detectCores]{parallel::detectCores}} function will be
used to maximize the number of cores available.}

\item{cluster}{Character vector indicating which type of parallelization to
use, if \code{nCores > 1}. Options include \code{"mclapply"} and
\code{"SOCK"}.}

\item{fixedPar}{Numeric. If provided, then this will be set as the
coefficient value for all edges in the network. Provides a way to
standardize the parameter values while varying the sparsity of the network.
If \code{length(fixedPar) == 1}, then the same value will be used for all
parameters. If \code{length(fixedPar) == 2}, then the first value will be
for pairwise relationships, and the second value will be for interaction
terms.}

\item{V2}{If \code{V2 = 1} and \code{m2} is between 0 and 1, the number of
interaction terms in the model will be determined by multiplying \code{m2}
with the number of elements in the interaction matrix and taking the
\code{ceiling}.}

\item{...}{Additional arguments.}
}
\value{
Power simulation results
}
\description{
Samples data based on several parameters, mainly used to see how different
sample sizes perform given various parameterizations when simulating from
network models, especially moderated networks. See \code{\link{simNet}} for
more details about arguments as well as the warning about simulations that
fail.
}
\details{
Evaluates how closely an estimated network is with the true network with
regards to metrics such as sensitivity, specificity, and precision, among
others. Doesn't calculate values for power, but can be used to serve a
similar function as a traditional power analysis based on simulated datasets.
}
\examples{
\donttest{
x <- mnetPowerSim(niter = 10, N = c(100, 200))
summary(x)
plot(x)
}
}
\seealso{
\code{\link{summary.mnetPower}, \link{plotPower}, \link{simNet},
  \link{mlGVARsim}}
}
