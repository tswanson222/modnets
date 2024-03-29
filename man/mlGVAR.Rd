% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlGVAR.R
\name{mlGVAR}
\alias{mlGVAR}
\title{Fit GVAR models with multilevel data}
\usage{
mlGVAR(
  data,
  m = NULL,
  selectFUN = NULL,
  subjectNets = FALSE,
  idvar = "ID",
  exogenous = TRUE,
  center = TRUE,
  scale = TRUE,
  fixedType = "g",
  betweenType = "g",
  centerWithin = TRUE,
  scaleWithin = FALSE,
  rule = "OR",
  threshold = "none",
  verbose = TRUE,
  pcor = FALSE,
  fixedArgs = NULL,
  betweenArgs = NULL,
  bm = FALSE,
  beepno = NULL,
  dayno = NULL,
  deleteMissing = TRUE,
  ...
)
}
\arguments{
\item{data}{\code{n x k} dataframe or matrix}

\item{m}{Character vector or numeric vector indicating the moderator(s), if
any. Can also specify \code{"all"} to make every variable serve as a
moderator, or \code{0} to indicate that there are no moderators. If the
length of \code{m} is \code{k - 1} or longer, then it will not be possible
to have the moderators as exogenous variables. Thus, \code{exogenous} will
automatically become \code{FALSE}.}

\item{selectFUN}{Choose a variable selection function. Can specify either
\code{"varSelect"} or \code{"resample"} to use the corresponding functions.
If you want to use the \code{\link{resample}} function though, then it is
recommended to specify \code{selectFUN} as one of: \code{"stability",
  "split", "bootstrap"} in order to identify the specific method. If
\code{selectFUN = "resample"}, then it is recommended to add the
\code{sampMethod} argument to the call to \code{\link{mlGVAR}}.}

\item{subjectNets}{If \code{TRUE}, then subject-specific networks are fit for
all subjects and returned in the final output. Otherwise, can specify a
single value or a vector of values to represent which subjects to return
individual networks for -- specifically, the SUR network. One caveat is
that variable selection methods are not applied to these subject-specific
networks. Further modeling could be done using the output, however.}

\item{idvar}{Character string to indicate which variable contains the
participant identification numbers.}

\item{exogenous}{Logical. Indicates whether moderator variables should be
treated as exogenous or not. If they are exogenous, they will not be
modeled as outcomes/nodes in the network. If the number of moderators
reaches \code{k - 1} or \code{k}, then \code{exogenous} will automatically
be \code{FALSE}.}

\item{center}{Logical. Determines whether to mean-center the variables.}

\item{scale}{Logical. Determines whether to standardize the variables.}

\item{fixedType}{If logical, then any variable selection procedure specified
by \code{selectFUN} will not be applied to the SUR network. Alternatively,
a variable selection result, such as the output from either
\code{\link{varSelect}} or \code{\link{modSelect}}, can be supplied to
choose a specific constrained model in advance.}

\item{betweenType}{If logical, then any variable selection procedure
specified by \code{selectFUN} will not be applied to the SUR network.
Alternatively, a variable selection result, such as the output from either
\code{\link{varSelect}} or \code{\link{modSelect}}, can be supplied to
choose a specific constrained model in advance.}

\item{centerWithin}{Following the application of \code{center} and
\code{scale}, this determines whether to center variables within individual
subjects to create subject-centered values.}

\item{scaleWithin}{Following the application of \code{center} and
\code{scale}, this determines whether to scale variables within individual
subjects to create subject-standardized values.}

\item{rule}{Only applies to the between-subject network when a threshold is
supplied. The \code{"AND"} rule will only preserve edges when both
corresponding coefficients have p-values below the threshold, while the
\code{"OR"} rule will preserve an edge so long as one of the two
coefficients have a p-value below the supplied threshold.}

\item{threshold}{Logical or numeric. If \code{TRUE}, then a default value of
.05 will be set. Indicates whether a threshold should be placed on the
models at each iteration of the sampling. A significant choice by the
researcher.}

\item{verbose}{Logical. Determines whether to output progress bars and
messages in the console during the fitting process.}

\item{pcor}{See corresponding argument in the \code{\link{fitNetwork}}
function}

\item{fixedArgs}{A named list of arguments for the variable selection
function can be provided here, specifically those that are meant to be
applied to the SUR network estimation.}

\item{betweenArgs}{A named list of arguments for the variable selection
function can be provided for the between-subjects network.}

\item{bm}{Logical. Determines whether the same moderators are applied in the
between-subjects network. By default, the value of \code{m} only applies to
the SUR network. This allows one to decide whether or not to apply those
moderators in the between-subject network.}

\item{beepno}{Character string or numeric value to indicate which variable
(if any) encodes the survey number within a single day. Must be used in
conjunction with \code{dayno} argument.}

\item{dayno}{Character string or numeric value to indiciate which variable
(if any) encodes the survey number within a single day. Must be used in
conjunction with \code{beepno} argument.}

\item{deleteMissing}{Logical. Determines whether to automatically perform
listwise deletion if there are any missing values in the dataset.}

\item{...}{Additional arguments.}
}
\value{
\code{mlGVAR} objects
}
\description{
Fits a graphical vector autoregressive model to data containing multiple time
points measured for multiple individuals.
}
\details{
Uses a pseudo-mixed effects approach, wherein fixed effects are estimated and
random effects are approximated. See the work of Epskamp et al. (2018) for
more details on how these types of effects are estimated.
}
\examples{
\donttest{
fit1 <- mlGVAR(mlgvarDat, 'M')

fit2 <- mlGVAR(mlgvarDat, 'M', bm = TRUE) # Fit the same moderator in the between-subjects network

fit3 <- mlGVAR(mlgvarDat, 'M', selectFUN = 'varSelect')
}
}
\references{
Epskamp, S., Waldorp, L. J., Mottus, R., & Borsboom, B. (2018).
The gaussian graphical model in cross-sectional and time-series data.
Multivariate Behavioral Research. 53, 453-580.
}
\seealso{
\code{\link{mlGVARsim}, \link{lmerVAR}, \link{fitNetwork},
  \link{varSelect}, \link{resample}, \link{net}, \link{netInts},
  \link{plotNet}}
}
