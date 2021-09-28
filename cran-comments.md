## Resubmission
This is a resubmission. The following changes have been made:
\itemize{
  \item{In the \code{lmerVAR} function, the call to \code{options(warn = -1)} has been removed. This also entailed removing the \code{warnings} argument from the function.}
  \item{In the internal \code{margCIs} function, the \code{set.seed} input has been edited to be an argument of the function.}
  \item{Removed all instances of \\dontrun, and updated examples for all functions.}
}

## Test environments
* ubuntu 21.04, R 4.1.1
* win-builder (devel and release)
* local Windows 10 install, R 4.1.0

## R CMD check results
There were no ERRORs, WARNINGs.

There was 1 NOTE: New submission

## Downstream dependencies
There are currently no downstream dependencies for this package.
