## Resubmission
This is a resubmission. The following changes have been made:
* In the `lmerVAR` function, the call to `options(warn = -1)` has been removed. This also entailed removing the `warnings` argument.
* In the internal `margCIs` function, the `set.seed` input has been edited to be an argument of the function.
* Removed all instances of `\dontrun`, and updated examples for all functions.
* Added `\donttest` to all components of examples that take longer than 5s to run.
* Edited Description line of DESCRIPTION to include a link to the published article that describes the methods in this package.
* In the internal `simCor` function, the call to `options(warn = 2)` has been removed.
* In the README, the `nCores` argument of the `simNet` function has been removed, such that only a single core is used for simulating data.


## Test environments
* ubuntu 21.04, R 4.1.1
* win-builder (devel and release)
* local Windows 10 install, R 4.1.0

## R CMD check results
There were no ERRORs, WARNINGs.

There was 1 NOTE: New submission

## Downstream dependencies
There are currently no downstream dependencies for this package.
