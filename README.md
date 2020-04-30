# modnets
R package for moderated network models (UNDER CONSTRUCTION).

Designed to afford exploratory and confirmatory estimation of 3 types of moderated networks:

1. Cross-sectional moderated networks
	* Involves nodewise estimation of a GGM with higher-order interactions
2. Idiographic (temporal) moderated networks
	* Involves GLS estimation of a SUR time series model, as well as a concentration graph of the residuals.
3. Multi-level moderated networks
	* Uses one of two methods for estimation.
	* One is a two-step multilevel model, where fixed/random effects are estimated separately from between-subject effects
	* The other uses a formal multilevel moderated vector autoregressive model with `lmer`

## Downloading and using the package
Currently, the way to use the R package is as follows:

1. Download the folder `modnets` and place it in your working directory
2. Run: `source('modnets/functions.R')`
3. You now have access to all functions in the package!

The repository will be made accessible to download in R via `devtools` soon.

# Package Functions

## Primary functions
* The primary function used for the first two types of models is: `fitNetwork`. There are a variety of arguments and options that can be used for, e.g., penalized estimation, model specification, etc. These options are all contained within the `varSelect` algorithm.
* The primary functions used for the third model are: `mlGVAR` and `lmerVAR`, depending on which approach you wish to use.

## Model selection
* For model selection, you can use `varSelect` to employ either best-subset selection, the LASSO, ridge regression, or elastic net via `glmnet`, or the hierarchical LASSO via `glinternet`. These methods supports various information criteria as well as cross-validation for model selection.
* Additionally, you can use the `resample` function to use repeated subsampling or bootstrapping with the `varSelect` algorithm built in. 

## Stability \& power analyses
* Currently, these methods are not supported in the multilevel setting.
* For bootstrapping/edge-weight accuracy analysis, you can use the `bootNet` function (the name will be changed, given the popularity of the `bootnet` function by Sacha Epskamp that my code is based on).
* For case-dropping stability analysis, you can use `bootNet` while setting `caseDrop = TRUE`.
* For power analysis, you can use: `mnetPowerSim` to simulate data based on expected network structure(s).


Please contact trevorswanson222@gmail.com with any questions.
