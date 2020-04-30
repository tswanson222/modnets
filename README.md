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

# Examples

## Cross-sectional moderated network
```
source('modnets/functions.R')
# Package is now loaded! Make sure 'modnets' folder is in working directory.

### SIMULATE MODERATED NETWORK DATA
set.seed(1)
x <- simNet(N = 100, p = 5, m = TRUE, sparsity = .5)

head(x$data)
# A 100x6 dataset, where the variable 'M' is the moderator, and there are 5 outcomes.
# The 'x' object also contains elements that describe the true model parameters


### FITTING MODELS
# First, lets fit an unmoderated network, leaving out 'M' entirely
fit0 <- fitNetwork(data = x$data[, -6]) 

# Next, lets fit a model that only includes 'M' as a covariate
fit1 <- fitNetwork(data = x$data, covariates = 6) 

# Now, lets fit the saturated model where 'M' moderates all edges in the network
fit2 <- fitNetwork(data = x$data, moderators = 6) 

# Create a list of models we want to compare
fits <- list(fit0 = fit0, fit1 = fit1, fit2 = fit2)

plot(fit0)
plot(fit1)
plot(fit2)
# We can plot each of these models to see the resultant undirected network

plot(fit0, threshold = TRUE)
plot(fit1, threshold = TRUE)
plot(fit2, threshold = .05)
# Plot only significant edges of the network. threshold = TRUE defaults to threshold = .05

modTable(fits = fits)
# Performs likelihood ratio tests comparing the three models.
# Right-most column displays how many times each model was selected out of all the comparisons

modTable(fits = fits, nodes = TRUE)
# This does the same thing as above but at the nodewise level

```

More examples to be added soon.


Please contact t092s958@ku.edu with any questions.
