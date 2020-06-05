# modnets
R package for moderated network models (UNDER CONSTRUCTION).

Designed to afford exploratory and confirmatory estimation of 3 types of moderated networks:

1. Cross-sectional moderated networks
	* Involves nodewise estimation of a GGM with higher-order interactions
	* Can accomodate any combination of continuous and binary variables.
	* Nodewise regressions are fit using either OLS or logistic regression, depending on variable types.
2. Idiographic (temporal) moderated networks
	* Involves generalized least squares (GLS) estimation of multivariate time series model, as well as the inverse-covariance structure of the residuals. 
	* Currently only works for continuous variables, although exogenous moderators can be binary. 
	* Default estimation is seemingly unrelated regressions (SUR) via `systemfit`, but OLS is also available (unconstrained SUR estimates are equivalent to OLS estimates).
3. Multi-level moderated networks
	* Uses one of two methods for estimation.
	* One is a two-step multilevel model, where fixed/random effects are estimated separately from between-subject effects
	* The other uses a formal multilevel moderated vector autoregressive model with `lmer`
	* Only works for continous variables, although exogenous moderators can be binary.

Penalized estimators for each of these models are also available, such as the LASSO, ridge regression, elastic net, the (overlapping) group LASSO, and the hierarchical LASSO. Hyperparameter selection will be performed automatically based on either the AIC, BIC, EBIC, or cross-validation depending upon user input.


## Downloading and using the package
Currently, the way to use the R package is as follows:

1. Download the folder [modnets](https://github.com/tswanson222/modnets/tree/master/modnets) and place it in your working directory
2. In R, run: `source('modnets/functions.R')`
3. You now have access to all functions in the package!

The repository will be made accessible to download in R via `devtools` soon. Also, to access the datasets with the `settings` function you will need to download the [data](https://github.com/tswanson222/modnets/tree/master/data) folder and save it inside the [modnets](https://github.com/tswanson222/modnets/tree/master/modnets) folder.

Although we are loading the package via the `source` command, this command will not disturb the global environment. All components of the package will be stored in a separate environment called `.modnets`. This still gives the user access to all functions without cluttering the global environment.


# Package Functions

## Primary functions
* The primary function used for the first two types of models is: `fitNetwork`. There are a variety of arguments and options that can be used for, e.g., penalized estimation, model specification, etc.
* The primary functions used for the third model are: `mlGVAR` and `lmerVAR`, depending on which approach you wish to use.

## Model selection
* For model selection, you can use `varSelect` to employ either best-subset selection (via `leaps`), the LASSO, ridge regression, or elastic net (via `glmnet`), or the hierarchical LASSO (via `glinternet`). These methods support various information criteria as well as cross-validation for model selection, and are embedded within the `varSelect` function.
* As a note, all of the model selection procedures in `varSelect` operate on a sequential, nodewise basis.
* Additionally, you can use the `resample` function to use repeated subsampling or bootstrapping with the `varSelect` algorithm built in. 
* This latter method will take into account the actual model-fit values (such as those obtained in the GLS-driven SUR for temporal networks)

## Stability \& power analyses
* Currently, these methods are not supported in the multilevel setting.
* For bootstrapping/edge-weight accuracy analysis, you can use the `bootNet` function (the name will be changed, given the popularity of the `bootnet` function by Sacha Epskamp that my code is based on).
* For case-dropping stability analysis, you can use `bootNet` while setting `caseDrop = TRUE`.
* For power analysis, you can use: `mnetPowerSim` to simulate data based on expected network structure(s).


# Examples

## Cross-sectional moderated network
```r
source('modnets/functions.R')
# Package is now loaded! Make sure 'modnets' folder is in working directory.

### ================================================ ###
### ======= SIMULATE MODERATED NETWORK DATA ======== ###
### ================================================ ###
# Can simulate data with no moderators, or with one exogenous moderator
set.seed(123)
x <- simNet(N = 100, p = 5, m = TRUE, sparsity = .5)

str(x)
### Contents:
# x$data -------- 100x6 dataset, where 'M' is the moderator
# x$b1 ---------- true regression coefficients, where columns --> rows
# x$b2 ---------- true interaction coefficients, where (M * columns) --> rows
# x$intercepts -- true intercepts; defaults to 0
# x$m ----------- true mean of 'M'
# x$m1 ---------- coefficents for main effects of M on outcomes; default to 0

head(x$data)
print(x$b1)
print(x$b2)
print(x$intercepts)
print(x$m)
print(x$m1)


### ================================================ ###
### =============== FITTING MODELS ================= ###
### ================================================ ###
# First, lets fit an unmoderated network, leaving out 'M' entirely
fit0 <- fitNetwork(data = x$data[, -6]) 

# Next, lets fit a model that only includes 'M' as a covariate
fit1 <- fitNetwork(data = x$data, covariates = 6) 

# Now, lets fit the saturated model where 'M' moderates all edges in the network
fit2 <- fitNetwork(data = x$data, moderators = 6) 


### ================= PLOTTING ===================== ###
plot(fit0)
plot(fit1)
plot(fit2)
# We can plot each of these models to see the resultant undirected network

plot(fit0, threshold = .05)
plot(fit1, threshold = .05)
plot(fit2, threshold = .05)
# Plot only significant edges (p < threshold) of the network.

plot(fit0, threshold = TRUE, predict = TRUE)
plot(fit1, threshold = TRUE, predict = TRUE)
plot(fit2, threshold = TRUE, predict = TRUE)
plot(fit2, threshold = TRUE, predict = TRUE, con = 'adjR2')
# Using 'threshold = TRUE' is the same as 'threshold = .05'
# 'predict = TRUE' plots R2 values for each regression
# 'con = "adjR2"' uses adjusted R2

plot(fit2, mnet = TRUE)
plot(fit2, threshold = TRUE, mnet = TRUE)
# 'mnet = TRUE' plots the exogenous moderator


### ============== MODEL COMPARISON ================ ###
# Create a list of models we want to compare
fits <- list(fit0 = fit0, fit1 = fit1, fit2 = fit2)

modTable(fits = fits)
# Performs likelihood ratio tests comparing the three models

modTable(fits = fits, nodes = TRUE)
# This does the same thing as above but at the nodewise level
```

More examples to be added soon.


Please contact t092s958@ku.edu with any questions.
