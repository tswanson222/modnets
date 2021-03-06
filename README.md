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

A formal version of the package will be made accessible to download in R via `devtools` soon. The full repository may be found [here](https://github.com/tswanson222/modnets).

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
# Make sure 'modnets' folder is in the working directory
source('modnets/functions.R')
# Package is now loaded! 

### ================================================ ###
### ======= SIMULATE MODERATED NETWORK DATA ======== ###
### ================================================ ###
# Can simulate data with no moderators, or with one exogenous moderator
set.seed(999)
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

dat0 <- x$data[, -6]
dat1 <- x$data
# First, lets save an object that excludes the moderator (dat0)
# and save a second with the moderator (dat1)


### ================================================ ###
### =============== FITTING MODELS ================= ###
### ================================================ ###
# First, lets fit an unmoderated network, leaving out 'M' entirely
fit0 <- fitNetwork(data = dat0) 

# Next, lets fit a model that only includes 'M' as a covariate
fit1 <- fitNetwork(data = dat1, covariates = 6) 

# Now, lets fit the saturated model where 'M' moderates all edges in the network
fit2 <- fitNetwork(data = dat1, moderators = 6) 


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
plot(fit1, threshold = TRUE, predict = 'R2')
plot(fit2, threshold = TRUE, predict = 'adjR2')
# Using 'threshold = TRUE' is the same as 'threshold = .05'
# 'predict = TRUE' plots R2 values for each regression
# This can also be specified as a string, as shown

plot(fit2, threshold = TRUE, predict = fit0)
# This can also be used to visually compare networks
# Here, the light blue ring around each node shows
# the R2 for 'fit0', while the slightly darker piece 
# shows the increase in R2 that we see with 'fit2'

predictNet(fit2)
predictNet(fit2, fit1)
# We can extract these values using this function
# And can take the differences by supplying two networks
# Values for the second model are substracted by those from the first

plot(fit2, mnet = TRUE)
plot(fit2, threshold = TRUE, mnet = TRUE)
# 'mnet = TRUE' plots the exogenous moderator


### ============= VARIABLE SELECTION =============== ###
# These methods demonstrate the two-stage process for variable selection
# In the first stage, we use the data to select the active set of predictors
# In the second stage, we use those predictors to re-fit the models using OLS

### UNMODERATED NETWORKS
vars0 <- varSelect(data = dat0, criterion = 'BIC', method = 'glmnet')
vfit0 <- fitNetwork(data = dat0, type = vars0)
vfit1 <- fitNetwork(data = dat0, type = 'varSelect', criterion = 'BIC')
predictNet(vfit0, vfit1)
# In the first method, we use glmnet to perform variable selection for 
# each of the nodewise models. Then, we can subsequently include this in the
# 'fitNetwork' function. In the second approach, we can simply include everything
# in one command. We see that these produce the exact same models

vfit2 <- fitNetwork(data = dat0, type = 'varSelect', criterion = 'BIC', method = 'subset')
# We can also use best-subsets selection instead of the LASSO

predictNet(vfit2, vfit1)
# In this case, we see that best-subsets produced lower R2 for three nodes

vfit3 <- fitNetwork(data = dat0, type = 'varSelect', criterion = 'CV', seed = 1)
vfit3.1 <- fitNetwork(data = dat0, type = 'varSelect', criterion = 'CV', seed = 1)
vfit3.2 <- fitNetwork(data = dat0, type = 'varSelect', criterion = 'CV', seed = 99)
# We can also use cross-validation with glmnet (but not best-subsets)

predictNet(vfit3, vfit3.1)
predictNet(vfit3, vfit3.2)
# We see that setting a seed leads to reproducible results


### MODERATED NETWORKS
vars1 <- varSelect(data = dat1, m = 6, criterion = 'BIC', method = 'glinternet')
mfit1 <- fitNetwork(data = dat1, moderators = 6, type = vars1)
mfit2 <- fitNetwork(data = dat1, moderators = 6, type = 'varSelect', criterion = 'BIC')
predictNet(mfit1, mfit2)
# Again, we see that both methods produce the same model
# Creating the 'vars1' object separately can be useful when we wish
# to analyze the results from the variable selection process; plot outputs, obtain coefficients, etc.
# Also, all moderated networks use 'glinternet' as the selection method, and so it does not need to be specified

mfit2 <- fitNetwork(data = dat1, moderators = 6, type = 'varSelect', criterion = 'CV', seed = 1)
# We can use cross-validation with the glinternet algorithm as well


### ============== MODEL COMPARISON ================ ###
# Create a list of models we want to compare
fits <- list(fit0 = fit0, fit1 = fit1, fit2 = fit2, 
             vfit1 = vfit1, vfit2 = vfit2, vfit3 = vfit3,
             mfit1 = mfit1, mfit2 = mfit2)

modTable(fits = fits)
# Performs likelihood ratio tests comparing each model with every other

modTable(fits)$omnibus
# This shows us the final results. The 'LRT' column indicates
# the total number of times each model was selected across all tests
# We can see that 'fit2' (the saturated MNM) was selected across all tests
# The second-most selected was 'mfit2', which used glinternet with CV selection

modTable(fits = fits, nodes = TRUE)
# This does the same thing as above but at the nodewise level
```

More examples to be added soon.


Please contact t092s958@ku.edu with any questions.
