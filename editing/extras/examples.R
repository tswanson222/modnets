# Make sure 'modnets' folder is in the working directory
source('modnets/functions.R')
# Package is now loaded! 

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


### ============== MODEL COMPARISON ================ ###
# Create a list of models we want to compare
fits <- list(fit0 = fit0, fit1 = fit1, fit2 = fit2)

modTable(fits = fits)
# Performs likelihood ratio tests comparing the three models

modTable(fits = fits, nodes = TRUE)
# This does the same thing as above but at the nodewise level