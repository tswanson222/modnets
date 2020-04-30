# modnets
R package for moderated network models (under construction).

Designed to afford exploratory and confirmatory estimation of 3 types of moderated networks:

1. Cross-sectional moderated networks
	* Involves nodewise estimation of a GGM with higher-order interactions
2. Idiographic (temporal) moderated networks
	* Involves GLS estimation of a SUR time series model, as well as a concentration graph of the residuals.
3. Multi-level moderated networks
	* Uses one of two methods for estimation.
	* First method is a two-step multilevel model, where fixed/random effects are estimated separately from between-subject effects
	* Second method uses a formal multilevel moderated vector autoregressive model with `lmer`

The repository will be made accessible to download in R via `devtools` soon.

Please contact trevorswanson222@gmail.com with any questions.
