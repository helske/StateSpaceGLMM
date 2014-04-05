StateSpaceGLMM: Inference of Generalized Linear Mixed Models Using State Space Representation
======

StateSpaceGLMM is an R package for analysis of generalized linear mixed models via state space approach. StateSpaceGLMM is a front end for R package KFAS, which focuses on GLMMs without the complexity of KFAS. Current features include

* Support for Gaussian, Poisson, binomial, gamma and negative binomial distributions, and combinations of these. For example, some groups can be modelled as Poisson and others as negative binomial. Note that at least for now only one grouping variable for random effects is supported.
* REML and ML estimation of unknown model parameters.
* Prediction of missing and "future" observations.
* Correlating and independent random effects.

In future, option to add (possibly correlating) structural time series components will be added.

