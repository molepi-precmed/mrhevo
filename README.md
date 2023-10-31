# MR-Hevo -- Inference of causal effects by Mendelian randomization, marginalizing over distribution of pleiotropic effects

Inference of causal effects can be tackled like any other statistical problem, by computing the likelihood (or posterior distribution) of the parameter of interest (the causal effect) while marginalizing over the distribution of nuisance variables (in this case the pleiotropic effects).  If we can compute the posterior distribution of the parameter of interest, we can obtain the likelihood by dividing by the prior on that parameter.  

As the form of the distribution of pleiotropic effects over loci is unknown, any realistic statistical model has to specify a prior that encompassses a broad family of symmetric distributions ranging from a spike-and-slab to a Gaussian.  An initial implementation of this approach has been been described by [Berzuini et al (2020)](https://doi.org/10.1093/biostatistics/kxy027).  They specify a horseshoe prior for the pleiotropic effects, and generate the posterior distribution of all model parameters, including the causal effect parameter, by Markov chain Monte Carlo sampling.  An implementation that uses only summary statistics has been described by [Grant and Burgess (2020)](https://www.biorxiv.org/content/10.1101/2023.05.30.542988v1)

This method extends the likelihood-based approach 

1. to two-step Mendelian randomization, where step 1 uses only summary statistics for the effects of genetic instruments on exposure, and step 2 uses individual-level data to test the effects of these instruments on the outcome. 

2. to use a regularized horseshoe prior on shrinkage coefficients, as described by Piironen and Vehtari (2017).  This prior, widely known as the "Finnish horseshoe", has better computational properties than the original horseshoe.  In recognition of their contribution, the method is named `MR-Hevo` (hevo is the Finnish word for a horse). 

3. to generate classical maximum likelihood estimates and _p_-values for the causal effect. 


The motivation for this work was to develop a method to test formally for causality in [genome-wide aggregated _trans_- effects analysis](https://doi.org/10.1016/j.ajhg.2023.04.003), which aims to detect core genes for a disease or trait by testing for association with predicted _trans_- effects of SNPs on gene expression, aggregated over multiple QTLs.  With this approach, the genetic instruments are clumps of SNPs with trans- effects on the expression of a gene as transcript or circulating protein.  

This is work in progress.  We shall upload example datasets, vignettes and eventually an R package. 

