# MR-Hevo -- Inference of causal effects by Mendelian randomization, marginalizing over distribution of pleiotropic effects

Mendelian randomization has been widely used to study causal effects of exposures (broadly defined to include behavioural traits, biomarkers and gene expression levels) on diseases.  The biggest methodological challenge is how to infer causality when some of the genetic instruments have direct (pleiotropic) effects on the outcome that are not mediated through the exposure under study.  These pleiotropic effects are not directly observed, and their distribution over the instruments is unknown.

Inference of causal effects can be tackled like any other statistical problem, by computing the likelihood (or posterior distribution) of the parameter of interest (the causal effect) while marginalizing over the distribution of nuisance variables (in this case the pleiotropic effects).  If we can compute the posterior distribution of the parameter of interest, we can obtain the likelihood by dividing by the prior on that parameter.

As the form of the distribution of pleiotropic effects over loci is unknown, any realistic statistical model has to specify a prior that encompassses a broad family of symmetric distributions ranging from a spike-and-slab to a Gaussian.  An initial implementation of this approach has been been described by [Berzuini et al (2020)](https://doi.org/10.1093/biostatistics/kxy027).  They specify a horseshoe prior for the pleiotropic effects, and generate the posterior distribution of all model parameters, including the causal effect parameter, by Markov chain Monte Carlo sampling.  An implementation that uses only summary statistics has been described by [Grant and Burgess (2020)](https://www.biorxiv.org/content/10.1101/2023.05.30.542988v1)

This method extends the likelihood-based approach

1. to two-step Mendelian randomization, where step 1 uses only summary statistics for the effects of genetic instruments on exposure, and step 2 uses individual-level data to test the effects of these instruments on the outcome.

2. to use a regularized horseshoe prior on shrinkage coefficients, as described by [Piironen and Vehtari (2017)](https://doi.org/10.1214/17-EJS1337SI).  This prior, known as the "Finnish horseshoe", has better computational properties than the original horseshoe.  On this basis, the method is named `MR-Hevo` (hevo is Finnish for a horse).

3. to generate classical maximum likelihood estimates and _p_-values for the causal effect.

The motivation for this work was to develop a method to test formally for causality in [genome-wide aggregated _trans_- effects analysis](https://doi.org/10.1016/j.ajhg.2023.04.003), which aims to detect core genes for a disease or trait by testing for association with predicted _trans_- effects of SNPs on gene expression, aggregated over multiple QTLs.  With this approach, the genetic instruments are clumps of SNPs with trans- effects on the expression of a gene as transcript or circulating protein.

## Guide

- A description of the statistical model is available on [theory page](https://github.com/molepi-precmed/mrhevo/blob/main/theorymethods.pdf)
- Also refer to the package [vignette](https://htmlpreview.github.io/?https://github.com/molepi-precmed/mrhevo/blob/main/vignette.html)

## Installation

To install current development version of the package from GitHub use:

1. For the production version:

```r
library(devtools)
devtools::install_github(repo="molepi-precmed/mrhevo", ref="package")
```

2. For the development version (recommended if facing Stan compilation issues):

```r
library(devtools)
devtools::install_github(repo="molepi-precmed/mrhevo", ref="main")
```

If you want to run the examples: 

```r
## get the source code
devtools::install_github("molepi-precmed/mrhevo", ref="main", build=FALSE)
setwd("path/to/local/clone/of/mrhevo")

## load the package
devtools::load_all()

## run the examples
devtools::run_examples()
```

## Example: Using summary statistics

This example demonstrates how to run MR-Hevo using summary statistics only (no individual-level data required).

```r
library(mrhevo)

## Load example summary statistics dataset (included in package)
coeffs <- readRDS(system.file("data/coeffs.RDS", package = "mrhevo"))

## Required inputs:
## - alpha_hat: effect of each instrument on the exposure
## - se.alpha_hat: standard error of alpha_hat
## - gamma_hat: effect of each instrument on the outcome  
## - se.gamma_hat: standard error of gamma_hat
alpha_hat <- coeffs$alpha_hat
se.alpha_hat <- coeffs$se.alpha_hat
gamma_hat <- coeffs$gamma_hat
se.gamma_hat <- coeffs$se.gamma_hat

cat("Number of genetic instruments:", length(alpha_hat), "\n")

## Path to Stan models (included in package)
model.dir <- system.file("stan", package = "mrhevo")

## Run MR-Hevo with summary statistics
## Using default priors: fraction_pleio = 0.5 (prior guess that 50% of instruments are pleiotropic)
fit <- run_mrhevo.sstats(
  alpha_hat = alpha_hat,
  se.alpha_hat = se.alpha_hat,
  gamma_hat = gamma_hat,
  se.gamma_hat = se.gamma_hat,
  fraction_pleio = 0.5,
  slab_scale = 0.2,
  priorsd_theta = 1,
  model.dir = model.dir
)

## Extract posterior samples
posterior_samples <- rstan::extract(fit)

## Get posterior summary
print(summary(fit, pars = "theta")$summary)

## Compute MLE and p-value from posterior
mle_result <- mle.se.pval(posterior_samples$theta, rep(1, length(posterior_samples$theta)))
print(mle_result)

## Plot IV estimates with MLE as slope line through origin
theta_mle <- mle_result$Estimate
p <- plot_iv_estimates(alpha_hat, se.alpha_hat, gamma_hat, se.gamma_hat, theta_mle, coeffs$qtlname)
print(p)
```

### Example results

Running the analysis on the included dataset (`data/coeffs.RDS`) with 43 genetic instruments yields:

**Conventional MR Estimators:**

| Estimator | Estimate | SE | z | p-value |
|-----------|----------|------|--------|---------|
| Inverse Variance Weighted (IVW) | -0.300 | 0.0196 | -15.31 | 6.1e-53 |

**MR-Hevo Bayesian Analysis:**

The MR-Hevo model with regularized horseshoe prior on pleiotropic effects provides a posterior distribution for the causal effect. After running the Stan model:

```
             mean     se_mean        sd       2.5%        50%      97.5%    n_eff    Rhat
theta     -0.337     0.0014    0.051     -0.436     -0.338     -0.233    1390     1.00
```

**MLE from posterior:**

```
   Estimate      SE        z    pvalue pvalue.formatted
1:  -0.326  0.0578 -5.652 1.58e-08           2e-08
```

The MR-Hevo analysis gives a posterior mean estimate of **-0.337** (95% CI: -0.436 to -0.233), indicating a strong causal effect of the exposure on the outcome. The MLE p-value of 1.6e-08 provides strong evidence for a causal effect, consistent with the conventional IVW estimator but with appropriate uncertainty quantification that accounts for potential pleiotropy.

**Plot of IV estimates:**

![IV estimates with MLE slope](./iv_estimates_plot.png)

### Interpreting the results

- **Posterior distribution**: The `theta` parameter represents the causal effect of the exposure on the outcome
- **MLE and p-value**: The `mle.se.pval()` function computes a maximum likelihood estimate and associated p-value by fitting a quadratic approximation to the log-posterior
- **Conventional MR estimators**: For comparison, the package also computes the inverse-variance weighted (IVW) estimator

### Choosing prior parameters

- `fraction_pleio`: Prior guess for the proportion of instruments with pleiotropic effects (0.05 to 0.95). Default is 0.5.
- `slab_scale`: Scale parameter for the regularized horseshoe slab component. Default is 0.2.
- `priorsd_theta`: Prior standard deviation for the causal effect theta. Default is 1 (weakly informative).

## Troubleshooting

This package is in early beta testing and things can go wrong! Proceed with caution.

Here we have listed some common errors that you may experience when installing or using the package. Please note this list is not exhaustive.

### System tmp directory is not writable

 `error in sink(type = "output") : invalid connection`

This is a Stan error related to system `/tmp` directory being not writable. This is a common case on HPC systems and linux servers.

* It is recommended to create a temporary directory in user's home directory (e.g. `/home/$USER/tmp`) and point R to it by creating `~/.Renviron` file with the following content:

    ```sh
    TMPDIR=/home/<username>/tmp
    TMP=/home/<username>/tmp
    TEMP=/home/<username>/tmp
  ```

### Stan cannot compile the model

* Rare Stan errors related to C compiler were also reported by users. This package was tested on Ubuntu 20.04 with `gcc version 9.4.0 (Ubuntu 9.4.0-1ubuntu1~20.04.2)` and R4.3. This package was not tested on MacOS or Windows. In case of any errors related to compilation of Stan models refer to Stan and Rstan documentation, help and troubleshooting guides.
* We noted that specifying custom C++ compilers and flags in `Makevars` often leads to errors, hence, refer to Stan and Rstan documentation if you need to use custom compilers and proceed with caution.

# Copyright

This code was developed by [Paul McKeigue](https://precmed.cphs.mvm.ed.ac.uk/pmckeigue), [Andrii Iakovliev](https://whimsial.github.io), [Buddhiprabha Erabadda](https://www.linkedin.com/in/buddhiprabha/) and [Athina Spiliopoulou](https://precmed.cphs.mvm.ed.ac.uk/athina/) and licensed under [GPL-3 license](https://www.gnu.org/licenses/gpl-3.0.txt).
