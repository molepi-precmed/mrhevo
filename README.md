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

```r
library(devtools)
devtools::install_github(repo="molepi-precmed/mrhevo", ref="main")
```

You can then run the example script:

```r
devtools::run_examples()
```

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
