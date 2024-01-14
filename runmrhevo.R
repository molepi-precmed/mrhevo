
## this script runs a complete MRHevo analysis of the exemplar dataset coeffs.dt, with various tables and plots reporting the results

## the functions are in ./functions.mrhevo.R with documentation in correct format for R package
## I suggest making a package so that running the Stan model creates a a class mrhevo.stanfit. 
## This should contain the stanfit object hevo.stanfit and the data.table object coeffs.dt

## all the chunks of code called to generate tables and plots should be methods for this object
## sampler diagnostics, pairs plot, trace plot, max likelihood estimate, plot posterior density, posterior summaries, plot shrinkage coefficients, plot coeffs

library(data.table)
library(rstan)
library(ggrepel)
library(ggplotify)

source("../mrhevo/functions.mrhevo.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## exemplar dataset to be included in package
## dataset should have columns named qtlname, alpha_hat, se.alpha_hat, gamma_hat, se.gamma_hat
## other columns are optional
coeffs.dt <- readRDS("./coeffs_ADIPOQ_UKBB.RDS")
coeffs.dt <- readRDS("../mrhevo/coeffs_CLC_asthma_UKBB.RDS")
#coeffs.dt <- coeffs.dt[minpvalue < 1E-8]
## calculate coefficient ratios for each instrument
coeffs.dt <- get_coeffratios(coeffs.dt, use.delta=TRUE)

## set priors for Bayesian analysis
slab_scale <- 0.05
slab_df <- 2
fraction_pleio <- 0.2

priorsd_theta <- 1  # prior doesn't matter as we will divide posterior by it to get likelihood

## run Bayesian analysis to generate object of class stanfit
options(warn=1)
hevo.stanfit <-
    run_mrhevo.sstats(fraction_pleio=fraction_pleio,
   # run_mrhevo.fixedtau(tau=1E-6,
                      alpha_hat=coeffs.dt$alpha_hat,
                      se.alpha_hat=coeffs.dt$se.alpha_hat,
                      gamma_hat=coeffs.dt$gamma_hat,
                      se.gamma_hat=coeffs.dt$se.gamma_hat,
                      slab_scale=slab_scale,
                      slab_df=slab_df, 
                      priorsd_theta=1)
options(warn=2)

## sampler diagnostics
num.divergent <- get_num_divergent(hevo.stanfit)
num.maxtreedepth <- get_num_max_treedepth(hevo.stanfit)

## pairs plot of posterior samples
p.pairs <- ggplotify::as.ggplot(function() {
    pairs(hevo.stanfit,
          pars=c("theta", "log_c", "log_tau", "f", "lp__"))
})
p.pairs

## trace plot of posterior samples
p.traceplot <-  traceplot(hevo.stanfit,
                          pars=c("theta", "log_c", "log_tau", "f", "lp__"),
                          inc_warmup=FALSE)

## maximum likelihood estimate of theta from posterior density
theta.samples <-  unlist(extract(hevo.stanfit, pars="theta"))
prior.theta <- dnorm(theta.samples, mean=0, sd=priorsd_theta)
mle.theta <- mle.se.pval(x=theta.samples, prior=prior.theta)
mle.theta[, Estimator := "Marginalize over direct effects"]

## plot posterior and log-likelihood
options(warn=1)
p.bayesian.loglik <-  mle.se.pval(x=theta.samples, prior=prior.theta, return.asplot=TRUE)
options(warn=2)

## tabulate posterior summaries for global parameters
pars.forsummary <- c("theta", "log_c", "log_tau", "f", "lp__", "beta", "kappa")

fit.coeffs <- summary(hevo.stanfit,
                      pars=pars.forsummary, 
                      probs=c(0.1, 0.9))$summary 
fit.coeffs <- as.data.table(round(fit.coeffs, 3), keep.rownames="variable")
fit.coeffs <- fit.coeffs[, .(variable, mean, `10%`, `90%`, n_eff, Rhat)]
fit.coeffs[, n_eff := round(n_eff, 2)]

qtlnames <- coeffs.dt$qtlname
fit.coeffs[grep("beta\\[", variable), qtlname := ..qtlnames]
fit.coeffs[grep("kappa\\[", variable), qtlname := ..qtlnames]
fit.coeffs[qtlname=="", qtlname := variable]  
fit.coeffs

## plot shrinkage coefficients
p.shrinkage <- ggplot(data=fit.coeffs[grep("kappa", variable)],
                      aes(y=qtlname, x=mean, xmin=`10%`, xmax=`90%`)) +
    geom_pointrange()
p.shrinkage

p.beta <- ggplot(data=fit.coeffs[grep("beta", variable)],
                      aes(y=qtlname, x=mean, xmin=`10%`, xmax=`90%`)) +
    geom_pointrange() + 
    geom_vline(xintercept=0, linetype="dotted")
p.beta

## get MR "estimators": weighted mean, weighted median, penalized weighted median
## append MRHevo estimate
estimators <- get_estimatorsMR(coeffs.dt)
estimators <- rbind(estimators, mle.theta, fill=TRUE)

## plot coefficients and show estimators as slopes of lines through origin
p.coeffs <- ggplot(coeffs.dt,
                       aes(x=alpha_hat, y=gamma_hat)) + 
    geom_point(aes(size=size.theta_IV), alpha=0.8) +
    scale_size(guide="none") + 
    ggrepel::geom_text_repel(aes(label=qtlname), force=5, size=2.5, fontface="italic", color="blue") + 
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) + 
    scale_y_continuous(limits = c(min(c(coeffs.dt$gamma_hat, 0)),
                                  max(c(coeffs.dt$gamma_hat, 0))),
                       expand =  expansion(mult = c(0.1, 0.1))) + 
    geom_abline(data=estimators,
                aes(slope=Estimate, intercept=rep(0, nrow(estimators)), 
                    linetype=Estimator)) +
    labs(linetype="Estimate of causal effect as slope of line through origin") + 
    theme(legend.position="top") +
    theme(legend.direction="vertical") +
    theme(legend.box="horizontal") +
    xlab(paste("Effect of genetic instrument on exposure")) +
    ylab(paste("Effect on outcome"))
options(warn=1)
p.coeffs
options(warn=2)
