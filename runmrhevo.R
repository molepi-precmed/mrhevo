
## this script runs a complete MRHevo analysis of the exemplar dataset coeffs.dt, with various tables and plots reporting the results

## the functions are in ./functions.mrhevo.R with documentation in correct format for R package
## I suggest making a package so that running the Stan model creates a a class mrhevo.stanfit. 
## This should contain the stanfit object hevo.stanfit and the data.table object coeffs.dt

## all the chunks of code called to generate tables and plots should be methods for this object
## sampler diagnostics, pairs plot, trace plot, max likelihood estimate, plot posterior density, posterior summaries, plot shrinkage coefficients, plot coeffs

library(data.table)
library(ggrepel)

source("../mrhevo/functions.mrhevo.R")

## exemplar dataset to be included in package
## dataset should have columns named qtlname, alpha_hat, se.alpha_hat, gamma_hat, se.gamma_hat
## other columns are optional
#coeffs.dt <- readRDS("./coeffs_ADIPOQ_UKBB.RDS")
coeffs.dt <- readRDS("./coeffs_LPL_t2d.diag_1704571.RDS")[qtl_type=="trans"]
#coeffs.dt <- coeffs.dt[minpvalue < 1E-8]
## calculate coefficient ratios for each instrument
coeffs.dt <- get_coeffratios(coeffs.dt, use.delta=TRUE)

## set priors for Bayesian analysis
slab_scale <- 0.05
slab_df <- 2
fraction_pleio <- 0.2

priorsd_theta <- 1  # prior doesn't matter as we will divide posterior by it to get likelihood

f.errorsq <- function(tau, info, f.target) {
    asq <- info * tau^2
    f.expected <- 1 - mean(1 / (1 + asq))
    return((f.target - f.expected)^2)
}

tau0 <- optimize(f.errorsq, interval=c(0, 1), info=coeffs.dt$inv.var,
                 f.target=fraction_pleio, 
                 lower=0, upper=1)$minimum

alpha_hat <- coeffs.dt$alpha_hat
se.alpha_hat <- coeffs.dt$se.alpha_hat
gamma_hat <- coeffs.dt$gamma_hat
se.gamma_hat <- coeffs.dt$se.gamma_hat
info <- coeffs.dt$inv.var 
tau0 <- tau0 
slab_scale <- slab_scale
slab_df <- slab_df 
priorsd_theta <- 1.0

save(alpha_hat,
     se.alpha_hat,
     gamma_hat,
     se.gamma_hat,
     info, 
     tau0, 
     slab_scale,
     slab_df, 
     priorsd_theta,
     file="mrhevo_example.RData")

regularized=TRUE

## run Bayesian analysis to generate object of class stanfit
dense_mass=FALSE
target_accept_prob=0.95
num_warmup=500L
source("mrhevo_pyro.R")

## sampler diagnostics
#num.divergent <- get_num_divergent(hevo.stanfit)
#num.maxtreedepth <- get_num_max_treedepth(hevo.stanfit)

if(regularized) {
    pars.global <- c("theta", "log_c", "log_tau", "f")
} else {
    pars.global <- c("theta", "log_tau", "f")
}

global.summarystats <- mrhevo.summaries.dt[param %in% pars.global]
global.summarystats[, index := NULL]
global.summarystats[, lapply(.SD, round, 3), by=param, .SDcols=2:4]

## maximum likelihood estimate of theta from posterior density
theta.samples <-  mrhevo_samples.list[["theta"]]
prior.theta <- dnorm(theta.samples, mean=0, sd=priorsd_theta)
mle.theta <- mle.se.pval(x=theta.samples, prior=prior.theta)
mle.theta[, Estimator := "Marginalize over direct effects"]

## plot posterior and log-likelihood
options(warn=1)
p.bayesian.loglik <-  mle.se.pval(x=theta.samples, prior=prior.theta, return.asplot=TRUE)
options(warn=2)

J <- nrow(coeffs.dt)

## plot shrinkage coefficients
kappa.coeffs <- mrhevo.summaries.dt[param=="kappa"]
kappa.coeffs <- data.table(coeffs.dt[, .(qtlname)], kappa.coeffs)
kappa.coeffs[qtlname=="", qtlname := paste0("QTL_", .I)]
p.kappa <- ggplot(data=kappa.coeffs,
                      aes(y=qtlname, x=`50%`, xmin=`10%`, xmax=`90%`)) +
    geom_pointrange()
p.kappa

beta.coeffs <- mrhevo.summaries.dt[param=="beta"]
beta.coeffs <- data.table(coeffs.dt[, .(qtlname)], beta.coeffs)
beta.coeffs[qtlname=="", qtlname := paste0("QTL_", .I)]
p.beta <- ggplot(data=beta.coeffs,
                      aes(y=qtlname, x=`50%`, xmin=`10%`, xmax=`90%`)) +
    geom_pointrange()
p.beta

## get MR "estimators": weighted mean, weighted median, penalized weighted median
## append MRHevo estimate
estimators <- get_estimatorsMR(coeffs.dt, use.delta=TRUE)
estimators <- rbind(estimators, mle.theta, fill=TRUE)

## plot raw coefficients and show estimators as slopes of lines through origin
p.coeffs <- ggplot(coeffs.dt,
                       aes(x=alpha_hat, y=gamma_hat)) + 
    geom_point(aes(size=sqrt(max(inv.var) / inv.var)), alpha=0.8) +
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
    theme(legend.direction="vertical") 
    theme(legend.box="horizontal") +
    xlab(paste("Summary statistic for effect of genetic instrument on exposure")) +
    ylab(paste("Summary statistic for effect on outcome"))
options(warn=1)
p.coeffs
options(warn=2)

## get lp for each model
