#' Example of a complete MRHevo analysis of the exemplar dataset coeffs.dt,
#' with various tables and plots reporting the results. The example results will
#' be saved in mrhevo_example subdirectory of the current working directory.
#'
#' The package is loaded with the example dataset containing required columns:
#' qtlname, alpha_hat, se.alpha_hat, gamma_hat, se.gamma_hat. Other columns are
#' optional.
## -----------------------------------------------------------------------------
library(data.table)
library(rstan)
library(ggrepel)
library(ggplotify)
library(bayesplot)
library(cowplot)

info <- crayon::reset
note <- crayon::green
warn <- crayon::yellow
bold <- crayon::bold

header("Running mrhevo example script.")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

example.dir <- file.path(getwd(), "mrhevo_example")
dir.create(example.dir, showWarnings=FALSE)

## calculate coefficient ratios for each instrument
coeffs.dt <- get_coeffratios(coeffs.dt, use.delta=TRUE)

## set priors for Bayesian analysis
slab_scale <- 0.05
slab_df <- 2
fraction_pleio <- 0.2

## prior doesn't matter as we will divide posterior by it to get likelihood
priorsd_theta <- 1

## run Bayesian analysis to generate object of class stanfit
options(warn=1)
hevo.stanfit <-
    run_mrhevo.sstats(fraction_pleio=fraction_pleio,
                      alpha_hat=coeffs.dt$alpha_hat,
                      se.alpha_hat=coeffs.dt$se.alpha_hat,
                      gamma_hat=coeffs.dt$gamma_hat,
                      se.gamma_hat=coeffs.dt$se.gamma_hat,
                      slab_scale=slab_scale,
                      slab_df=slab_df,
                      priorsd_theta=1,
                      model.dir=devtools::package_file())
options(warn=2)

## get sampler diagnostics
num.divergent <- get_num_divergent(hevo.stanfit)
num.maxtreedepth <- get_num_max_treedepth(hevo.stanfit)

## create and save pairs plot of posterior samples
p.pairs <- mcmc_pairs(hevo.stanfit,
                      pars=c("theta", "log_c", "log_tau", "f", "lp__"),
                      off_diag_args=list(size = 0.75))
ggsave(file.path(example.dir, "pairsplot.png"), p.pairs)

## create and save trace plot of posterior samples
p.traceplot <-  traceplot(hevo.stanfit,
                          pars=c("theta", "log_c", "log_tau", "f", "lp__"),
                          inc_warmup=FALSE)
ggsave(file.path(example.dir, "traceplot.png"), p.traceplot)

## get maximum likelihood estimate of theta from posterior density
theta.samples <-  unlist(extract(hevo.stanfit, pars="theta"))
prior.theta <- dnorm(theta.samples, mean=0, sd=priorsd_theta)
options(warn=1)
mle.theta <- mle.se.pval(x=theta.samples, prior=prior.theta)
options(warn=2)
mle.theta[, Estimator := "Marginalize over direct effects"]

## plot posterior and log-likelihood
options(warn=1)
p.bayesian.loglik <- mle.se.pval(x=theta.samples, prior=prior.theta,
                                 return.asplot=TRUE)
options(warn=2)
ggsave(file.path(example.dir, "posterior_loglik.png"), p.bayesian.loglik)

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

msg(info, "Table of the posterior summaries for parameters.\n")
print(fit.coeffs)

## plot shrinkage coefficients
p.shrinkage <- ggplot(data=fit.coeffs[grep("kappa", variable)],
                      aes(y=qtlname, x=mean, xmin=`10%`, xmax=`90%`)) +
    geom_pointrange()
ggsave(file.path(example.dir, "shrinkageplot.png"), p.shrinkage)

p.beta <- ggplot(data=fit.coeffs[grep("beta", variable)],
                      aes(y=qtlname, x=mean, xmin=`10%`, xmax=`90%`)) +
    geom_pointrange() +
    geom_vline(xintercept=0, linetype="dotted")
ggsave(file.path(example.dir, "betaplot.png"), p.beta)

## get MR "estimators": weighted mean, weighted median,
## penalized weighted median and append MRHevo estimate
estimators <- get_estimatorsMR(coeffs.dt)
estimators <- rbind(estimators, mle.theta, fill=TRUE)

## plot coefficients and show estimators as slopes of lines through origin
p.coeffs <- ggplot(coeffs.dt,
                       aes(x=alpha_hat, y=gamma_hat)) +
    geom_point(aes(size=size.theta_IV), alpha=0.8) +
    scale_size(guide="none") +
    ggrepel::geom_text_repel(aes(label=qtlname),
                             force=5, size=2.5, fontface="italic", color="blue") +
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
ggsave(file.path(example.dir, "coeffsplots.png"), p.coeffs)

note.msg <- sprintf("Example script has finished.\n Results were written to %s",
                    example.dir)
msg(note, note.msg)
