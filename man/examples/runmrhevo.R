#' Example of a complete MRHevo analysis of the exemplar dataset coeffs.dt,
#' with various tables and plots reporting the results. The example results will
#' be saved in mrhevo_example subdirectory of the current working directory.
#'
#' The package is loaded with the example dataset containing required columns:
#' qtlname, alpha_hat, se.alpha_hat, gamma_hat, se.gamma_hat. Other columns are
#' optional.
## -----------------------------------------------------------------------------
library(data.table)
library(ggrepel)
library(cowplot)

info <- crayon::reset
note <- crayon::green
warn <- crayon::yellow
bold <- crayon::bold

header("Running mrhevo example script.")

example.dir <- file.path(getwd(), "mrhevo_example")
dir.create(example.dir, showWarnings=FALSE)

## calculate coefficient ratios for each instrument
coeffs.dt <- get_coeffratios(coeffs.dt, use.delta=TRUE)

## set priors for Bayesian analysis
fraction_pleio <- 0.2

## prior doesn't matter as we will divide posterior by it to get likelihood
priorsd_theta <- 1

## Path to Python model (included in package)
model_path <- system.file("python", "mrhevo_pyro.py", package="mrhevo")

## run Bayesian analysis using NumPyro (default backend)
options(warn=1)
hevo.fit <- run_mrhevo.numpyro(
    alpha_hat=coeffs.dt$alpha_hat,
    se.alpha_hat=coeffs.dt$se.alpha_hat,
    gamma_hat=coeffs.dt$gamma_hat,
    se.gamma_hat=coeffs.dt$se.gamma_hat,
    fraction_pleio=fraction_pleio,
    slab_scale=0.05,
    slab_df=2,
    priorsd_theta=priorsd_theta,
    model_path=model_path,
    num_warmup=500,
    num_samples=1000,
    num_chains=4)
options(warn=2)

## get maximum likelihood estimate of theta from posterior density
theta.samples <- hevo.fit$posterior$theta
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

## create and save pairs plot of posterior samples
p.pairs <- plot_posterior_pairs(hevo.fit,
                                pars=c("log_tau", "log_eta", "f"))
ggsave(file.path(example.dir, "pairsplot.png"), p.pairs)

## create and save histogram of kappa shrinkage coefficients
p.kappa <- plot_kappa_hist(hevo.fit)
ggsave(file.path(example.dir, "kappa_hist.png"), p.kappa)

msg(info, "NumPyro posterior summary:\n")
cat("theta mean:", mean(theta.samples), "\n")
cat("theta sd:  ", sd(theta.samples), "\n")
cat("f (pleiotropy fraction) mean:", mean(hevo.fit$posterior$f), "\n")
print(mle.theta)

## get MR "estimators": inverse-variance weighted and append MRHevo estimate
estimators <- get_estimatorsMR(coeffs.dt)
estimators <- rbind(estimators, mle.theta, fill=TRUE)

## plot coefficients and show estimators as slopes of lines through origin
theta_mle <- mle.theta$Estimate
p.coeffs <- plot_iv_estimates(
    alpha_hat=coeffs.dt$alpha_hat,
    se.alpha_hat=coeffs.dt$se.alpha_hat,
    gamma_hat=coeffs.dt$gamma_hat,
    se.gamma_hat=coeffs.dt$se.gamma_hat,
    theta=theta_mle,
    qtlname=coeffs.dt$qtlname)
ggsave(file.path(example.dir, "coeffsplot.png"), p.coeffs)

## --- Optional: Stan backend (requires install_mrhevo_stan()) ----------------
## model.dir <- system.file("stan", package="mrhevo")
## options(warn=1)
## hevo.stanfit <- run_mrhevo.sstats(
##     alpha_hat=coeffs.dt$alpha_hat,
##     se.alpha_hat=coeffs.dt$se.alpha_hat,
##     gamma_hat=coeffs.dt$gamma_hat,
##     se.gamma_hat=coeffs.dt$se.gamma_hat,
##     fraction_pleio=fraction_pleio,
##     slab_scale=0.05,
##     slab_df=2,
##     priorsd_theta=priorsd_theta,
##     model.dir=model.dir)
## options(warn=2)

note.msg <- sprintf("Example script has finished.\n Results were written to %s",
                    example.dir)
msg(note, note.msg)
