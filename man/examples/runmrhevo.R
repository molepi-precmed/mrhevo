## Complete MRHevo analysis of the exemplar dataset coeffs.dt.
## Results are written to a temporary directory.
## The NumPyro MCMC call is wrapped in \donttest{} because it requires a
## Python environment and takes more than 5 seconds to run.

library(data.table)

## calculate coefficient ratios for each instrument
coeffs.dt <- get_coeffratios(coeffs.dt, use.delta=TRUE)

## get IVW and other classical MR estimators
estimators <- get_estimatorsMR(coeffs.dt)
print(estimators)

## plot IV estimates (fast; always run)
p.coeffs <- plot_iv_estimates(
    alpha_hat    = coeffs.dt$alpha_hat,
    se.alpha_hat = coeffs.dt$se.alpha_hat,
    gamma_hat    = coeffs.dt$gamma_hat,
    se.gamma_hat = coeffs.dt$se.gamma_hat,
    theta        = estimators[Estimator == "IVW", Estimate],
    qtlname      = coeffs.dt$qtlname)

\donttest{
## run Bayesian analysis using NumPyro (requires install_mrhevo_python())
example.dir  <- file.path(tempdir(), "mrhevo_example")
dir.create(example.dir, showWarnings = FALSE)
model_path   <- system.file("python", "mrhevo_pyro.py", package = "mrhevo")

hevo.fit <- run_mrhevo.numpyro(
    alpha_hat    = coeffs.dt$alpha_hat,
    se.alpha_hat = coeffs.dt$se.alpha_hat,
    gamma_hat    = coeffs.dt$gamma_hat,
    se.gamma_hat = coeffs.dt$se.gamma_hat,
    fraction_pleio = 0.2,
    slab_scale   = 0.05,
    slab_df      = 2,
    priorsd_theta = 1,
    model_path   = model_path,
    num_warmup   = 500,
    num_samples  = 1000,
    num_chains   = 4)

## MLE and SE from posterior
theta.samples <- hevo.fit$posterior$theta
prior.theta   <- dnorm(theta.samples, mean = 0, sd = 1)
mle.theta     <- mle.se.pval(x = theta.samples, prior = prior.theta)
mle.theta[, Estimator := "Marginalise over direct effects"]
print(mle.theta)

## save plots to temp directory
ggsave(file.path(example.dir, "posterior_loglik.png"),
       mle.se.pval(x = theta.samples, prior = prior.theta, return.asplot = TRUE))
ggsave(file.path(example.dir, "pairsplot.png"),
       plot_posterior_pairs(hevo.fit, pars = c("log_tau", "log_eta", "f")))
ggsave(file.path(example.dir, "kappa_hist.png"),
       plot_kappa_hist(hevo.fit))

## --- Optional: Stan backend (requires install_mrhevo_stan()) ----------------
## hevo.stanfit <- run_mrhevo_stan(
##     alpha_hat    = coeffs.dt$alpha_hat,
##     se.alpha_hat = coeffs.dt$se.alpha_hat,
##     gamma_hat    = coeffs.dt$gamma_hat,
##     se.gamma_hat = coeffs.dt$se.gamma_hat,
##     fraction_pleio = 0.2,
##     slab_scale   = 0.05,
##     slab_df      = 2,
##     priorsd_theta = 1)
}
