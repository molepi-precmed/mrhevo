## run_mrhevo_pdcd1_t1d.R
## Run MR-Hevo (NumPyro backend) for the PDCD1 -> T1D analysis.
## Saves fit object and derived summaries to mrhevo_pdcd1_t1d_fit.rds.
## After this script completes, re-render pdcd1_t1d_mrhevo.Rmd to include results.

library(data.table)
devtools::load_all(".", quiet = TRUE)

## ---- Load and merge coefficients ----
alpha_dt <- readRDS("alpha_pdcd1_OID00791.rds")
gamma_dt  <- readRDS("gamma_t1d.rds")

coeffs.dt <- merge(
    alpha_dt[, .(qtlname, chr, locus_start,
                 n_snps, eff_rank, alpha_hat, se.alpha_hat)],
    gamma_dt[, .(qtlname, n_snps_gamma, gamma_hat, se.gamma_hat)],
    by = "qtlname"
)
setorder(coeffs.dt, chr, locus_start)
coeffs.dt <- get_coeffratios(coeffs.dt, use.delta = TRUE)

message(sprintf("Running MR-Hevo on %d instruments", nrow(coeffs.dt)))

## ---- Run NumPyro sampler ----
priorsd_theta  <- 1.0
fraction_pleio <- 0.5

fit.numpyro <- run_mrhevo.numpyro(
    alpha_hat      = coeffs.dt$alpha_hat,
    se.alpha_hat   = coeffs.dt$se.alpha_hat,
    gamma_hat      = coeffs.dt$gamma_hat,
    se.gamma_hat   = coeffs.dt$se.gamma_hat,
    fraction_pleio = fraction_pleio,
    slab_scale     = 0.2,
    slab_df        = 2,
    priorsd_theta  = priorsd_theta
)

## ---- Derive MLE of theta ----
theta.samples <- as.vector(fit.numpyro$posterior$theta)
prior.theta   <- dnorm(theta.samples, mean = 0, sd = priorsd_theta)
mle.theta     <- mle.se.pval(x = theta.samples, prior = prior.theta)

message(sprintf("theta MLE = %.4f (SE %.4f, p = %s)",
                mle.theta$Estimate, mle.theta$SE, mle.theta$pvalue.formatted))

## ---- Save ----
out <- list(
    fit         = fit.numpyro,
    coeffs.dt   = coeffs.dt,
    mle.theta   = mle.theta,
    priorsd_theta = priorsd_theta
)
saveRDS(out, "mrhevo_pdcd1_t1d_fit.rds")
message("Saved: mrhevo_pdcd1_t1d_fit.rds")
message("Now re-render pdcd1_t1d_mrhevo.Rmd to include posterior results.")
