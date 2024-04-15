library(data.table)
library(ggplot2)
library(reticulate)
## source /home/pmckeigue/.virtualenvs/r-reticulate/bin/activate

summary.param <- function(param.samples) {
    if(length(dim(param.samples)) > 1) {
        param.summary <- as.data.table(t(apply(param.samples, 2, quantile,
                                               probs=c(0.5, 0.1, 0.9))))
    } else {
        param.summary <- as.data.table(as.list(quantile(param.samples,
                                                        probs=c(0.5, 0.1, 0.9))))
    }
    return(param.summary)
}

np <- import("numpy", as="np")

#load("mrhevo_example.RData")
alpha_hat <- np$array(alpha_hat)
gamma_hat <- np$array(gamma_hat)
se_alpha_hat <- np$array(se.alpha_hat)
se_gamma_hat <- np$array(se.gamma_hat)
info <- np$array(info)
scale_global <- tau0

source_python("../mrhevo/mrhevo_pyro.py")

run.nullhevo <- FALSE
if(run.nullhevo) { # fit null model
    nullhevo_model <- nullhevo(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info)
    nullhevo_mcmc <- setup_sampler(nullhevo_model, dense_mass=dense_mass, target_accept_prob=target_accept_prob)

    start_time <- Sys.time()
    nullhevo_mcmc <- run_sampler(nullhevo_mcmc, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info)
    elapsed.time <- Sys.time() - start_time
    nullhevo_idata <- get_idata(nullhevo_mcmc)
    nullhevo_diagnostics <- get_diagnostics(nullhevo_idata)
    nullhevo_diagnostics.names <- sapply(nullhevo_diagnostics, function(x) np$asarray(x[[1]]))
    ## names are diverging, lp
    num_divergent <- length(which(np$asarray(nullhevo_diagnostics[[1]][[2]])))
    cat("Fitting null model took ", elapsed.time, "minutes with", num_divergent, "divergences\n")

    nullhevo_samples <- get_null_mcmc_samples(nullhevo_mcmc)
    nullhevo_samples.names <- sapply(nullhevo_samples, function(x) np$asarray(x[[1]]))
    nullhevo_samples.list <- sapply(nullhevo_samples, function(x) np$asarray(x[[2]]))
    names(nullhevo_samples.list) <- nullhevo_samples.names
}

run.mrhorse <- FALSE
if(run.mrhorse) {
    mrhorse_model <- mrhorse(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info)
    mrhorse_mcmc <- setup_sampler(mrhorse_model)
    start_time <- Sys.time()
    mrhorse_mcmc <- run_sampler(mrhorse_mcmc, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info)
    elapsed.time <- Sys.time() - start_time
    mrhorse_idata <- get_idata(mrhorse_mcmc)
    mrhorse_diagnostics <- get_diagnostics(mrhorse_idata)
    mrhorse_diagnostics.names <- sapply(mrhorse_diagnostics, function(x) np$asarray(x[[1]]))
    num_divergent <- length(which(np$asarray(mrhorse_diagnostics[[1]][[2]])))
    lp <- np$asarray(mrhorse_diagnostics[[2]][[2]]) # dimension chains x numdraws
    cat("Fitting unregularized horseshoe took ", elapsed.time, "minutes with", num_divergent, "divergences: mean log posterior", round(mean(lp), 1), "\n")
} else {    
    mrhevo_model <- mrhevo(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info)
    mrhevo_mcmc <- setup_sampler(mrhevo_model, dense_mass=dense_mass,
                                 target_accept_prob=target_accept_prob, num_warmup=num_warmup)

    mrhevograph_model <- mrhevoforgraph(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info)

    render_graph(mrhevograph_model, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info, filename="mrhevograph.pdf")

    render_graph(mrhevograph_model, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info, filename="mrhevograph.dot")

    render_graph(mrhevograph_model, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info, filename="mrhevograph.svg")

    start_time <- Sys.time()
    mrhevo_mcmc <- run_sampler(mrhevo_mcmc, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info)
    elapsed.time <- Sys.time() - start_time
    mrhevo_idata <- get_idata(mrhevo_mcmc)
    mrhevo_diagnostics <- get_diagnostics(mrhevo_idata)
    mrhevo_diagnostics.names <- sapply(mrhevo_diagnostics, function(x) np$asarray(x[[1]]))
    ## names are diverging, lp
    num_divergent <- length(which(np$asarray(mrhevo_diagnostics[[1]][[2]])))
    lp <- np$asarray(mrhevo_diagnostics[[2]][[2]]) # dimension chains x numdraws
    cat("Fitting with target acceptance rate", target_accept_prob, "took", elapsed.time, "minutes with", num_divergent, "divergences: mean log posterior", round(mean(lp), 1), "\n")
}

## save trace plots and pairs plot as pdfs 
save_plots(mrhevo_idata)

mrhevo_samples <- get_mcmc_samples(mrhevo_mcmc)
mrhevo_samples.names <- sapply(mrhevo_samples, function(x) np$asarray(x[[1]]))
mrhevo_samples.list <- sapply(mrhevo_samples, function(x) np$asarray(x[[2]]))
names(mrhevo_samples.list) <- mrhevo_samples.names

## table of medians and 80% credible intervals
mrhevo.summaries.list <- lapply(seq_along(mrhevo_samples.list),
                                function(i) {
                                    x <- mrhevo_samples.list[[i]]
                                    param = names(mrhevo_samples.list)[i]
                                    if(length(dim(x)) > 1)
                                        index <- 1:ncol(x) else index <- 1
                                    return(data.table(param, index, summary.param(x)))
                                })
mrhevo.summaries.dt <- rbindlist(mrhevo.summaries.list)

