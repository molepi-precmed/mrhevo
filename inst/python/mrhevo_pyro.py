# source /home/pmckeigue/.virtualenvs/r-reticulate/bin/activate
# import mrhevo_pyro to run in a session

import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"]="0.20"

import numpy as np
import jax.numpy as jnp
from jax import random

import numpyro as npyr
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS
#npyr.set_platform("cpu")
npyr.set_host_device_count(4)
npyr.enable_x64()

import arviz as az
import matplotlib.pyplot as plt
import graphviz as gv

def mrhevo(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, hierarchical_alpha=1):
    """
    Regularized horseshoe prior for MR with optional hierarchical prior on alpha.
    
    Parameters:
    -----------
    hierarchical_alpha : int
        1 = use hierarchical prior on alpha (default)
        0 = use alpha_hat directly (no shrinkage on alpha)
    """
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    nu_local = 1
    nu_global = 1
    J = alpha_hat.shape[0]
    
    # Global shrinkage parameters
    aux1_global = npyr.sample("aux1_global", dist.HalfNormal(1.0), rng_key=rng_key)
    aux2_global = npyr.sample("aux2_global", dist.InverseGamma(0.5 * nu_global, 0.5 * nu_global), rng_key=rng_key) 
    caux = npyr.sample("caux", dist.InverseGamma(0.5 * slab_df, 0.5 * slab_df), rng_key=rng_key) 
    τ = npyr.deterministic("tau", aux1_global * jnp.sqrt(aux2_global) * scale_global) 
    eta = npyr.deterministic("eta", slab_scale * jnp.sqrt(caux)) 
    θ = npyr.sample("theta", dist.Normal(0, priorsd_theta), rng_key=rng_key)
    
    # Hierarchical prior on alpha
    if hierarchical_alpha == 1:
        mu_alpha = npyr.sample("mu_alpha", dist.Normal(0, 1.0), rng_key=rng_key)
        sigma_alpha = npyr.sample("sigma_alpha", dist.HalfNormal(1.0), rng_key=rng_key)
    
    with npyr.plate("J instruments", J):
        # Alpha: hierarchical or direct
        if hierarchical_alpha == 1:
            alpha = npyr.sample("alpha", dist.Normal(mu_alpha, sigma_alpha), rng_key=rng_key)
        else:
            # Non-hierarchical: use alpha_hat as observed
            alpha = npyr.sample("alpha", dist.Normal(0, 1.0), rng_key=rng_key)
        
        z = npyr.sample("z", dist.Normal(0, 1.0), rng_key=rng_key)
        aux1_local = npyr.sample("aux1_local", dist.HalfNormal(1.0), rng_key=rng_key)
        aux2_local = npyr.sample("aux2_local", dist.InverseGamma(0.5 * nu_local, 0.5 * nu_local), rng_key=rng_key) 
        lambda_ = npyr.deterministic("lambda_", jnp.multiply(aux1_local, jnp.sqrt(aux2_local))) 
        lambda_tilde = npyr.deterministic("lambda_tilde",
                                          jnp.sqrt(jnp.divide(eta**2 * jnp.square(lambda_),
                                                              eta**2 + τ**2 * jnp.square(lambda_))))
        beta = npyr.deterministic("beta", jnp.multiply(z, lambda_tilde * τ))
        gamma = npyr.deterministic("gamma", beta + θ * alpha)
        
        # Kappa calculation: match Stan model (no info term)
        kappa = npyr.deterministic("kappa", jnp.divide(1.0 , 1.0 + jnp.square(lambda_tilde)))
        
        # Observe data
        if hierarchical_alpha == 1:
            alpha_hat_obs = npyr.sample("alpha_hat", dist.Normal(alpha, se_alpha_hat), obs=alpha_hat)
        else:
            # In non-hierarchical, treat alpha as directly observed
            alpha_hat_obs = npyr.sample("alpha_hat", dist.Normal(alpha, se_alpha_hat), obs=alpha_hat)
        
        gamma_hat = npyr.sample("gamma_hat", dist.Normal(gamma, se_gamma_hat), obs=gamma_hat)
    
    f = npyr.deterministic("f", jnp.sum( 1.0 - kappa ) / J) 
    log_eta = npyr.deterministic("log_eta", jnp.log(eta))
    log_τ = npyr.deterministic("log_tau", jnp.log(τ))
    return mrhevo


def mrhorse(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, hierarchical_alpha=1):
    """Unregularized horseshoe (no slab)"""
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    nu_local = 1
    nu_global = 1
    J = alpha_hat.shape[0]
    
    aux1_global = npyr.sample("aux1_global", dist.HalfNormal(1.0), rng_key=rng_key)
    aux2_global = npyr.sample("aux2_global", dist.InverseGamma(0.5 * nu_global, 0.5 * nu_global), rng_key=rng_key) 
    τ = npyr.deterministic("tau", aux1_global * jnp.sqrt(aux2_global) * scale_global) 
    θ = npyr.sample("theta", dist.Normal(0, priorsd_theta), rng_key=rng_key)
    
    if hierarchical_alpha == 1:
        mu_alpha = npyr.sample("mu_alpha", dist.Normal(0, 1.0), rng_key=rng_key)
        sigma_alpha = npyr.sample("sigma_alpha", dist.HalfNormal(1.0), rng_key=rng_key)
    
    with npyr.plate("J instruments", J):
        if hierarchical_alpha == 1:
            alpha = npyr.sample("alpha", dist.Normal(mu_alpha, sigma_alpha), rng_key=rng_key)
        else:
            alpha = npyr.sample("alpha", dist.Normal(0, 1.0), rng_key=rng_key)
        z = npyr.sample("z", dist.Normal(0, 1.0), rng_key=rng_key)
        aux1_local = npyr.sample("aux1_local", dist.HalfNormal(1.0), rng_key=rng_key)
        aux2_local = npyr.sample("aux2_local", dist.InverseGamma(0.5 * nu_local, 0.5 * nu_local), rng_key=rng_key) 
        lambda_ = npyr.deterministic("lambda_", jnp.multiply(aux1_local, jnp.sqrt(aux2_local))) 
        lambda_tilde = npyr.deterministic("lambda_tilde", lambda_)
        beta = npyr.deterministic("beta", jnp.multiply(z, lambda_tilde * τ))
        gamma = npyr.deterministic("gamma", beta + θ * alpha)
        kappa = npyr.deterministic("kappa", jnp.divide(1.0 , 1.0 + jnp.square(lambda_tilde)))
        alpha_hat = npyr.sample("alpha_hat", dist.Normal(alpha, se_alpha_hat), obs=alpha_hat)
        gamma_hat = npyr.sample("gamma_hat", dist.Normal(gamma, se_gamma_hat), obs=gamma_hat)
    f = npyr.deterministic("f", jnp.sum( 1.0 - kappa ) / J) 
    log_τ = npyr.deterministic("log_tau", jnp.log(τ))
    return mrhorse


def mrhevoforgraph(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info):
    """Model for graph rendering"""
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    nu_local = 1
    nu_global = 1
    J = alpha_hat.shape[0]
    μ_α = npyr.sample("μ_α", dist.Normal(0, 1.0), rng_key=rng_key)
    σ_α = npyr.sample("σ_α", dist.HalfNormal(1.0), rng_key=rng_key)
    η_aux = npyr.sample("η_aux", dist.InverseGamma(0.5 * slab_df, 0.5 * slab_df), rng_key=rng_key) 
    τ = npyr.sample("τ", dist.HalfCauchy(1.0), rng_key=rng_key)
    η = npyr.deterministic("η", slab_scale * jnp.sqrt(η_aux)) 
    θ = npyr.sample("θ", dist.Normal(0, priorsd_theta), rng_key=rng_key)
    with npyr.plate("J instruments", J):
        α = npyr.sample("α", dist.Normal(μ_α, σ_α), rng_key=rng_key)
        λ = npyr.sample("λ", dist.HalfCauchy(1.0), rng_key=rng_key)
        λ_tilde = npyr.deterministic("λ~", 
                                     jnp.sqrt(jnp.divide(η**2 * jnp.square(λ),
                                                         η**2 + τ**2 * jnp.square(λ))))
        β = npyr.sample("β", dist.Normal(0, λ_tilde * τ), rng_key=rng_key)
        γ = npyr.deterministic("γ", β + θ * α)
        κ = npyr.deterministic("κ", jnp.divide(1.0 , 1.0 + jnp.square(λ_tilde)))
        α_hat = npyr.sample("α^", dist.Normal(α, se_alpha_hat), obs=alpha_hat)
        γ_hat = npyr.sample("γ^", dist.Normal(γ, se_gamma_hat), obs=gamma_hat)
    f = npyr.deterministic("f", jnp.sum( 1.0 - κ ) / J) 
    return mrhevoforgraph


def nullhevo(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, hierarchical_alpha=1):
    """Null model (no causal effect, theta = 0)"""
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    nu_local = 1
    nu_global = 1
    J = alpha_hat.shape[0]
    
    if hierarchical_alpha == 1:
        mu_alpha = npyr.sample("mu_alpha", dist.Normal(0, 1.0), rng_key=rng_key)
        sigma_alpha = npyr.sample("sigma_alpha", dist.HalfNormal(1.0), rng_key=rng_key)
    
    aux1_global = npyr.sample("aux1_global", dist.HalfNormal(1.0), rng_key=rng_key)
    aux2_global = npyr.sample("aux2_global", dist.InverseGamma(0.5 * nu_global, 0.5 * nu_global), rng_key=rng_key) 
    caux = npyr.sample("caux", dist.InverseGamma(0.5 * slab_df, 0.5 * slab_df), rng_key=rng_key) 
    τ = npyr.deterministic("tau", aux1_global * jnp.sqrt(aux2_global) * scale_global) 
    eta = npyr.deterministic("eta", slab_scale * jnp.sqrt(caux)) 
    θ = 0.0  # Null model: no causal effect
    
    with npyr.plate("J instruments", J):
        if hierarchical_alpha == 1:
            alpha = npyr.sample("alpha", dist.Normal(mu_alpha, sigma_alpha), rng_key=rng_key)
        else:
            alpha = npyr.sample("alpha", dist.Normal(0, 1.0), rng_key=rng_key)
        z = npyr.sample("z", dist.Normal(0, 1.0), rng_key=rng_key)
        aux1_local = npyr.sample("aux1_local", dist.HalfNormal(1.0), rng_key=rng_key)
        aux2_local = npyr.sample("aux2_local", dist.InverseGamma(0.5 * nu_local, 0.5 * nu_local), rng_key=rng_key) 
        lambda_ = npyr.deterministic("lambda_", jnp.multiply(aux1_local, jnp.sqrt(aux2_local))) 
        lambda_tilde = npyr.deterministic("lambda_tilde",
                                          jnp.sqrt(jnp.divide(eta**2 * jnp.square(lambda_),
                                                              eta**2 + τ**2 * jnp.square(lambda_))))
        beta = npyr.deterministic("beta", jnp.multiply(z, lambda_tilde * τ))
        gamma = npyr.deterministic("gamma", beta)  # No causal effect
        kappa = npyr.deterministic("kappa", jnp.divide(1.0 , 1.0 + jnp.square(lambda_tilde)))
        alpha_hat = npyr.sample("alpha_hat", dist.Normal(alpha, se_alpha_hat), obs=alpha_hat)
        gamma_hat = npyr.sample("gamma_hat", dist.Normal(gamma, se_gamma_hat), obs=gamma_hat)
    f = npyr.deterministic("f", jnp.sum( 1.0 - kappa ) / J) 
    log_eta = npyr.deterministic("log_eta", jnp.log(eta))
    log_τ = npyr.deterministic("log_tau", jnp.log(τ))
    return nullhevo


# Module-level cache: reuse compiled MCMC objects across calls so JAX does not
# recompile the model on every invocation.
_mcmc_cache = {}

def get_samples_numpy(mcmc, sample_names):
    """Extract named posterior samples as a dict of plain NumPy arrays.

    Converts all requested parameters from JAX arrays to NumPy in one Python-
    level pass, avoiding the overhead of multiple R→Python round-trips.

    Parameters
    ----------
    mcmc : MCMC
        Completed MCMC object.
    sample_names : list of str
        Parameter names to extract.

    Returns
    -------
    dict
        Mapping name → numpy.ndarray for each name found in get_samples().
    """
    samples = mcmc.get_samples()
    return {name: np.asarray(samples[name]) for name in sample_names if name in samples}


def get_cached_mcmc(target_accept_prob=0.95, num_warmup=500, num_samples=1000, num_chains=4):
    """Return a cached MCMC object for the mrhevo model.

    On the first call with a given set of sampler hyperparameters the MCMC
    object is created and JAX compiles the model.  Subsequent calls with the
    same hyperparameters return the existing object, skipping recompilation.

    Parameters
    ----------
    target_accept_prob : float
        Target acceptance probability for NUTS (default 0.95).
    num_warmup : int
        Number of warmup steps (default 500).
    num_samples : int
        Number of posterior samples per chain (default 1000).
    num_chains : int
        Number of parallel chains (default 4).

    Returns
    -------
    MCMC
        Configured (and possibly already compiled) MCMC object.
    """
    key = (target_accept_prob, num_warmup, num_samples, num_chains)
    if key not in _mcmc_cache:
        nuts = NUTS(mrhevo, target_accept_prob=target_accept_prob)
        _mcmc_cache[key] = MCMC(
            nuts,
            num_warmup=num_warmup,
            num_samples=num_samples,
            num_chains=num_chains,
            chain_method="parallel",
            progress_bar=True,
        )
    return _mcmc_cache[key]


def setup_sampler(model, dense_mass=False, target_accept_prob=0.95, num_warmup=500):
    """Set up a NumPyro NUTS sampler.

    Parameters
    ----------
    model : callable
        NumPyro probabilistic model function.
    dense_mass : bool, optional
        Whether to use a dense mass matrix (default False).
    target_accept_prob : float, optional
        Target acceptance probability for the NUTS kernel (default 0.95).
    num_warmup : int, optional
        Number of warmup (burn-in) iterations (default 500).

    Returns
    -------
    MCMC
        Configured MCMC object ready to run.
    """
    # https://num.pyro.ai/en/latest/mcmc.html
    nuts_kernel = NUTS(model, target_accept_prob=target_accept_prob, dense_mass=dense_mass)
    mcmc = MCMC(nuts_kernel, num_warmup=num_warmup, num_samples=1000, num_chains=4,
                chain_method="parallel", progress_bar=True)
    
    return mcmc


def render_graph(mcmc, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info, filename):
    """Render the probabilistic graphical model to a file.

    Parameters
    ----------
    mcmc : callable
        NumPyro model function (e.g. mrhevo).
    alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat : array-like
        Summary statistics passed to the model.
    slab_scale : float
        Slab scale for the regularized horseshoe prior.
    slab_df : float
        Degrees of freedom for the slab component.
    scale_global : float
        Scale for the global shrinkage parameter tau.
    priorsd_theta : float
        Standard deviation of the prior on the causal effect theta.
    info : float
        Mean Fisher information for the IV estimates.
    filename : str
        Output filename for the rendered graph.
    """
    mrhevograph = npyr.render_model(mcmc, model_args=(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info), render_distributions=True, filename=filename)

    
def run_sampler(mcmc, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info):
    """Run the NUTS sampler and return the MCMC object.

    Parameters
    ----------
    mcmc : MCMC
        Configured MCMC object (from setup_sampler).
    alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat : array-like
        Summary statistics passed as model arguments.
    slab_scale : float
        Slab scale for the regularized horseshoe prior.
    slab_df : float
        Degrees of freedom for the slab component.
    scale_global : float
        Scale for the global shrinkage parameter tau.
    priorsd_theta : float
        Standard deviation of the prior on the causal effect theta.
    info : float
        Mean Fisher information for the IV estimates.

    Returns
    -------
    MCMC
        MCMC object with posterior samples.
    """
    np.random.seed(42)
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    mcmc.run(rng_key, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info, extra_fields=('potential_energy',))
    return mcmc

def get_idata(mcmc):
    """Convert MCMC output to an ArviZ InferenceData object.

    Parameters
    ----------
    mcmc : MCMC
        Completed MCMC object.

    Returns
    -------
    arviz.InferenceData
        ArviZ inference data object with posterior samples and diagnostics.
    """
    idata = az.from_numpyro(mcmc)
    return(idata)

def get_diagnostics(idata):
    """Extract sampler diagnostics from an ArviZ InferenceData object.

    Parameters
    ----------
    idata : arviz.InferenceData
        ArviZ inference data object.

    Returns
    -------
    list
        List of (key, value) pairs from the sample statistics group.
    """
    samplestats = idata.sample_stats
    samplestats_keyvaluepairs = samplestats.items()
    samplestats_list = list(samplestats_keyvaluepairs)
    return samplestats_list

def save_plots(idata):
    """Save trace and pairs plots to PDF files in the current directory.

    Parameters
    ----------
    idata : arviz.InferenceData
        ArviZ inference data object with posterior samples.
    """
    plt.rcParams['figure.constrained_layout.use'] = True
    az.plot_trace(idata, var_names=["theta", "tau", "eta", "f", "beta", "kappa"], compact=True)
    plt.savefig("mrhevo_traceplot.pdf")
    az.plot_pair(idata, var_names=["theta", "log_tau", "log_eta", "f"], marginals=True,
                 textsize=32, divergences=True)
    plt.savefig("mrhevo_pairsplot.pdf")

def get_mcmc_samples(mcmc):
    """Extract a curated set of posterior samples from a completed MCMC run.

    Returns samples for the parameters theta, tau, log_tau, eta, log_eta, f,
    b, alpha, beta, kappa, and lambda_tilde.

    Parameters
    ----------
    mcmc : MCMC
        Completed MCMC object.

    Returns
    -------
    list
        List of (name, array) pairs for the requested parameters.
    """
    posterior_samples = mcmc.get_samples()
    posterior_samples = {key: posterior_samples[key] for key in posterior_samples.keys()
                         & {"theta", "tau", "log_tau", "eta", "log_eta", "f", "b", "alpha", "beta", "kappa", "lambda_tilde"}} 
    samples_keyvaluepairs = posterior_samples.items()
    samples_list = list(samples_keyvaluepairs)
    return samples_list

def get_null_mcmc_samples(mcmc):
    """Extract posterior samples for the null model (theta fixed to zero).

    Returns samples for tau, f, alpha, and lambda_tilde only.

    Parameters
    ----------
    mcmc : MCMC
        Completed MCMC object from the null model.

    Returns
    -------
    list
        List of (name, array) pairs for the requested parameters.
    """
    posterior_samples = mcmc.get_samples()
    posterior_samples = {key: posterior_samples[key] for key in posterior_samples.keys()
                         & {"tau", "f", "alpha", "lambda_tilde"}} 
    samples_keyvaluepairs = posterior_samples.items()
    samples_list = list(samples_keyvaluepairs)
    return samples_list

def main():
    """Run a complete MR-Hevo analysis from an R data file.

    Reads model inputs from ``mrhevo_example.RData`` in the current directory,
    fits the regularized horseshoe model with NumPyro, prints diagnostics, and
    saves trace and pairs plots.
    """
    data = pyr.read_rda('mrhevo_example.RData')
    alpha_hat = data["alpha_hat"]
    se_alpha_hat = data["se.alpha_hat"]
    gamma_hat = data["gamma_hat"]
    se_gamma_hat = data["se.gamma_hat"]
    slab_scale = data["slab_scale"]
    slab_df = data["slab_df"]
    scale_global = data["tau0"]
    priorsd_theta = data["priorsd_theta"]
    info = data["info"]
    np.random.seed(42) #probably don't need next 2 lines
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    mrhevo_model = mrhevo(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info)
    mrhevo_sampler = setup_sampler(mrhevo_model)
    mrhevo_mcmc = run_sampler(mrhevo_sampler, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info)
    mrhevo_idata = get_idata(mrhevo_mcmc)
    mrhevo_diagnostics = get_diagnostics(mrhevo_idata)
    save_plots(mrhevo_idata)
    mrhevo_mcmc.print_summary()
    pe = mrhevo_mcmc.get_extra_fields()['potential_energy']
    print('Expected log joint density: {:.2f}'.format(np.mean(-pe)))
    mrhevo_samples = get_mcmc_samples(mrhevo_mcmc)
    #from numpyro.contrib.nested_sampling import NestedSampler
    #mrhevo_ns = NestedSampler(mrhevo)
    #mrhevo_ns.run(random.PRNGKey(42), alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info)
    #samples = ns.get_samples(random.PRNGKey(3), num_samples=1000)


#if __name__ == '__main__':
#    main()

