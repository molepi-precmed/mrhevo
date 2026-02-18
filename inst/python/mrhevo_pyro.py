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

def mrhevo(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info):
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    nu_local = 1
    nu_global = 1
    J = alpha_hat.shape[0]
    mu_alpha = npyr.sample("mu_alpha", dist.Normal(0, 1.0), rng_key=rng_key)
    sigma_alpha = npyr.sample("sigma_alpha", dist.HalfNormal(1.0), rng_key=rng_key)
    aux1_global =  npyr.sample("aux1_global", dist.HalfNormal(1.0), rng_key=rng_key)
    aux2_global = npyr.sample("aux2_global", dist.InverseGamma(0.5 * nu_global, 0.5 * nu_global), rng_key=rng_key) 
    caux = npyr.sample("caux", dist.InverseGamma(0.5 * slab_df, 0.5 * slab_df), rng_key=rng_key) 
    τ = npyr.deterministic("tau", aux1_global * jnp.sqrt(aux2_global) * scale_global) 
    eta = npyr.deterministic("eta", slab_scale * jnp.sqrt(caux)) 
    θ = npyr.sample("theta", dist.Normal(0, priorsd_theta), rng_key=rng_key)
    with npyr.plate("J instruments", J):
        alpha = npyr.sample("alpha", dist.Normal(mu_alpha, sigma_alpha), rng_key=rng_key)
        z = npyr.sample("z", dist.Normal(0, 1.0), rng_key=rng_key)
        aux1_local = npyr.sample("aux1_local", dist.HalfNormal(1.0), rng_key=rng_key)
        aux2_local = npyr.sample("aux2_local", dist.InverseGamma(0.5 * nu_local, 0.5 * nu_local), rng_key=rng_key) 
        lambda_ = npyr.deterministic("lambda_", jnp.multiply(aux1_local, jnp.sqrt(aux2_local))) 
        lambda_tilde = npyr.deterministic("lambda_tilde",
                                          jnp.sqrt(jnp.divide(eta**2 * jnp.square(lambda_),
                                                              eta**2 + τ**2 * jnp.square(lambda_))))
        beta = npyr.deterministic("beta", jnp.multiply(z, lambda_tilde * τ))
        gamma = npyr.deterministic("gamma", beta + θ * alpha)
        kappa = npyr.deterministic("kappa", jnp.divide(1.0 , 1.0 + jnp.multiply(jnp.square(lambda_) * τ**2, info)))
        alpha_hat = npyr.sample("alpha_hat", dist.Normal(alpha, se_alpha_hat), obs=alpha_hat)
        gamma_hat = npyr.sample("gamma_hat", dist.Normal(gamma, se_gamma_hat), obs=gamma_hat)
        b = npyr.deterministic("b", jnp.divide(1.0, (1.0 + eta**2 * info))) 
    f = npyr.deterministic("f", jnp.sum( 1.0 - kappa ) / J) 
    log_eta = npyr.deterministic("log_eta", jnp.log(eta))
    log_τ = npyr.deterministic("log_tau", jnp.log(τ))
    return mrhevo

def mrhorse(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info):
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    nu_local = 1
    nu_global = 1
    J = alpha_hat.shape[0]
    mu_alpha = npyr.sample("mu_alpha", dist.Normal(0, 1.0), rng_key=rng_key)
    sigma_alpha = npyr.sample("sigma_alpha", dist.HalfNormal(1.0), rng_key=rng_key)
    aux1_global =  npyr.sample("aux1_global", dist.HalfNormal(1.0), rng_key=rng_key)
    aux2_global = npyr.sample("aux2_global", dist.InverseGamma(0.5 * nu_global, 0.5 * nu_global), rng_key=rng_key) 
    τ = npyr.deterministic("tau", aux1_global * jnp.sqrt(aux2_global) * scale_global) 
    θ = npyr.sample("theta", dist.Normal(0, priorsd_theta), rng_key=rng_key)
    with npyr.plate("J instruments", J):
        alpha = npyr.sample("alpha", dist.Normal(mu_alpha, sigma_alpha), rng_key=rng_key)
        z = npyr.sample("z", dist.Normal(0, 1.0), rng_key=rng_key)
        aux1_local = npyr.sample("aux1_local", dist.HalfNormal(1.0), rng_key=rng_key)
        aux2_local = npyr.sample("aux2_local", dist.InverseGamma(0.5 * nu_local, 0.5 * nu_local), rng_key=rng_key) 
        lambda_ = npyr.deterministic("lambda_", jnp.multiply(aux1_local, jnp.sqrt(aux2_local))) 
        lambda_tilde = npyr.deterministic("lambda_tilde", lambda_)
        beta = npyr.deterministic("beta", jnp.multiply(z, lambda_tilde * τ))
        gamma = npyr.deterministic("gamma", beta + θ * alpha)
        kappa = npyr.deterministic("kappa", jnp.divide(1.0 , 1.0 + jnp.multiply(jnp.square(lambda_) * τ**2, info)))
        alpha_hat = npyr.sample("alpha_hat", dist.Normal(alpha, se_alpha_hat), obs=alpha_hat)
        gamma_hat = npyr.sample("gamma_hat", dist.Normal(gamma, se_gamma_hat), obs=gamma_hat)
    f = npyr.deterministic("f", jnp.sum( 1.0 - kappa ) / J) 
    log_τ = npyr.deterministic("log_tau", jnp.log(τ))
    return mrhorse

def mrhevoforgraph(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info):
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    nu_local = 1
    nu_global = 1
    J = alpha_hat.shape[0]
    #μ_α = npyr.sample("<μ<sub>α</sub>>", dist.Normal(0, 1.0), rng_key=rng_key)
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
        κ = npyr.deterministic("κ", jnp.divide(1.0 , 1.0 + jnp.multiply(jnp.square(λ) * τ**2, info)))
        α_hat = npyr.sample("α^", dist.Normal(α, se_alpha_hat), obs=alpha_hat)
        γ_hat = npyr.sample("γ^", dist.Normal(γ, se_gamma_hat), obs=gamma_hat)
    f = npyr.deterministic("f", jnp.sum( 1.0 - κ ) / J) 
    return mrhevoforgraph

def nullhevo(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info):
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    nu_local = 1
    nu_global = 1
    J = alpha_hat.shape[0]
    mu_alpha = npyr.sample("mu_alpha", dist.Normal(0, 1.0), rng_key=rng_key)
    sigma_alpha = npyr.sample("sigma_alpha", dist.HalfNormal(1.0), rng_key=rng_key)
    aux1_global =  npyr.sample("aux1_global", dist.HalfNormal(1.0), rng_key=rng_key)
    aux2_global = npyr.sample("aux2_global", dist.InverseGamma(0.5 * nu_global, 0.5 * nu_global), rng_key=rng_key) 
    caux = npyr.sample("caux", dist.InverseGamma(0.5 * slab_df, 0.5 * slab_df), rng_key=rng_key) 
    τ = npyr.deterministic("tau", aux1_global * jnp.sqrt(aux2_global) * scale_global) 
    eta = npyr.deterministic("eta", slab_scale * jnp.sqrt(caux)) 
    with npyr.plate("J instruments", J):
        alpha = npyr.sample("alpha", dist.Normal(mu_alpha, sigma_alpha), rng_key=rng_key)
        z = npyr.sample("z", dist.Normal(0, 1.0), rng_key=rng_key)
        aux1_local = npyr.sample("aux1_local", dist.HalfNormal(1.0), rng_key=rng_key)
        aux2_local = npyr.sample("aux2_local", dist.InverseGamma(0.5 * nu_local, 0.5 * nu_local), rng_key=rng_key) 
        lambda_ = npyr.deterministic("lambda_", jnp.multiply(aux1_local, jnp.sqrt(aux2_local))) 
        lambda_tilde = npyr.deterministic("lambda_tilde",
                                          jnp.sqrt(jnp.divide(eta**2 * jnp.square(lambda_),
                                                              eta**2 + τ**2 * jnp.square(lambda_))))
        beta = npyr.deterministic("beta", jnp.multiply(z, lambda_tilde * τ))
        gamma = npyr.deterministic("gamma", beta)
        kappa = npyr.deterministic("kappa", jnp.divide(1.0 , 1.0 + jnp.multiply(jnp.square(lambda_) * τ**2, info)))
        alpha_hat = npyr.sample("alpha_hat", dist.Normal(alpha, se_alpha_hat), obs=alpha_hat)
        gamma_hat = npyr.sample("gamma_hat", dist.Normal(gamma, se_gamma_hat), obs=gamma_hat)
        b = npyr.deterministic("b", jnp.divide(1.0, (1.0 + eta**2 * info))) 
    f = npyr.deterministic("f", jnp.sum( 1.0 - kappa ) / J) 
    log_eta = npyr.deterministic("log_eta", jnp.log(eta))
    log_τ = npyr.deterministic("log_tau", jnp.log(τ))
    return nullhevo

def setup_sampler(model, dense_mass=False, target_accept_prob=0.95, num_warmup=500):
    # https://num.pyro.ai/en/latest/mcmc.html
    nuts_kernel = NUTS(model, target_accept_prob=target_accept_prob, dense_mass=dense_mass)
    mcmc = MCMC(nuts_kernel, num_warmup=num_warmup, num_samples=1000, num_chains=4,
                chain_method="parallel", progress_bar=True)
    
    return mcmc

def render_graph(mcmc, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info, filename):
    mrhevograph = npyr.render_model(mcmc, model_args=(alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info), render_distributions=True, filename=filename)

    
def run_sampler(mcmc, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info):
    np.random.seed(42) 
    rng_key = random.PRNGKey(42)
    rng_key, rng_key_ = random.split(rng_key)
    mcmc.run(rng_key, alpha_hat, se_alpha_hat, gamma_hat, se_gamma_hat, slab_scale, slab_df, scale_global, priorsd_theta, info, extra_fields=('potential_energy',))
    return mcmc

def get_idata(mcmc):
    idata = az.from_numpyro(mcmc)
    return(idata)

def get_diagnostics(idata):
    samplestats = idata.sample_stats 
    samplestats_keyvaluepairs = samplestats.items()
    samplestats_list = list(samplestats_keyvaluepairs)
    return samplestats_list

def save_plots(idata):
    plt.rcParams['figure.constrained_layout.use'] = True
    az.plot_trace(idata, var_names=["theta", "tau", "eta", "f", "beta", "kappa"], compact=True)
    plt.savefig("mrhevo_traceplot.pdf")
    az.plot_pair(idata, var_names=["theta", "log_tau", "log_eta", "f"], marginals=True,
                 textsize=32, divergences=True)
    plt.savefig("mrhevo_pairsplot.pdf")

def get_mcmc_samples(mcmc):
    posterior_samples = mcmc.get_samples() 
    posterior_samples = {key: posterior_samples[key] for key in posterior_samples.keys()
                         & {"theta", "tau", "log_tau", "eta", "log_eta", "f", "b", "alpha", "beta", "kappa", "lambda_tilde"}} 
    samples_keyvaluepairs = posterior_samples.items()
    samples_list = list(samples_keyvaluepairs)
    return samples_list

def get_null_mcmc_samples(mcmc):
    posterior_samples = mcmc.get_samples() 
    posterior_samples = {key: posterior_samples[key] for key in posterior_samples.keys()
                         & {"tau", "f", "alpha", "lambda_tilde"}} 
    samples_keyvaluepairs = posterior_samples.items()
    samples_list = list(samples_keyvaluepairs)
    return samples_list

def main():
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

