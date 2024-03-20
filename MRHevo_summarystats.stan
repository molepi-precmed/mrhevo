/* MRHevo_summarystats.stan */
functions {
  // elementwise inverse of a vector
  vector inv_vec(vector x) {
    vector[dims(x)[1]] res;
    for (m in 1:dims(x)[1]) {
      res[m] = 1 / x[m];
    }
    return res;
  }
}

data {
  int <lower=0> J; // number of instruments
  real < lower=0 > scale_global; // scale for the half-t prior on tau
  real < lower =1 > nu_global ; // degrees of freedom for the half-t prior on global scale param tau
  real < lower =1 > nu_local ; // degrees of freedom for the half-t priors on local scale parameters
  real < lower=0 > slab_scale; // slab scale for regularized horseshoe
  real < lower=0 > slab_df; // slab degrees of freedom for regularized horseshoe
  real < lower=0 > priorsd_theta; // sd of prior on theta
  vector[J] gamma_hat; // estimated coeffs for regression of Y on Z
  vector <lower=0> [J] sd_gamma_hat; // standard error of estimated coeffs gamma_hat
  vector[J] alpha_hat; // estimated coeffs for effects of Z on X 
  vector <lower=0> [J] sd_alpha_hat; // standard error of estimated coeffs alpha_hat
  vector <lower=0> [J] info; // inverse variance of ratio estimator
}

parameters {
  real < lower=0 > aux1_global;
  real < lower=0 > aux2_global;
  real < lower=0 > caux ;
  vector[J] z; // pleiotropic effects, before global and local scaling
  real theta; // causal effect of X on Y
  vector[J] alpha; // effect of each instrument Z on X, distributed as N(alpha_hat, sd.alpha_hat)
  vector < lower=0 > [J] aux1_local;
  vector < lower=0 > [J] aux2_local;
}

transformed parameters {
  vector < lower =0 > [J] lambda ; // local shrinkage parameters
  real<lower=0> tau;
  vector <lower=0> [J] lambda_tilde;// 'truncated' local shrinkage parameters 
  real < lower =0 > c; // slab scale
  vector[J] gamma; // effect of each instrument Z on Y, distributed as N(gamma_hat, sd.gamma. hat)
  vector[J] beta; // pleiotropic effects after scaling

  lambda   = aux1_local .* sqrt(aux2_local); // local scale parameters
  //  t distribution: standard Gaussian scaled by inverse gamma with parameters 0.5 * \nu, 0.5 * \nu.  This is scaled again by scale_global 
  tau = aux1_global * sqrt(aux2_global) * scale_global;  // global scale parameter for pleiotropic effects
  c = slab_scale * sqrt(caux); // manuscript uses s_slab for slab_scale, \eta for c^2
  lambda_tilde = sqrt( c^2 * square(lambda) ./ (c^2 + tau^2 * square(lambda) ));
  beta = z .* lambda_tilde * tau; // scaled vector of pleiotropic effects
  gamma = beta + theta * alpha;
}

model {
  alpha_hat ~ normal(alpha, sd_alpha_hat);
  gamma_hat ~ normal(gamma, sd_gamma_hat);
  theta ~ normal(0, priorsd_theta);
  z ~ std_normal();

  // half Student-t priors for local and global scale parameters (nu = 1 corresponds to Cauchy)
  aux1_local ~ std_normal(); // aux1_local
  aux2_local ~ inv_gamma(0.5 * nu_local, 0.5 * nu_local); 
  aux1_global ~ normal(0, 1);  // aux1_global
  aux2_global ~ inv_gamma(0.5 * nu_global, 0.5 * nu_global); 
  // regularized horseshoe has an extra parameter caux
  caux ~ inv_gamma(0.5 * slab_df , 0.5 * slab_df ); // this will be rescaled by slab_scale 
}

generated quantities {
  vector<lower=0>[J] b; // lower limit of each regularized shrinkage coefficient
  vector<lower=0>[J] kappa; // unregularized shrinkage factor for each coefficient
  real<lower=0> f; // effective fraction of nonzero coefficients
  real log_c; // log of regularization parameter
  real log_tau; // log of global scale parameter
  b = 1.0 ./ (1.0 + c^2 * info); 
  kappa = 1 ./ (1.0 + (square(lambda) * tau^2) .* info);
  f = sum( 1.0 - kappa ) / J; 
  log_c = log(c);
  log_tau = log(tau);
}
