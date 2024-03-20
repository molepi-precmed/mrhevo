/* MRHevo_logistic.stan */
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
  int <lower=0, upper=1> logistic;
  int <lower=0> N;
  int <lower=0> U; // number of unpenalized covariates
  int <lower=0> J; // number of instruments
  real < lower=0 > scale_intercept_y; // prior sd for Y intercept
  real < lower=0 > scale_beta_u; // prior sd for effects of unpenalized covariates
  real < lower=0 > scale_global; // scale for the half-t prior on tau
  real < lower =1 > nu_global ; // degrees of freedom for the half-t prior on global scale param tau
  real < lower =1 > nu_local ; // degrees of freedom for the half-t priors on local scale parameters
  real < lower=0 > slab_scale; // slab scale for regularized horseshoe
  real < lower=0 > slab_df; // slab degrees of freedom for regularized horseshoe
  real < lower=0 > priorsd_theta; // sd of prior on theta
  matrix[N, U] X_u; // unpenalized covariates
  matrix[N, J] Z; // instruments
  //  vector[N] Y; // uncomment for linear regresssion
  int < lower=0, upper=1 > Y[N]; 
  vector[J] alpha_hat; // estimated coeffs for effects of Z on X 
  vector <lower=0> [J] sd_alpha_hat; // standard error of estimated coeffs alpha_hat
}

transformed data {
  matrix[N, U] Q_ast;
   // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(X_u) * sqrt(N - 1);
 }

parameters {
  //  real <lower=0> sigma_y;
  real < lower=0 > aux1_global;
  real < lower=0 > aux2_global;
  real < lower=0 > caux ;
  real intercept_y; // intercept for Y
  vector[U] beta_u; // coefficients for QR-transformed unpenalized covariates
  vector[J] alpha; // effect of each instrument Z on X, distributed as N(alpha_hat, sd.alpha_hat)
  vector[J] z; // pleiotropic effects, before global and local scaling
  real theta; // causal effect of X on Y
  vector < lower=0 > [J] aux1_local;
  vector < lower=0 > [J] aux2_local;
}

transformed parameters {
  vector[N] Xpred; // predicted value of X given Z, alpha
  vector < lower =0 > [J] lambda ; // local shrinkage parameters
  real<lower=0> tau;
  vector <lower=0> [J] lambda_tilde;// 'truncated' local shrinkage parameters 
  vector[J] beta;
  real < lower =0 > c; // slab scale
  lambda   = aux1_local .* sqrt(aux2_local); // local scale parameters

  //  standard t distribution: standard Gaussian scaled by inverse gamma with parameters 0.5 * \nu, 0.5 * \nu.  This is scaled again by scale_global 
  tau = aux1_global * sqrt(aux2_global) * scale_global;  // global scale parameter for pleiotropic effects
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt( c^2 * square(lambda) ./ (c^2 + tau^2 * square(lambda) ));
  beta = z .* lambda_tilde * tau; // scaled vector of pleiotropic effects
  Xpred = Z * alpha; // N x J matrix postmultiplied by column vector of length J gives column vector of length N 
}

model {
  // theta should be a column vector of length 1
  // Y ~ normal(intercept_y + Z * beta + Xpred * theta, sigma_y); 
  Y ~ bernoulli_logit(intercept_y + Q_ast * beta_u + Z * beta + Xpred * theta);
  alpha ~ normal(alpha_hat, sd_alpha_hat);
  intercept_y ~ normal(0, scale_intercept_y);
  // half Student-t priors for local and global scale parameters (nu = 1 corresponds to horseshoe)
  beta_u ~ normal(0, scale_beta_u); 
  theta ~ normal(0, priorsd_theta);
  z ~ std_normal();

  aux1_local ~ std_normal(); // aux1_local
  aux2_local ~ inv_gamma(0.5 * nu_local, 0.5 * nu_local); 
  aux1_global ~ normal(0, 1);  // aux1_global
  aux2_global ~ inv_gamma(0.5 * nu_global, 0.5 * nu_global); 
  // regularized horseshoe has an extra parameter caux
  caux ~ inv_gamma(0.5 * slab_df , 0.5 * slab_df );
}

generated quantities {
  vector[J] kappa; // shrinkage factors for each coefficient
  real f; // effective fraction of nonzero coefficients
  real log_c; // log of regularization parameter
  real log_tau; // log of global scale parameter
  kappa = inv_vec(1.0 + lambda_tilde .* lambda_tilde);
  f = sum(1 - kappa) / J; 
  log_c = log(c);
  log_tau = log(tau);
}
