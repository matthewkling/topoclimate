// "minimal model": a single-species coupled topoclimate-niche model, with invariant topoclimate effects

data {
  int<lower=0> K; // number of niche dimensions
  int<lower=0> N; // number of data points
  int<lower=1> D; // number of deltas
  matrix[N,D] z; // microclimate predictors
  matrix[N,K] m; // macroclimate
  int<lower=0> nn[N]; // number of (sub)plots
  int<lower=0> y[N]; // number of presences
}

parameters {
  matrix[D,K] delta; // topoclimate regression parameters
  row_vector[K] mu; // niche mean
  vector<lower=0>[K] tau; // niche scale
  cholesky_factor_corr[K] Lomega;
  real<lower=0,upper=1> alpha; // niche height
}

transformed parameters {
  cov_matrix[K] sigma = quad_form_diag(tcrossprod(Lomega), tau); // niche covariance matrix
}

model {
        
  // likelihood //
  matrix[N,K] x;
  x = m + z * delta - rep_matrix(mu, N); // topoclimate re-centered on mu
  target += binomial_lpmf(y | nn, alpha * exp(-.5 * rows_dot_product(x / sigma, x)));
  
  // priors //
  target += normal_lpdf(to_vector(delta) | 0, .1);
  target += normal_lpdf(to_vector(mu) | 0, 1);
  target += gamma_lpdf(to_vector(tau) | 1.5, 1);
  target += beta_lpdf(alpha | 1.5, 3);
  target += lkj_corr_cholesky_lpdf(Lomega | 10);
}
