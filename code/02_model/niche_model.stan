// "niche model": a single-species niche model

data {
  int<lower=0> K; // number of niche dimensions
  int<lower=0> N; // number of data points
  matrix[N,K] m; // macroclimate
  array[N] int<lower=0> nn; // number of (sub)plots
  array[N] int<lower=0> y; // number of presences
}

parameters {
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
  vector[N] p;
  x = m - rep_matrix(mu, N); // climate re-centered on mu
  p = alpha * exp(-.5 * rows_dot_product(x / sigma, x));
  target += binomial_lpmf(y | nn, p);
  
  // priors //
  target += normal_lpdf(to_vector(mu) | 0, 1);
  target += gamma_lpdf(to_vector(tau) | 1.5, 1);
  target += beta_lpdf(alpha | 1.5, 3);
  target += lkj_corr_cholesky_lpdf(Lomega | 10);
}

generated quantities {
  real test = binomial_lpmf(2 | 7, .1);
  real loglik_sum;
  vector[N] log_lik;
  vector[N] pp;
  for(n in 1:N) {
          row_vector[K] m0 = m[n] - mu;
          pp[n] = alpha * exp(-.5 * m0 * inverse(sigma) * m0');
          log_lik[n] = binomial_lpmf(y[n] | nn[n], pp[n]);
  }
  loglik_sum = sum(log_lik);
}
