data {
  int<lower=0> K; // number of niche dimensions
  int<lower=0> N; // number of data points
  int<lower=1> D; // number of deltas per climate var
  int<lower=0> S; // number of species
  array[N] int<lower=0> ss; // species index
  array[N] int<lower=0> nn; // number of (sub)plots
  array[N] int<lower=0> y; // number of presences
  matrix[N,D] z; // microclimate predictors
  int<lower=1> Nadj; // number of adjacencies (per topo*macro var)
  array[Nadj] int<lower=0> adj2; // adjacency indices of spline bases
  array[Nadj] int<lower=0> adj1; // adjacency indices of spline bases
  array[Nadj] int<lower=0> adj0; // adjacency indices of spline bases
  real lambda; // smoothing parameter
  matrix[N,K] m; // macroclimate
  
  array[S] int<lower=0> s0; // index of first record for each species
  array[S] int<lower=0> s1; // index of last record for each species
  
  matrix[S,K] m_span; // max - min climate at presences
  matrix[S,K] m_min; // min ditto
}

parameters {
  matrix[D,K] delta; // topoclimate regression parameters
  matrix<lower=0,upper=1>[S,K] mu_unit; // niche mean, on 0-1 scale
  matrix<lower=0>[S,K] tau; // niche scale
  array[S] cholesky_factor_corr[K] Lomega;
  vector<lower=0,upper=1>[S] alpha; // niche height
}

transformed parameters {
  array[S] cov_matrix[K] sigma; // niche covariance matrix
  matrix[S,K] mu; // niche mean, rescaled
  matrix<lower=0>[Nadj,K] diffs; // differences between adjacent spline bases
  for(s in 1:S) sigma[s] = quad_form_diag(tcrossprod(Lomega[s]), tau[s]);
  mu = mu_unit .* m_span + m_min;
  diffs = delta[adj0,] - (2 * delta[adj1,]) + delta[adj2,];
  diffs = diffs .* diffs;
}

model {
  matrix[N,K] x;
  int start;
  int end;
  
  target += normal_lpdf(to_vector(delta) | 0, .1);
  target += beta_lpdf(to_vector(mu_unit) | 1.1, 1.1);
  target += gamma_lpdf(to_vector(tau) | 1.5, 1);
  target += beta_lpdf(alpha | 1.5, 3);
  for(s in 1:S) target += lkj_corr_cholesky_lpdf(Lomega[s] | 10);
  target += exponential_lpdf(to_vector(diffs) | pow(10, lambda));
  
  x = m + z * delta - mu[ss]; // topoclimate re-centered on mu
  
  for(s in 1:S){
    start = s0[s];
    end = s1[s];
    
    target += binomial_lpmf(
      y[start:end] | 
      nn[start:end], 
      alpha[s] * exp(-.5 * rows_dot_product(x[start:end] / sigma[s], x[start:end])));
  }
}
