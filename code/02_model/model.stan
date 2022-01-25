data {
  int<lower=0> K; // number of niche dimensions
  int<lower=0> N; // number of data points
  int<lower=1> D; // number of deltas
  int<lower=0> S; // number of species
  int<lower=0> ss[N]; // species index
  int<lower=0> nn[N]; // number of trials for each obs
  int<lower=0> y[N]; // species presence/absence
  matrix[N,D] z; // microclimate predictors
  matrix[N,K] m; // macroclimate
  
  int<lower=0> s0[S]; // index of first record for each species
  int<lower=0> s1[S]; // index of last record for each species
}

parameters {
  matrix[D,K] delta; // topoclimate regression parameters
  matrix[S,K] mu; // niche mean
  matrix<lower=0>[S,K] tau; // niche scale
  cholesky_factor_corr[K] Lomega[S];
  vector<lower=0,upper=1>[S] alpha; // niche height
}

transformed parameters {
  cov_matrix[K] sigma[S]; // niche covariance matrix
  for(s in 1:S) sigma[s] = quad_form_diag(tcrossprod(Lomega[s]), tau[s]);
}

model {
  matrix[N,K] x;
  int start;
  int end;
  
  target += normal_lpdf(to_vector(delta) | 0, .1);
  target += normal_lpdf(to_vector(mu) | 0, 1);
  target += normal_lpdf(to_vector(tau) | 0, 10);
  target += beta_lpdf(alpha | 1.5, 3);
  for(s in 1:S) target += lkj_corr_cholesky_lpdf(Lomega[s] | 3);
  
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
