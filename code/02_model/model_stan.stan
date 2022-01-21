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
  vector[N] p; // occurrence probability
  matrix[N,K] xc; // topoclimate re-centered on mu
  
  to_vector(delta) ~ normal(0, .1); // topoclimate priors
  to_vector(mu) ~ normal(0, 1);
  to_vector(tau) ~ normal(0, 10);
  alpha ~ beta(1.5, 3);
  
  for(s in 1:S){
    Lomega[s] ~ lkj_corr_cholesky(3);
  }
  
  xc = m + z * delta - mu[ss]; 
  for(i in 1:N) p[i] = alpha[ss[i]] * exp(-.5 * (xc[i] / sigma[ss[i]] * xc[i]'));
  y ~ binomial(nn, p);
}
