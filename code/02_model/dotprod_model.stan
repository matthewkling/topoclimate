functions {
  real partial_sum(
    int[] si, // species indices
    int start, 
    int end,
    
    int[] y, // all the occurrences
    int[] nn, // all the sample sizes
    matrix xc, // all the predictor data
    
    matrix[] sigma, // all the sigmas
    vector alpha, // all the alphas
    
    int[] sp1, // index of first record for each species, length S
    int[] spn, // number of records per species, length S
    
    int K
    ) {
      int sp1s;
      real alph;
      matrix[K,K] inv_sigma;
      
      real subtarget = 0.0;
      
      for(s in start:end){
        vector[spn[s]] p;
        alph = alpha[s];
        inv_sigma = inverse(sigma[s]);
        sp1s = sp1[s];
        
        // for(i in 1:spn[s]) p[i] = alph * exp(-.5 * (xc[sp1s + i - 1] * inv_sigma * xc[sp1s + i - 1]'));
        p = alph * exp(-.5 * rows_dot_product(xc[sp1s:(sp1s + spn[s] - 1)] * inv_sigma, xc[sp1s:(sp1s + spn[s] - 1)]'));
        subtarget += binomial_lpmf(y[sp1s:(sp1s + spn[s] - 1)] | nn[sp1s:(sp1s + spn[s] - 1)], p);
      }
      return subtarget;
    }
}

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
  
  
  int<lower=0> si[S]; // species index, i.e. 1:S
  int<lower=0> sp1[S]; // index of first record of each species
  int<lower=0> spn[S]; // number of records for each species
  
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
    start = sp1[s];
    end = start + spn[s] - 1;
    
    target += binomial_lpmf(
      y[start:end] | 
      nn[start:end], 
      alpha[s] * exp(-.5 * rows_dot_product(x[start:end] / sigma[s], x[start:end])));
  }
  // target += reduce_sum(partial_sum, si, 1, y, nn, xc, sigma, alpha, sp1, spn, K);
}
