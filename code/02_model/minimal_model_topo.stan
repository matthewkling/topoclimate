// "minimal model": a single-species coupled topoclimate-niche model, with invariant topoclimate effects

data {
        int<lower=0> K; // number of niche dimensions
        int<lower=0> N; // number of data points
        int<lower=1> D; // number of deltas
        matrix[N,D] z; // microclimate predictors
        matrix[N,K] m; // macroclimate
        array[N] int<lower=0> nn; // number of (sub)plots
        array[N] int<lower=0> y; // number of presences
        int apply_beta;
        int apply_delta;
        //int apply_epsilon;
}

parameters {
        matrix[D,K] delta; // topoclimate regression parameters
        row_vector[K] mu; // niche mean
        vector<lower=0>[K] tau; // niche scale
        cholesky_factor_corr[K] Lomega;
        real alpha; // niche intercept
        vector[D] beta; // topo niche params
        //matrix[N,K] epsilon;
        //real<lower=0> epsilon_scale; // hyperparameter
}

transformed parameters {
        cov_matrix[K] sigma = quad_form_diag(tcrossprod(Lomega), tau); // niche covariance matrix
}

model {
        
        // likelihood //
        matrix[N,K] x;
        vector[N] eta;
        x = m + z * delta * apply_delta - rep_matrix(mu, N);// + epsilon; // topoclimate re-centered on mu
        // target += normal_lpdf(to_vector(epsilon) | 0, epsilon_scale);
        // target += gamma_lpdf(epsilon_scale | 1, 10 + (1 - apply_epsilon) * 1000); // force scale to nearly 0 when apply_epsilon is false (actual 0 breaks optimization )
        eta = alpha - rows_dot_product(x / sigma, x) + z * beta * apply_beta;
        target += binomial_logit_lpmf(y | nn, eta);
        
        
        // priors //
        target += normal_lpdf(to_vector(delta) | 0, .1);
        target += normal_lpdf(to_vector(mu) | 0, 1);
        target += gamma_lpdf(to_vector(tau) | 1.5, 1);
        target += normal_lpdf(alpha | 0, 1);
        target += normal_lpdf(beta | 0, 1);
        target += lkj_corr_cholesky_lpdf(Lomega | 10);
}
