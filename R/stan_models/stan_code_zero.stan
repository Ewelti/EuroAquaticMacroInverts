// generated with brms 2.15.0
functions {
  /* zero-inflated poisson log-PDF of a single response 
  * Args: 
    *   y: the response value 
  *   lambda: mean parameter of the poisson distribution
  *   zi: zero-inflation probability
  * Returns:  
    *   a scalar to be added to the log posterior 
  */ 
    real zero_inflated_poisson_lpmf(int y, real lambda, real zi) { 
      if (y == 0) { 
        return log_sum_exp(bernoulli_lpmf(1 | zi), 
                           bernoulli_lpmf(0 | zi) + 
                             poisson_lpmf(0 | lambda)); 
      } else { 
        return bernoulli_lpmf(0 | zi) +  
          poisson_lpmf(y | lambda); 
      } 
    }
  /* zero-inflated poisson log-PDF of a single response 
  * logit parameterization of the zero-inflation part
  * Args: 
    *   y: the response value 
  *   lambda: mean parameter of the poisson distribution
  *   zi: linear predictor for zero-inflation part 
  * Returns:  
    *   a scalar to be added to the log posterior 
  */ 
    real zero_inflated_poisson_logit_lpmf(int y, real lambda, real zi) { 
      if (y == 0) { 
        return log_sum_exp(bernoulli_logit_lpmf(1 | zi), 
                           bernoulli_logit_lpmf(0 | zi) + 
                             poisson_lpmf(0 | lambda)); 
      } else { 
        return bernoulli_logit_lpmf(0 | zi) +  
          poisson_lpmf(y | lambda); 
      } 
    }
  /* zero-inflated poisson log-PDF of a single response
  * log parameterization for the poisson part
  * Args: 
    *   y: the response value 
  *   eta: linear predictor for poisson distribution
  *   zi: zero-inflation probability
  * Returns:  
    *   a scalar to be added to the log posterior 
  */ 
    real zero_inflated_poisson_log_lpmf(int y, real eta, real zi) { 
      if (y == 0) { 
        return log_sum_exp(bernoulli_lpmf(1 | zi), 
                           bernoulli_lpmf(0 | zi) + 
                             poisson_log_lpmf(0 | eta)); 
      } else { 
        return bernoulli_lpmf(0 | zi) +  
          poisson_log_lpmf(y | eta); 
      } 
    }
  /* zero-inflated poisson log-PDF of a single response 
  * log parameterization for the poisson part
  * logit parameterization of the zero-inflation part
  * Args: 
    *   y: the response value 
  *   eta: linear predictor for poisson distribution
  *   zi: linear predictor for zero-inflation part 
  * Returns:  
    *   a scalar to be added to the log posterior 
  */ 
    real zero_inflated_poisson_log_logit_lpmf(int y, real eta, real zi) { 
      if (y == 0) { 
        return log_sum_exp(bernoulli_logit_lpmf(1 | zi), 
                           bernoulli_logit_lpmf(0 | zi) + 
                             poisson_log_lpmf(0 | eta)); 
      } else { 
        return bernoulli_logit_lpmf(0 | zi) +  
          poisson_log_lpmf(y | eta); 
      } 
    }
  // zero-inflated poisson log-CCDF and log-CDF functions
  real zero_inflated_poisson_lccdf(int y, real lambda, real zi) { 
    return bernoulli_lpmf(0 | zi) + poisson_lccdf(y | lambda); 
  }
  real zero_inflated_poisson_lcdf(int y, real lambda, real zi) { 
    return log1m_exp(zero_inflated_poisson_lccdf(y | lambda, zi));
  }
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  vector[N] cYear;  // population-level design matrix
  // data needed for ARMA correlations
  int<lower=0> Kar;  // AR order
  // number of lags per observation
  int<lower=0> J_lag[N];
  real meanResponse; //mean value of data
  real sdResponse; //sd value of data
}
transformed data {
  int max_lag = 1;
}
parameters {
  real b;  // trend
  real Intercept;  // temporary intercept for centered predictors
  vector[N] zerr;  // unscaled residuals
  real<lower=0> sderr;  // SD of residuals
  vector<lower=-1,upper=1>[Kar] ar;  // autoregressive coefficients
  real<lower=0,upper=1> zi;  // zero-inflation probability
}
transformed parameters {
  vector[N] err;  // actual residuals
  // compute ctime-series residuals
  err = sderr * zerr;
}
model {
    // matrix storing lagged residuals
    matrix[N, max_lag] Err = rep_matrix(0, N, max_lag);
    // initialize linear predictor term
    vector[N] mu = Intercept + cYear * b +  err;
    // include ARMA terms
    for (n in 1:N) {
      for (i in 1:J_lag[n]) {
        Err[n + 1, i] = err[n + 1 - i];
      }
      mu[n] += Err[n, 1:Kar] * ar;
    }
    for (n in 1:N) {
      target += zero_inflated_poisson_log_lpmf(Y[n] | mu[n], zi);
    }
  // priors including constants
  target += normal_lpdf(b | 0,5);
  target += student_t_lpdf(Intercept | 3, meanResponse, sdResponse);
  target += student_t_lpdf(sderr | 3, 0, sdResponse);
  target += std_normal_lpdf(zerr);
  target += beta_lpdf(zi | 1, 1);
}
